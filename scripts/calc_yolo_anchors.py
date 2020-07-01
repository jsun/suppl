import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandabox import imgUtils


def iou(box, clusters):
    x = np.minimum(clusters[:, 0], box[0])
    y = np.minimum(clusters[:, 1], box[1])
    intersection = x * y
    box_area = box[0] * box[1]
    cluster_area = clusters[:, 0] * clusters[:, 1]
    iou_ = intersection / (box_area + cluster_area - intersection)
    return iou_


def kmeans(boxes, k, dist=np.median,seed=2020):
    rows = boxes.shape[0]
    distances     = np.empty((rows, k)) ## N row x N cluster
    last_clusters = np.zeros((rows,))
    np.random.seed(seed)
    clusters = boxes[np.random.choice(rows, k, replace=False)]
    while True:
        for icluster in range(k): # I made change to lars76's code here to make the code faster
            distances[:,icluster] = 1 - iou(clusters[icluster], boxes)
        nearest_clusters = np.argmin(distances, axis=1)
        if (last_clusters == nearest_clusters).all():
            break
        for cluster in range(k):
            clusters[cluster] = dist(boxes[nearest_clusters == cluster], axis=0)
        last_clusters = nearest_clusters
    return clusters,nearest_clusters,distances


def plot_cluster_result(plt,clusters,nearest_clusters,WithinClusterSumDist,wh):
    current_palette = list(sns.xkcd_rgb.values())
    for icluster in np.unique(nearest_clusters):
        pick = nearest_clusters==icluster
        c = current_palette[icluster]
        plt.rc('font', size=8)
        plt.plot(wh[pick,0],wh[pick,1],'p',
                 color=c,
                 alpha=0.5,label='cluster = {}, N = {:6.0f}'.format(icluster,np.sum(pick)))
        plt.text(clusters[icluster,0],
                 clusters[icluster,1],
                 'c{}'.format(icluster),
                 fontsize=20,color='red')
        plt.title('Clusters')
        plt.xlabel('width')
        plt.ylabel('height')
    plt.legend(title='Mean IoU = {:5.4f}'.format(WithinClusterSumDist))




def calc_anchors(input_dpath, k_max):

    imgutils = imgUtils()
    wh = []


    for xml_fpath in glob.glob(input_dpath + '/*.xml'):
        xmlinfo = imgutils.parse_PascalVOC(xml_fpath)
        for obj in xmlinfo['objects']:
            w = (obj['xmax'] - obj['xmin']) / xmlinfo['shape'][0]
            h = (obj['ymax'] - obj['ymin']) / xmlinfo['shape'][1]
            wh.append([w, h])
    wh = np.array(wh)
    img_w, img_h, img_ch = xmlinfo['shape']

    results = {}
    for k in range(2, k_max):
        clusters, nearest_clusters, distances = kmeans(wh, k, seed=2020, dist=np.mean)
        WithinClusterMeanDist = np.mean(distances[np.arange(distances.shape[0]),nearest_clusters])
        result = {'clusters':             clusters,
                  'nearest_clusters':     nearest_clusters,
                  'distances':            distances,
                  'WithinClusterMeanDist': WithinClusterMeanDist}
        print('{:2.0f} clusters: mean IoU = {:5.4f}'.format(k,1-result['WithinClusterMeanDist']))
        results[k] = result


    figsize = (15,35)
    count = 1
    fig = plt.figure(figsize=figsize)
    for k in range(k_max - 2,k_max):
        result               = results[k]
        clusters             = result['clusters']
        nearest_clusters     = result['nearest_clusters']
        WithinClusterSumDist = result['WithinClusterMeanDist']

        ax = fig.add_subplot(k_max / 2, 2, count)
        plot_cluster_result(plt, clusters, nearest_clusters, 1 - WithinClusterSumDist, wh)
        count += 1
    #plt.show()

    anchors = results[k_max - 1]["clusters"] * np.array([img_w, img_h])
    print(anchors)


if __name__ == '__main__':

    input_dpath = sys.argv[1]
    calc_anchors(input_dpath, 20)


