import os
import sys
import glob
import shutil






if __name__ == '__main__':
    
    from_dpath = sys.argv[1]
    to_dpath = sys.argv[2]
    subset_id = int(sys.argv[3])
    
    for d in sorted(glob.glob(os.path.join(from_dpath, '*'))):
        print(os.path.basename(d))
        
        d_output_path = os.path.join(to_dpath, os.path.basename(d))
        
        if not os.path.exists(d_output_path):
            os.mkdir(d_output_path)
        
        for f in glob.glob(os.path.join(d, '*')):
            fid = int(os.path.splitext(os.path.basename(f))[0].split('_')[-1])
            
            if subset_id == 1:
                if fid < 101:
                    shutil.copy(f, d_output_path)
            else:
                if fid > 100:
                    shutil.copy(f, d_output_path)
            
            
    
    
    

