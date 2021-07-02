import os
import sys
import argparse
import datetime
import numpy as np
import pandas as pd
import sklearn
import sklearn.svm
import sklearn.preprocessing
import sklearn.ensemble
import sklearn.tree
import sklearn.neighbors
import sklearn.linear_model
import sklearn.pipeline
import utils





def simulate(X, y, pipe, params, use_random_seed=True):
    output_matrix = None
    
    if use_random_seed:
        # same datasets will give the same random seed
        random_seed = int(np.sum(y))
    else:
        # random seed is depending on date and time and never be same even if the datasets are same
        random_seed = int(datetime.datetime.now().strftime('%Y%m%d%H%M%S')) + int(np.sum(y))
    
    # 10-fold cross validation for model validation
    i = 0
    kf = sklearn.model_selection.KFold(n_splits=10, random_state=random_seed, shuffle=True)
    for train_index, test_index in kf.split(X):
        i += 1
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        # 2-fold cross validation to determine hyper-paramaters of the model
        kf_ = sklearn.model_selection.KFold(n_splits=2, random_state=random_seed * 2, shuffle=True)
        #gs = sklearn.model_selection.GridSearchCV(pipe, params, n_jobs=-1, cv=kf_)
        gs = sklearn.model_selection.GridSearchCV(pipe, params, n_jobs=64, cv=kf_)
        gs.fit(X_train, y_train)
       
        pred_values = gs.predict(X_test)
        ouotput_matrix_ = pd.DataFrame({'true': y_test,
                                        'cv': i,
                                        'predicted': pred_values})
        
        output_matrix = pd.concat([output_matrix, ouotput_matrix_], axis=0)
        
    return output_matrix




def generate_pipeline(algorithm, feature_type, test_run):
    estimators = []
    pipe = None
    params = None
    
    # normalize values if not one-hot expression
    if feature_type == 'decimal':
        estimators.append(('scaler', sklearn.preprocessing.StandardScaler()))
    
    
    if algorithm == 'svm':
        estimators.append(('reg', sklearn.svm.SVR(kernel='rbf')))
        params = {
            'reg__gamma':   np.logspace(-5, 10, 30),
            'reg__C':       np.logspace(-5, 10, 30),
            'reg__epsilon': np.logspace(-5, 10, 30),
        }
        if test_run:
            params = {
                'reg__gamma':   np.array([0.001, 0.1, 10]),
                'reg__C':       np.array([0.001, 0.1, 10]),
                'reg__epsilon': np.array([0.001, 0.1, 10]),
            }
        
    elif algorithm == 'rf':
        estimators.append(('reg', sklearn.ensemble.RandomForestRegressor(criterion='mse')))
        params = {
            'reg__n_estimators': np.arange(2, 50, 1),
            'reg__max_depth':    np.arange(2, 10, 1),
        }
        if test_run:
            params = {
                'reg__n_estimators': np.array([2, 4, 6]),
                'reg__max_depth':    np.array([2, 3, 4]),
            }
        
    elif algorithm == 'dc':
        estimators.append(('reg', sklearn.tree.DecisionTreeRegressor(criterion='mse')))
        params = {
            'reg__max_depth': np.arange(2, 10, 1),
        }
        if test_run:
            params = {
                'reg__max_depth': np.array([2, 3, 4]),
            }
    
    elif algorithm == 'knn':
        estimators.append(('reg', sklearn.neighbors.KNeighborsRegressor()))
        params = {
            'reg__n_neighbors': np.arange(2, 10, 1),
        }
        if test_run:
            params = {
                'reg__n_neighbors': np.array([2, 3, 4]),
            }
    
    elif algorithm == 'lasso':
        estimators.append(('reg', sklearn.linear_model.Lasso()))
        params = {
            'reg__alpha': np.linspace(0, 1, 20),
        }
        if test_run:
            params = {
                'reg__alpha': np.array([0, 0.5, 1.0]),
            }
        
    elif algorithm == 'elasticnet':
        estimators.append(('reg', sklearn.linear_model.ElasticNet()))
        params = {
            'reg__alpha': np.linspace(0, 1, 20),
            'reg__l1_ratio': np.linspace(0, 1, 20),
        }
        if test_run:
            params = {
                'reg__alpha': np.array([0, 0.5, 1.0]),
                'reg__l1_ratio': np.array([0, 0.5, 1.0]),
            }
        
    else:
        raise ValueError('set the correct algorithm!')

    pipe = sklearn.pipeline.Pipeline(estimators)
    
    
    return pipe, params






def main(algorithm, dataset, feature_type, randomize_type, output, test_run):
    
    pipe, params = generate_pipeline(algorithm, feature_type, test_run)
    data = utils.load_dataset(dataset)
    data = utils.randomize_dataset(data, randomize_type)
    X, y = utils.get_xy(data, feature_type)
    pred_values = simulate(X, y, pipe, params, use_random_seed=False)
    pred_values.to_csv(output, header=True, index=False, sep='\t')




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Model Generator')
    parser.add_argument('--algorithm', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--feature-type', required=True)
    parser.add_argument('--randomize-type', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--test-run', action='store_true')

    args = parser.parse_args()

    main(args.algorithm, args.dataset, args.feature_type, args.randomize_type, args.output, args.test_run)
    
    
    





