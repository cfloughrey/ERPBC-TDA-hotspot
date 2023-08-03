# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 13:01:17 2020

@author: ciara
"""

import numpy as np

def linear(data, n_features = 0, weights = None, feature_list = None):
    """Return a linear combination of features for each vector"""

    #define the number of vectors and features in the dataframe
    #if a subset isn't selected, keep all features in calculation
    total_samples, total_features = data.shape

    #create an array of random generated numbers
    if weights is None:
        weights = np.random.uniform(low=-1, high=1, size=total_features)


    # if the user specifies a subset of features
    if n_features >= total_features:
        raise ValueError(f"linear: n_features must be less than {total_features}")

    if n_features != 0:
        #random int value beween 1,n for size of subset features
        if feature_list is None:
            rng = np.random.default_rng()
            feature_list = rng.choice(total_features, n_features, replace=False)

        #create an array of random generated numbers for the subset features
        if weights is None:
            weights = np.random.uniform(low=-1, high=1, size=n_features)


    #empty lens function
    lens = np.zeros(total_samples)

    #for each vector in the dataset
    for v in range(total_samples):
        x = data[v,:]
        
        #for each feature in the vector,
        if n_features != 0: #subset features
            for k in range(n_features):
                feature_i = feature_list[k]
                #find the sum of all features x random
                lens[v] += x[feature_i] * weights[feature_i]

        else:    #all features
            for f in range(total_features):
                print(total_features)
                #find the sum of all features x random
                lens[v] += x[f] * weights[f]


    lens_settings = {"lens": lens,
                    "weights": weights,
                    "feature_list": feature_list}

    return lens_settings

def Lens(data, nonzero_features = 100, weights = None, feature_list = None, weight_range = [-1,1]):
    """Return a linear combination of a subset of features for each vector"""
    
    #define the number of features samples used in the data
    total_samples, total_features = data.shape
    
    #randomly select the subset of features and their corresponding weights
    if weights is None and feature_list is None:

        #create an index for the subset of features sampled in the lens function
        rng = np.random.default_rng()
        feature_list = rng.choice(total_features, nonzero_features, replace=False)
        
        #create a list of corresponding weights for each feature
        weights = np.random.uniform(low=weight_range[0], high=weight_range[1], size=len(feature_list))
    
    else:
        nonzero_features = len(feature_list)
        
    
    #for each sample in the dataset calculate the lens function 
    lens = np.zeros(total_samples)
    
    #for each vector in the dataset
    for v in range(total_samples):
        x = data[v,:]
    
        for k in range(nonzero_features):
            feature_i = feature_list[k]
            #find the sum of all features x random
            lens[v] += x[feature_i] * weights[k]
            
    lens_settings = {"lens": lens,
                    "weights": weights,
                    "feature_list": feature_list}
    
    return lens_settings
