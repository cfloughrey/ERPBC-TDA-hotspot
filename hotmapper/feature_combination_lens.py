# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 13:01:17 2020

@author: ciara
"""

import numpy as np
from itertools import combinations_with_replacement


def linear(X, r = None):
    """Return a linear combination of features for each vector"""
    #define the number of vectors and features in the dataframe
    i, n = X.shape

    #create an array of random generated numbers
    if r is None:
        r = np.random.uniform(low=-1, high=1, size=n)

    #filter function
    ff = np.zeros(i)

    #for each vector in the dataset
    for v in range(i):
        x = X[v,:]
        #for each feature in the vector
        for w in range(n):
            #find the sum of all features x random
            ff[v] += x[w] * r[w]

    return ff, r


def non_linear(X, r = None):
    """Return a non-linear combination of features for each vector"""
    #define the number of vectors and features in the dataframe
    i, n = X.shape

    #(r is outside the for loop as it needs to be same for every vector)
    #create a vector of the number of every possible feature pair combination
    pair_n = len(list(combinations_with_replacement(range(n), 2)) )

    #create a array of random numbers of same length
    if r is None:
        r = np.random.uniform(low=-1, high=1, size=pair_n)

    #filter function
    ff = np.zeros(i)

    #for each vector
    for v in range(i):
        x = X[v,:]

        #return a list of all possible combinations
        pair = list(combinations_with_replacement(x, 2))

        #for each combination of features
        for w in range(len(pair)):
            #find the sum of each pair x random
            ff[v] += pair[w][0] * pair[w][1] * r[w]

    return ff, r


def linear_subset(X, g, r = None, nr = None):
    """Create a linear combination of a subset of features, specifiying
    g between range (0,1) """
    #define the number of vectors and features in the dataframe
    i, n = X.shape

    #define subset value m
    m = int(g * n)

    #random int value beween 1,n for size of subset features
    if nr is None:
        rng = np.random.default_rng()
        nr = rng.choice(n, m, replace=False)

    #create a array of random numbers of same length
    if r is None:
        r = np.random.uniform(low= -1, high=1, size=n)


    #filter function
    ff = np.zeros(i)

    #for each vector
    for v in range(i):
        x = X[v,:]
        #for each value in the range(m) (new subset of features)
        for k in range(m):
            #index the list that specifies which features to select
            nr_i = nr[k]
            #find linear combination of selected features
            ff[v] += x[nr_i] * r[nr_i]

    return ff, r, nr



def manual_linear_subset(X, g, nr, r):
   """Manually recreate a linear combination of a subset of features, limiting
   g between range (0,1). The random interger value (nr) and random numbers (r)
   should be specified"""

   #define the number of vectors and features in the dataframe
   i, n = X.shape

   #define subset value m
   m = int(g * n)


   #filter function
   ff = np.zeros(i)

   #for each vector
   for v in range(i):
       x = X[v,:]
       #for each value in the range(m) (new subset of features)
       for k in range(m):

           #index the list that specifies which features to select
           nr_i = nr[k]
           #find linear combination of selected features
           ff[v] += x[nr_i] * r[nr_i]

   return ff


def non_linear_subset(X,g, r = None, nr = None):
    """Create a non linear combination of a subset of features, specifiying
    g as the number of features to use"""
    #define the number of vectors and features in the dataframe
    i, n = X.shape

    #create a array of random numbers of length p
    if r is None:
        r = np.random.uniform(low=-1, high=1, size= g)

    #for the selected number of pairs
    if nr is None:
        nr = [(np.random.randint(0, n, 2)) for i in range(g)]

    #filter function
    ff = np.zeros(i)

    #for each vector in the dataset
    for v in range(i):
        x = X[v,:]

        #for the number of feature pairs
        for k in range(g):
            f1 = nr[k][0]
            f2 = nr[k][1]

            #find linear combination of selected features
            ff[v] += x[f1] * x[f2] * r[k]

    return ff, r, nr
