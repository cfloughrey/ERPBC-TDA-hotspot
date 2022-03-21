# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:00:39 2020

@author: ciara
"""

import networkx as nx
from sklearn.metrics import  silhouette_score
import numpy as np
from statistics import mean
from itertools import chain



def node_size(node_samples):
    """Find the number of patients in each node"""
    size_lst = [len(node_samples[i]) for i in node_samples]
    return size_lst


def range01(x):
    """Return the scaled transofmration of an array between 0 and 1"""
    r01 = np.array([(x- min(x))/ (max(x)- min(x))])
    return r01


def generate_artificial_hotspot(X, radius, random_seed = 42):
    #random sample from X
    np.random.seed(random_seed) 
    core_sample = X[np.random.randint(X.shape[0], size=1), :]

    #generate y labels
    y = [0] * X.shape[0]
    for i in X:
        dist = np.linalg.norm(i - core_sample)
        if dist < radius:
            l = np.where(X==i)
            index = l[0][1]
            y[index] = 1
    return y


def colour_by_y(node_index_dictionary, y_values, norm = False):
    """for each node, create a list of the y values for each index"""
    if norm == True:
        #normalise y values to range
        normal_y = range01(np.array(y_values))

        #then perform rest of algorithm
        node_values = []
        for node,l in node_index_dictionary.items():
            node_y = [normal_y[0][index] for index in l]
            node_values.append(np.mean(node_y))

    else:
        node_values = []
        for node,l in node_index_dictionary.items():
            node_y = [y_values[index] for index in l]
            node_values.append(np.mean(node_y))

    return node_values


def norm(power,mag,matrix):
    """Computing euclidean norm of the matrix"""
    m = np.abs(matrix)
    m = np.power(matrix,power)
    m = np.sum(m,axis=1)
    m = np.power(m,mag/power)
    return m


def flatten(listOfLists):
    "Convert a list of lists to a single list"
    return chain.from_iterable(listOfLists)


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def obtain_community_partition(G):
    """create subgraphs for each connected component in the graph"""
    S = {i:G.subgraph(c).copy() for i,c in enumerate(nx.connected_components(G))}
    return S


def build_neighbourhood(hotspot_dict):
    """Input a dictionary containing all the nodes in the component and
    the division of nodes into clusters. Output the global neighbourhood for
    each cluster"""
    neighbour = [list(np.setdiff1d(list(hotspot_dict["component"]),c)) for c in hotspot_dict["clusters"]]
    return neighbour


def identify_community_neighbours(community_nodes, original_graph):
    """Function identifies the neighbours for each hotspot community.
    Returns a list of the neighbour indexes"""
    neighbours = []
    #for each community
    for com in community_nodes:
        com_n =[]
        #for each node
        for node in com:
            #identify the neighbours in that node
            node_neighbour = [n for n in original_graph.neighbors(node)]
            com_n.append(list(node_neighbour))
        #create list of neighbours
        neighbours.append(list(com_n))
    #make a flat set of unique values for each component
    nghb_f = []
    for c in neighbours:
        flat = [item for sublist in c for item in sublist]
        flat_set = set(flat)
        nghb_f.append(flat_set)
    #obtain list of communities and their neighbours
    com_nghb = [[] for i in range(len(community_nodes))]

    #for each node in the list of neighbours
    for i, n in enumerate(nghb_f):
        #if that node is in one of the communities
        for j,com in enumerate(community_nodes):
            if any(item in n for item in com):
                com_nghb[i].append(j)
        #remove the hotspot itself
        if (i in com_nghb[i]):
            com_nghb[i].remove(i)

    return com_nghb



def community_size_samples(community_nodes, node_samples):
    """Function obtains the number of unique patients in each community
    partition in the Mapper graph. Node samples is a dictionary of sample
    indexes for each node"""
    com_ind = []
    #for each community
    for c in community_nodes:
        #for each node obtain the list of indexes
        com_ind_sum = [node_samples[n] for n in list(c)]
        #combine node index lists for each community
        flat_list = [item for sublist in com_ind_sum for item in sublist]
        com_ind.append(flat_list)
    #the size of unique indexes in each community
    ciu_size = [len(set(l)) for l in com_ind]
    return ciu_size




def community_filterval(community_nodes, node_filter_values):
    """Function obtains the mean filter value for each community in the Mapper graph"""
    com_mean = []
    #for each community
    for c in community_nodes:
        #for each node obtain the list of indexes
        com_values = [node_filter_values[i] for n,i in enumerate(c)]
        #combine node index lists for each community
        val_mean = np.mean(com_values)
        #append to list
        com_mean.append(round(val_mean,3))
    #the size of unique indexes in each community
    return com_mean



def cluster_filterval(cluster_nodes, node_filter_values):
    """Function obtains the mean filter value for each cluster in the Mapper graph"""
    #for each node obtain the list of indexes
    com_values = [node_filter_values[i] for n,i in enumerate(cluster_nodes)]
    #combine node index lists for each community
    val_mean = np.mean(com_values)
    #the size of unique indexes in each community
    return round(val_mean,3)
