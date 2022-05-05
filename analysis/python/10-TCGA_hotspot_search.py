# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:36:53 2022

@author: ciara
"""

import numpy as np
import pandas as pd
import hdbscan
import hotmapper.mapper as hm
import hotmapper.hotspot as hmh
import hotmapper.utils as hmu
from itertools import product
from lifelines.statistics import logrank_test



#------------------------------read in files----------------------------------# 
#tcga gene expression after DSGA transformation 
X = "tcga/tcga_dct.csv" 

#mapper graph will be coloured by survival time censored to 10 years
surv = "tcga/10_year_survival.csv" #survival event & time censored to 10 years
y = surv["Time"].to_list() 

#lens function is built on the random index (r) and random value (nr) files generated
#by the hotspot search on the METABRIC dataset
r = "metabric/r.txt"
nr = "metabric/nr.txt"




#------------------------------functions--------------------------------------# 
def evaluate_hotspot_survival(mapper_hotspots):
    """Function evaluates the survival outcome of each hotspot using log-rank test. 
    Input is the hotspot dataframes of each component in the graph generated by the 
    hotmapper.hotspot.Subgraphs().run_hotspot_search() function"""
    
    #the number of components in the graphs 
    components = [c for c in mapper_hotspots]
    
    #for each component in the graph
    survival_results = []


    for i in components:
        #find the list of nodes for all hotspots 
        hotspot_samples = mapper_hotspots[i]["hotspot samples"]
        #for each hotspot    
        for k,samples in enumerate(hotspot_samples):     
            #seperate the er+ cohort into hotspot and global neighbourhood
            y = pd.DataFrame([0] * len(X.index), columns = ["Hotspot"], index = X.index)
            X_H = X.iloc[samples]
            y.loc[X_H.index,"Hotspot"] = 1
            
            #apply log rank test to each hotspot division 
            hot = (y["Hotspot"] == 1)
            T = surv["Time"]
            E = surv["Event"]
        
            #obtain the p-value 
            lr = logrank_test(T[hot], T[~hot], event_observed_A=E[hot], event_observed_B=E[~hot], alpha=99)
            survival_results.append(lr.p_value)
            print(f"\n pvalue : {lr.p_value}\n")
            
    #if multiple hotspots in a graph, prioritise the hotspot with the lowest p-value
    survival_min = np.min(survival_results)
    return survival_min
            


#----------------------search options --------------------------------#
#identify all combinations of interval and overlap options for search 
intervals =  range(10,32,2) #intervals randing from 10 to 32
overlap = np.linspace(0.1,0.45,8) #overlap ranging from 10% to 45%
io_list = list(product(intervals, overlap))



#----------------------build mapper graph--------------------------------#
#the mapper graph is built with the lens function identified from the hotspot search on METABRIC 
#to save time we only need to initialise the class and project the lens function once 
#initialise the mapper class
mapper = hm.Mapper(np.array(X), text =False)

#Project the dataset across the lens function identified from the search
lens_train = mapper.lens_function(selection = "linear_subset", r = r, nr = nr)




#----------------------search for hotspot-------------------------------#
#collected all successful interval and overlap options in a dictionary
io_success = {}

#for each grid combination of the interval & overlap, search lens for hotspot
for i in io_list:
    #select the interval and overlap options 
    i_param = i[0]
    o_param = i[1]
    print(f"intervals: {i_param}, overlap: {o_param}")
    
    # Build a cover on the lens function - specify the number of intervals and the percentage overlap
    mapper.covering(intervals = i_param, overlap = o_param)

    
    # Run hdbcsan clustering algorithm and build the graph 
    mapper.cluster_data(algorithm = hdbscan.HDBSCAN())
    
    #build the graph with edges and nodes, coloured by survival time
    graph = mapper.build_graph(attribute = np.array(s))


    #----------------------confirm presence of hotspots--------------------------------#
    #search for the presence of any hotspots with >6 months survival difference and minimum sample size of 30 
    graph_components = hmh.Subgraphs(mapper)
    hotspots = graph_components.run_hotspot_search(attribute_threshold = 6, min_sample_size = 30)
    
    
    #for each component in the graph, if a hotspot is present, test the survival outcome 
    for c in hotspots:
        if hotspots[c]['hotspot nodes']: 
            s_hot = evaluate_hotspot_survival(hotspots)
    
    io_success[(i_param,o_param)] = s_hot 


#find interval and overlap option with the strongest p-value 
optimal = min(io_success, key=io_success.get)
print(f"\n optimal set of parameters {optimal}")