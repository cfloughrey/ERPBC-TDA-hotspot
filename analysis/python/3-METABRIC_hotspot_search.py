# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:32:56 2022

@author: ciara
"""

import hotmapper.parameter_search as hmps
import numpy as np
import pandas as pd
from sklearn import decomposition, manifold
import itertools
import hdbscan
from lifelines.statistics import logrank_test



#------------------------------read in files----------------------------------# 
X = "metabric_dct.csv" #metabric gene expression after DSGA transformation 
surv = "10_year_survival.csv" #survival event & time censored to 10 years

#mapper graph will be coloured by survival time censored to 10 years
time = surv["Time"].to_list() 



#----------------------set up parameters----------------------------------# 
#initialise the search class from the hotmapper module
search = hmps.Parameter_Search(np.array(X))

#select the parameter options for the search 
parameters = {"lens_option" : "linear_subset", #50% of features used in lens
              "interval_list" : range(10,32,2) , #intervals randing from 10 to 32
              "overlap_list" : np.linspace(0.1,0.45,8), #overlap ranging from 10% to 45%
              "clustering_algorithm" : hdbscan.HDBSCAN(), #hdbscan clustering algorithm
              "attribute_function" : time, #survival time colouring the graph
              "epsilon" : 12, #12 months difference between hotspot and neighbourhood
              'min_samples' : 30, #minimum 30 patients in a hotspot 
              'extreme': "lower"} #specify hotspots with low survival time




#-----------------------------run search----------------------------------#
#specify the number of times to run a search for a hotspot with signficant survival difference 
runs = 100 
signficance = False
count = 0

while count < runs: 
    if signficance == False:
        print("\nSearching lens space...")
        print(F"{count} / {runs}")
            
        #mapper graphs are built for each lens and tested for the presence of a hotspot
        search.build_graphs(parameters)
        
        #all information about hotspots contained in parameter_hs dictionary
        graph_info = search.parameter_hs
        
        #the sucessful parameters contained in parameter_io
        p_success = search.parameter_io

        #the attributes to build the lens function are in parameter_lens
        lens_info = search.parameter_lens

        
        #-----------------------------run survival analysis----------------------------------# 
        survival_results = []
        hotspot_id = []
        columns = ["graph", "component", "interval", "overlap", "size", "logrank"]
        
        #multiple graphs are generated from the different successful parameter options 
        #the hotspots in each graph are tested for significant survival 
        for i, graph in graph_info.items():
            #there can be multiple components in a graph
            for j, component in enumerate(graph.loc['hotspot samples']):
                #for each hotspot in a graph component
                for k, hotspot_samples in  enumerate(component):
                    
                    #seperate the er+ cohort into hotspot and global neighbourhood
                    y = pd.DataFrame([0] * len(X.index), columns = ["Hotspot"], index = X.index)
                    X_H = X.iloc[hotspot_samples]
                    y.loc[X_H.index,"Hotspot"] = 1
                    
                    #apply log rank test to each hotspot division 
                    hot = (y["Hotspot"] == 1)
                    T = surv["Time"]
                    E = surv["Event"]
            
                    lr = logrank_test(T[hot], T[~hot], event_observed_A=E[hot], event_observed_B=E[~hot], alpha=99)
                    pvalue = lr.p_value
                    
                    
                    #if any hotspots have lower p-value than 0.001 then results are saved 
                    if pvalue < 0.001:
                        #append results to list to build dataframe summarising survival analysis for each hotspot 
                        hotspot_id.append(str(i) + str(j) + str(k))
                        hotspot_results = [i, j, p_success[i][0], p_success[i][1], sum(y["Hotspot"]), pvalue]
                        survival_results.append(hotspot_results)
                        signficance = True
                    
        #if any hotspots have lower p-value than 0.001 then search ends
        hotspot_df = pd.DataFrame(survival_results, columns = columns, index = hotspot_id)
        if hotspot_df.empty:
            count += 1
            print("Hotspot not significant ... continue search") 
            
            
    #------------------------save results to file-----------------------------------#
    else:
        hotspot_df.to_csv("hotspot_dataframe.csv")
    
        #save linear random integers etc from filter function
        #to allow us to receate it later
        r = lens_info[0]
        nr = lens_info[1]
        np.savetxt("r.txt", r ,delimiter=",")
        np.savetxt("nr.txt", nr ,delimiter=",")
         
        print(F"Hotspot significantly impacts survival \n search ends.") 

        break
    
