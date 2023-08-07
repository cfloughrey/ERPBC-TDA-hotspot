# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:32:56 2022

@author: ciara
"""

import hotmapper as hm
import numpy as np
import pandas as pd
from sklearn import decomposition, manifold
import itertools
import hdbscan
from lifelines.statistics import logrank_test



#------------------------------read in files----------------------------------# 
X = "metabric/metabric_dct.csv" #metabric gene expression after DSGA transformation 
rfs = "metabric/10_year_rfs.csv" #relapse free event & time censored to 10 years. 

#define the attribute function as patients who have relapsed before 10 years
outcome = np.array((rfs['Time']<=120) & (rfs['Event'] == "1:Recurred")).astype(int)


#----------------------set up parameters----------------------------------# 
#initialise the search class from the hotmapper module
search = hm.automated_parameter_search.Search(np.array(X))

#select the parameter options for the search 
parameters = {"non_zero_lens_features" : int(X.shape[1]/2), #50\% of non-zero features in the lens fucntion
              "interval_list" : range(10,32,2) , #no. of interval options
              "overlap_list" : np.linspace(0.1,0.45,8), #percentage of overlap options
              "clustering_algorithm" : hdbscan.HDBSCAN(), #keep clustering algorithm consistent 
              "attribute_function" : outcome, #patients who have relapsed before 10 years
              "epsilon" : 0.1, #difference in attribute between hotspot and neighbourhood
              'min_samples' : 30, #minimum sample size for hotspot
              'extreme': "higher"} #hotspots with higher occurence of relapse


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
        
        #the sucessful parameters containing hotspots
        p_success = search.parameters

        #the attributes to build the lens function are in parameter_lens
        weights = search.parameter_lens["weights"]
        feature_list = search.parameter_lens["feature_list"]
        
        #-----------------------------run survival analysis----------------------------------# 
        survival_results = []
        hotspot_id = []
        columns = ["interval", "overlap", "nodes", "size", "logrank"]
        
        #multiple graphs are generated from the different successful parameter options 
        #the hotspots in each graph are tested for significant survival 
        for ps, collection in search.parameter_samples.items():
            for i, hotspot_samples in enumerate(collection):
                hotspot_results = []
                #seperate the er+ cohort into hotspot and global neighbourhood
                y = pd.DataFrame([0] * len(X.index), columns = ["Hotspot"], index = X.index)
                X_H = X.iloc[hotspot_samples]
                y.loc[X_H.index,"Hotspot"] = 1
                
                #apply log rank test to each hotspot division 
                hot = (y["Hotspot"] == 1)
                T = rfs["Time"]
                E = rfs["Event"]
        
                lr = logrank_test(T[hot], T[~hot], event_observed_A=E[hot], event_observed_B=E[~hot], alpha=99)
                pvalue = lr.p_value
                
                
                #if any hotspots have lower p-value than 0.001 then results are saved 
                if pvalue < 0.01:
                    #append results to list to build dataframe summarising survival analysis for each hotspot 
                    hotspot_id.append(str(ps[0]) + str(int(ps[1] * 100)) + str(i))
                    hotspot_results = [ps[0], ps[1], p_success[ps][i], sum(y["Hotspot"]), pvalue]
                    survival_results.append(hotspot_results)
                    signficance = True

                    
        #if any hotspots have lower p-value than 0.001 then search ends
        hotspot_df = pd.DataFrame(survival_results, columns = columns, index = hotspot_id)
        if hotspot_df.empty:
            count += 1
            print("Hotspot not significant ... continue search") 
            
            
    #------------------------save results to file-----------------------------------#
    else:
        hotspot_df.to_csv("metabric/hotspot_dataframe.csv")
    
        #save weights and order of features from randomly generated lens function to allow us to receate it later
        np.savetxt("metabric/weights.txt", weights ,delimiter=",")
        np.savetxt("metabric/feature_list.txt", feature_list ,delimiter=",")
         
        print("Hotspot significantly impacts survival \n search ends.") 

        break
    
