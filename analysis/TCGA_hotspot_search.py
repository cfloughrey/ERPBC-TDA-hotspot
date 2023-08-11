# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:36:53 2022

@author: ciara
"""

import hotmapper as hm
import numpy as np
import pandas as pd
import hdbscan
from lifelines.statistics import logrank_test




#------------------------------read in files----------------------------------# 
#tcga gene expression after DSGA transformation 
X = "tcga/tcga_dct.csv" 

#mapper graph will be coloured by survival time censored to 10 years
surv = "tcga/10_year_survival.csv" #survival event & time censored to 10 years
outcome = np.array((surv['Time']<=120) & (surv['Event'] == "1:DECEASED")).astype(int)

#lens function is built on the random index (r) and random value (nr) files generated
#by the hotspot search on the METABRIC dataset
weights = "metabric/weights.txt"
feature_list = "metabric/feature_list.txt"




#----------------------set up parameters----------------------------------# 
#initialise the search class from the hotmapper module
search = hm.automated_parameter_search.Search(np.array(X))

#select the parameter options for the search 
discovery_lens = hm.random_lens.Lens(np.array(X), nonzero_features = len(weights), weights = weights, feature_list = feature_list) #set the predefined lens

parameters = {"predefined_lens" : discovery_lens, #we keep the lens function identified by the discovery search
              "interval_list" : range(10,32,2) , #no. of interval options
              "overlap_list" : np.linspace(0.1,0.45,8), #percentage of overlap options
              "clustering_algorithm" : hdbscan.HDBSCAN(), #keep clustering algorithm consistent 
              "attribute_function" : outcome, #patients who have relapsed before 10 years
              "epsilon" : 0.1, #difference in attribute between hotspot and neighbourhood
              'min_samples' : 30, #minimum sample size for hotspot
              'extreme': "higher"} #hotspots with higher occurence of relapse


#-----------------------------run search----------------------------------#
print("\nSearching lens space...")
print(F"{count} / {runs}")
    
#mapper graphs are built for each lens and tested for the presence of a hotspot
search.build_graphs(parameters)

#the sucessful parameters containing hotspots
p_success = search.parameters

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
        T = surv["Time"]
        E = surv["Event"]

        lr = logrank_test(T[hot], T[~hot], event_observed_A=E[hot], event_observed_B=E[~hot], alpha=99)
        pvalue = lr.p_value
        
        #append results to list to build dataframe summarising survival analysis for each hotspot 
        hotspot_id.append(str(ps[0]) + str(int(ps[1] * 100)) + str(i))
        hotspot_results = [ps[0], ps[1], p_success[ps][i], sum(y["Hotspot"]), pvalue]
        survival_results.append(hotspot_results)

            
#if any hotspots have lower p-value than 0.001 then search ends
hotspot_df = pd.DataFrame(survival_results, columns = columns, index = hotspot_id)
print(hotspot_df)

#find optimal interval and overlap option according to log-rank test 
optimal = hotspot_df.loc[hotspot_df['logrank'].idxmin()]
print(f"\n optimal set of parameters\n {optimal}")
