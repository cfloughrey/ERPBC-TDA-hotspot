# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 09:38:33 2022

@author: ciara
"""

import numpy as np
import pandas as pd
import hdbscan
import matplotlib.pyplot as plt
import hotmapper.mapper as hm
import hotmapper.hotspot as hmh



#------------------------------read in files----------------------------------# 
X = "metabric/metabric_dct.csv" #metabric gene expression after DSGA transformation 

#mapper graph will be coloured by survival time censored to 10 years
surv = "metabric/10_year_survival.csv" #survival event & time censored to 10 years
time = surv["Time"].to_list() 

#lens function is built on the random index (r) and random value (nr) files generated
#by the hotspot search
r = "metabric/r.txt"
nr = "metabric/nr.txt"



#----------------------build mapper graph--------------------------------#
#initialise the mapper class
mapper = hm.Mapper(np.array(X))

#Project the dataset across the lens function identified from the search
lens_train = mapper.lens_function(selection = "linear_subset", r = r, nr = nr)

#Build a cover on the lens function using the intervals and overlap with the highest p-value
#from the logrank test 
mapper.covering(intervals = 24, overlap = 0.1)

# Run hdbcsan clustering algorithm and build the graph 
mapper.cluster_data(algorithm = hdbscan.HDBSCAN())

#build the graph with edges and nodes, coloured by survival time
graph = mapper.build_graph(attribute = time)

#visualise graph 
mapper.visualise(size = 10, style = 2, labels = True, col_legend_title = "survival time \n(months)")



#----------------------confirm presence of hotspots--------------------------------#
#confirm the presence of a hotspot in node[2]
graph_components = hmh.Subgraphs(mapper)
hotspots = graph_components.run_hotspot_search(attribute_threshold = 12, min_sample_size = 30)
print(hotspots)

#save hotspot class 
y_samples = hotspots[0]["hotspot samples"][0]
y = pd.DataFrame([0] * len(X.index), columns = ["Hotspot"], index = X.index)
y.loc[X.index[y_samples],"Hotspot"] = 1



#----------------------save files--------------------------------#
#save the 287 genes found in the lens function 
lens_genes = X.columns[nr].to_list()
pd.DataFrame(lens_genes).to_csv("metabric/lens_genes.csv")        

#save the hotspot class labels 
y.to_csv("metabric/hotspot_labels.csv")       

