# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 09:38:33 2022

@author: ciara
"""

import hotmapper as hm
import numpy as np
import pandas as pd
from hdbscan import HDBSCAN
from itertools import chain


#------------------------------read in files----------------------------------# 
X = "metabric/metabric_dct.csv" #metabric gene expression after DSGA transformation 
rfs = "metabric/10_year_rfs.csv" #relapse free event & time censored to 10 years. 

#define the attribute function as patients who have relapsed before 10 years
outcome = np.array((rfs['Time']<=120) & (rfs['Event'] == "1:Recurred")).astype(int)

#lens function is built on the weights assigned to the subset list of features generated by the hotspot search
weights = "metabric/weights.txt"
feature_list = "metabric/feature_list.txt"

#----------------------build mapper graph--------------------------------#
#construct the lens function
linear_lens = hm.random_lens.Lens(np.array(X), nonzero_features = len(weights), weights = weights, feature_list = feature_list)

# build the mapper graph according to the identified parameters
mapper = hm.mapper.MapperGraph(data = np.array(X), 
                                lens_function = linear_lens["lens"], 
                                intervals = 24, 
                                overlap = 0.1,
                                clustering_algorithm = HDBSCAN())
mapper.build_graph() 

# visualise the graph, colouring by relapse outcome
hm.visualisation.draw_graph(mapper_graph = mapper.graph, 
          attribute_function = outcome, 
          samples_in_nodes = mapper.samples_in_nodes,
          size = 10, 
          style =3, 
          col_legend_title = "Relapse before\n 10 years",
          labels = True,
          tick_labels = False)


#----------------------confirm presence of hotspots--------------------------------#
hotspot_search = hm.hotspot.HotspotSearch(mapper_graph = mapper.graph,
                                          attribute_function = outcome, 
                                          samples_in_nodes = mapper.samples_in_nodes)

#check the distribution of hotspot nodes 
hotspot_nodes = hotspot_search.search_graph(attribute_threshold = 0.1, 
                                            min_sample_size = 30, 
                                            attribute_extreme = "higher", 
                                            plot_dendrogram = True )
n = list(chain.from_iterable(hotspot_nodes))

#visualise again but highlight the nodes in the graph identified as hotspots
hm.visualisation.draw_graph(mapper_graph = mapper.graph, 
          attribute_function = outcome, 
          samples_in_nodes = mapper.samples_in_nodes,
          size = 10, 
          style = 3, 
          hotspot_nodes = n,
          labels = False,
          tick_labels = False,
          col_legend_title = "Relapse before\n 10 years")
print(f"hotspots... {n}")



#save hotspot class 
node_id = mapper.samples_in_nodes
node_id.index = X.index
y = pd.DataFrame(node_id[n].max(axis=1), columns = ["Hotspot"])

#----------------------save files--------------------------------#
#save the 287 genes found in the lens function 
lens_genes = X.columns[feature_list].to_list()
pd.DataFrame(lens_genes).to_csv("metabric/lens_genes.csv")        

#save the hotspot class labels 
y.to_csv("metabric/hotspot_labels.csv")       

