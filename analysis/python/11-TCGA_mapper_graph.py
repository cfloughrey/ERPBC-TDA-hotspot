# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:39:42 2022

@author: ciara
"""



from ciaratools import experiments, data_collections
import numpy as np
import pandas as pd
import hdbscan
import hotmapper.mapper as hm
import hotmapper.hotspot as hmh
import hotmapper.utils as hmu



#------------------------------read in files----------------------------------# 
#tcga gene expression after DSGA transformation 
X = "tcga/tcga_dct.csv" 

#mapper graph will be coloured by survival time censored to 10 years
surv = "tcga/10_year_survival.csv" #survival event & time censored to 10 years
y = surv["Time"].to_list() 

#list describing the distance of each tcga patient to the metabric centroid 
d = "tcga/distance.csv"

#lens function is built on the random index (r) and random value (nr) files generated
#by the hotspot search on the METABRIC dataset
r = "metabric/r.txt"
nr = "metabric/nr.txt"



#----------------------build mapper graph--------------------------------#
#initialise the mapper class
mapper = hm.Mapper(np.array(X))

#Project the dataset across the lens function identified from the search
lens_train = mapper.lens_function(selection = "linear_subset", r = r, nr = nr)

#Build a cover on the lens function using the intervals and overlap identified from the search on TCGA
mapper.covering(intervals = 30, overlap = 0.3)

# Run hdbcsan clustering algorithm and build the graph 
mapper.cluster_data(algorithm = hdbscan.HDBSCAN())

#build the graph with edges and nodes, coloured by survival time
graph = mapper.build_graph(attribute = y)

#visualise graph by survival 
mapper.visualise(size = 10, style = 2, labels = False, col_legend_title = "survival time (months)")



#----------------------confirm presence of hotspots--------------------------------#
#confirm the presence of a hotspot in node 1
graph_components = hmh.Subgraphs(mapper)
hotspots = graph_components.run_hotspot_search(attribute_threshold = 6, min_sample_size = 30)
print(hotspots)


#hotspot is in node 1, which is the first hotspot group of the first component
#obtain all the samples    
nlab = hotspots[0]["hotspot nodes"][0]
hnodes = [mapper.samples_in_node[i] for i in nlab]
hnodes_flat = np.unique(list(hmu.flatten(hnodes)))

#create new dataframe with hotspot group labels for tcga 
h_tcga = pd.DataFrame(index=surv.index, columns = ["Hotspot"])
h_tcga["Hotspot"] = 0
hotspot_patients = h_tcga.iloc[hnodes_flat].index
h_tcga.loc[hotspot_patients,"Hotspot"] = 1

#save tcga hotspot group labels
h_tcga.to_csv("tcga/hotspot_labels.csv")



#-------------------investigate hotspot distance to metabric centroid----------------#
#we can also visualise graph by the distance to metabric centroid 
graph = mapper.build_graph(attribute = d)
#visualise graph by survival 
mapper.visualise(size = 10, style = 2, labels = False, col_legend_title = "distance to centroid")
