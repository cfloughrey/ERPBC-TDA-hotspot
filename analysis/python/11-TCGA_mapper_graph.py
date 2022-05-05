# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:25:15 2022

@author: ciara
"""



from ciaratools import experiments, data_collections
import numpy as np
import pandas as pd
import hdbscan
import matplotlib.pyplot as plt
import hotmapper.mapper as hm
import hotmapper.hotspot as hmh
import hotmapper.utils as hmu


# #------------------------------read in files----------------------------------# 
# X = "metabric/metabric_dct.csv" #metabric gene expression after DSGA transformation 

# #mapper graph will be coloured by survival time censored to 10 years
# surv = "metabric/10_year_survival.csv" #survival event & time censored to 10 years
# time = surv["Time"].to_list() 

# #lens function is built on the random index (r) and random value (nr) files generated
# #by the hotspot search
# r = "metabric/r.txt"
# nr = "metabric/nr.txt"



#-----------set up experiment------------------------------# 

project_directory = "C:/Users/ciara/OneDrive/Documents/PhD File Transfer/Projects/Project 1/HD June 21 ER+/"
experiment_name = "HD_pipeline/final"
dataset = "TCGA"

#creating new folder to hold results for experiment in project path 
output_path = experiments.results_folder(project_directory, experiment_name)



#-----------read in files----------------------------------# 
bc = data_collections.BreastCancerData(l = "ciara", d=dataset)   
X = bc.load_file("DcT_ER")


#for testing each hotspot for survival 
surv = bc.load_file("Survival_10")
s = surv["Time"].to_list()


labels = pd.read_csv(filepath_or_buffer = F"{output_path}/tcga_new_labels.csv",
                    delim_whitespace = False,
                    header = 0, index_col = 0)

y = labels["Hotspot"].to_list()

distance = pd.read_csv(filepath_or_buffer = F"{output_path}/distance.csv",
                    delim_whitespace = False, index_col = 0)

d = distance["x"].to_list()


###FILTER FUNCTION
#read in random value 
r_linsub = pd.read_csv(filepath_or_buffer = F"{output_path}/r.txt",
                   delim_whitespace = False,
                   header = None)

r = r_linsub[0].tolist()

nr_linsub = pd.read_csv(filepath_or_buffer = F"{output_path}/nr.txt",
                   delim_whitespace = False,
                   header = None,
                   dtype =np.int32)

nr = nr_linsub[0].tolist()



#----------------------build mapper graph--------------------------------#
#initialise the mapper class
mapper = hm.Mapper(np.array(X))

#Project the dataset across the lens function identified from the search
lens_train = mapper.lens_function(selection = "linear_subset", r = r, nr = nr)

#Build a cover on the lens function using the intervals and overlap with the highest p-value
#from the logrank test 
#28 0.1 = distance search 
mapper.covering(intervals = 28, overlap = 0.1)

# Run hdbcsan clustering algorithm and build the graph 
mapper.cluster_data(algorithm = hdbscan.HDBSCAN())

#build the graph with edges and nodes, coloured by survival time
graph = mapper.build_graph(attribute = d)

#visualise graph 
mapper.visualise(size = 10, style = 2, labels = False, col_legend_title = "distance to centroid")



#----------------------confirm presence of hotspots--------------------------------#
#confirm the presence of a hotspot in node[0]
graph_components = hmh.Subgraphs(mapper)
thres = (np.max(y)-np.min(y))/10
hotspots = graph_components.run_hotspot_search(attribute_threshold = thres, min_sample_size = 30)
print(hotspots)



#find out how many training true samples are in the hotspot class 
nlab = hotspots[0]["hotspot nodes"][0]
hnodes = [mapper.samples_in_node[i] for i in nlab]
hnodes_flat = np.unique(list(hmu.flatten(hnodes)))


surv["Label"] = 0
hotspot_patients = surv.iloc[hnodes_flat].index
surv.loc[hotspot_patients,"Label"] = 1
surv.to_csv(F"C:/Users/ciara/OneDrive/Documents/PhD File Transfer/Projects/Project 1/check_TCGA/results/hotspot_search/tcga_surv_label.csv")
