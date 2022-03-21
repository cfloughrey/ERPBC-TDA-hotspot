# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:45:19 2022

@author: ciara
"""

import numpy as np 
import pandas as pd
from sklearn.impute import KNNImputer


def check_for_missing_data(list_of_datasets):
    sum_missing = [data[data.isnull().any(axis=1)].shape[0] for data in list_of_datasets]     
    return sum_missing 


def match_tumour_normal_data(list_of_tumour_datasets, normal_dataset):
    for dataset_name,tumour_dataset in list_of_tumour_datasets.items():
        #transpose data
        tumour_dataset_T = tumour_dataset.T
        normal_dataset_T = normal_dataset.T
        
        #obtain list of genes in both datasets 
        tumour_genes = tumour_dataset_T.columns
        normal_genes = normal_dataset_T.columns
    
        #match the names 
        matched_genes = [i for i in normal_genes if i in tumour_genes]
        
        #subset original datasets
        tumour_matched = tumour_dataset_T[matched_genes]
        normal_matched = normal_dataset_T[matched_genes]
        
        print(F"{dataset_name}: \n tumour shape ({tumour_matched.shape} \n normal shape {normal_matched.shape})")
        
        #save to file 
        tumour_matched.to_csv(F"{dataset_name}_tumour_matched.csv") 
        normal_matched.to_csv(F"{dataset_name}_normal_matched.csv") 





#-----------read in files----------------------------------# 
#read in the gene expression datasets 
metabric = "metabric/gene_expression.csv" #18930 genes, 1429 er+ bc patients
tcga = "tcga/gene_expression.csv" #19957 genes, 790 er+ bc patients
gtex = "gtex/gene_expression.csv" #36043 genes, 169 healthy patients



#-----------Impute missing data----------------------------------# 
#check for missing data in the datasets         
missing_data_no = check_for_missing_data([metabric, tcga, gtex])
print(F"number of missing genes: {missing_data_no}")

#There is missing data for 434 genes in the TCGA dataset, these are imputed using KNN
imputer = KNNImputer(n_neighbors=10)
tcga = pd.DataFrame(imputer.fit_transform(tcga), columns = tcga.columns, index = tcga.index)

#There is a large amount of missing data in GTEX (2594 genes). These genes are removed 
gtex = gtex.replace(-np.inf, np.nan).dropna()



#-----------Match genes between tumour and normal data-------------------------# 
#the matching genes between the tumour and normal data are identified and saved to new files 
#this has to be run seperately for metabric and tcga 

#create a dictionary of only the er+ bc tumour datasets
bc = {"metabric" : metabric, "tcga" : tcga}

#save tumour and normal data for each breast cancer dataset to new dataframe
#metabric: 17903 genes, tcga: 18406 genes 
match_tumour_normal_data(bc, gtex) 
