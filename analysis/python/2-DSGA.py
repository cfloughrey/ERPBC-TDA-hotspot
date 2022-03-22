# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:52:14 2022

@author: ciara
"""


from hotmapper.DSGA_transformation import DSGA
import numpy as np
import pandas as pd



#------------------------------read in files----------------------------------# 
#the matches tumour and normal gene expression data for each er+ data collection
metabric_t = "metabric_tumour_matched.csv" #1429 patients, 17903 genes
metabric_n = "metabric_normal_matched.csv" #168 patients, 17903 genes

tcga_t = "tcga_tumour_matched.csv" #790 patients, 18406 genes
tcga_n = "tcga_normal_matched.csv" #168 patients, 18406 genes


#-----------run disease specific genomic analysis -----------------------------# 
#run DSGA to build the disease component transformation of the dataset
#threshold METABRIC DcT to those that show a signfiicant deviation
metabric_dct = DSGA(df_normal = metabric_n, df_tumour = metabric_t, threshold = True) #575 genes

#do not threshold TCGA DcT, but match to the METABRIC DcT genes 
tcga_dct = DSGA(df_normal = tcga_n, df_tumour =  tcga_t, threshold = False) #18406 genes
gene_list = metabric_dct.columns
tcga_dct_T = tcga_dct.T
tcga_dct_thres = tcga_dct_T[gene_list] #575 genes



#--------------------------------save files----------------------------------# 
pd.Series(gene_list).to_csv(F"dct_genelist.csv")
metabric_dct.to_csv(F"metabric_dct.csv")
tcga_dct_thres.to_csv(F"tcga_dct.csv")

