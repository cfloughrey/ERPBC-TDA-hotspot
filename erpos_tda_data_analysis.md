# Detecting hotspots in ER+ breast cancer using TDA
The original datasets used in the analysis are publically available and described in the [main description](README.md). 

## Pre-process files for Disease Specific Genomic Analysis (DSGA)
The DSGA method uses a dataset of tumour samples and a dataset of normal tissue samples as input. Run [match_genes_in_datasets.py](ERPBC-TDA-hotspot/analysis/match_genes_in_datasets.py) to handle missing data and match the genes features used in the tumour and normal datasets. In our analysis we compile the input pair of tumour and normal datasets seperately for the METABRIC discovery cohort and the TCGA validation cohort. The GTEX dataset of gene expression from healthy breast tissue is used as the normal samples. 

The gene expression datasets are as follows:
+ METABRIC discovery dataset = 18930 genes for 1429 ER+ Breast Cancer (BC) patients
+ TCGA validation dataset = 19957 genes for 790 ER+ BC patients
+ GTEX normal dataset = 36043 genes for 169 healthy individuals

After preproccessing for DSGA, the tumour and normal dataset pairs are: 
+ Discovery pair = 17903 genes
+ Validation pair = 18406 genes

## Disease Specific Genomic Analysis (DSGA)
We run [DSGA.py](ERPBC-TDA-hotspot/analysis/DSGA.py) to build the disease component (DcT) of the tumour dataset by comparing it against the healthy state model representing as the dataset of normal tissue. In the discovery dataset we simultaneously perform feature reduction. In the validation DcT dataset we skip independent feature reduction and instead retain the gene features kept in the discovery DcT dataset. 

The output DcT gene expression datasets are: 
+ Discovery dataset (METABRIC) = 575 genes for 1429 ER+ Breast Cancer (BC) patients
+ Discovery dataset (TCGA) = 575 genes for 790 ER+ BC patients


## Searching for hotspots in the discovery dataset
Run [discovery_hotspot_search.py](ERPBC-TDA-hotspot/analysis/discovery_hotspot_search.py) to search for a hotspot of patients experiencing relapse before 10 years. We search for a hotspot in Mapper graphs generated from the discovery METABRIC DcT dataset across a predefined range of parameter options. If a hotspot fitting the conditions is found for any combination of Mapper parameters, we investigate the survival outcome for the group of patients contained in the hotspot against the rest of the cohort using log rank tests. If a hotspot group exists that has significantly higher occurence of relapse, we save the settings used to generate the lens function so results can be replicated.

## Build the Mapper graph - Discovery dataset
Run [discovery_mapper_graph.py](ERPBC-TDA-hotspot/analysis/discovery_mapper_graph) to look at the results of the hotspot search. Using the METABRIC discovery DcT dataset, we construct and visualise the Mapper graph and the hotspot using the parameters identified during the algorithm search. 

## Investigate the clinical features of the discovery dataset
Investigate associations between the hotspot group and clinical features (e.g. PAM50 subtype) in the discovery dataset using [discovery_investigating_clinical_features.R](ERPBC-TDA-hotspot/analysis/discovery_investigating_clinical_features.R)

## Perform survival analysis on the discovery hotspot group
[discovery_survival.R](ERPBC-TDA-hotspot/analysis/discovery_survival.R) plots the kaplan-meier curve comparing 10-year relapse-free survival outcome for the hotspot group against the neighbourhood in the METABRIC discovery dataset.  

## Differential expression analysis 
The scripts for differential expression analysis on the discovery and validation datasets are available on request. 

## Searching for hotspots in the validation dataset
In [validation_hotspot_search.R](ERPBC-TDA-hotspot/analysis/validation_discovery_survival.R) we search for hotspots in the validation set across the interval and overlap range, searching for patients experiencing event of death before 10 years. We keep the lens function constant, only considering the successful lens identified from the discovery dataset hotspot search. We compare the 10-year overall survival outcome for hotspots identified to select the final interval and overlap parameters. 

## Build the Mapper graph - Validation dataset
Run [validation_mapper_graph.py](ERPBC-TDA-hotspot/analysis/validation_mapper_graph) to look at the results of the hotspot search on the TCGA validation dataset. We visualise the overall survival outcome in the Mapper graph. We also compare the distance from the TCGA nodes to the METABRIC hotspot group centroid. 

## Perform survival analysis on the validation hotspot group
[validation_survival.R](ERPBC-TDA-hotspot/analysis/validation_survival.R) plots the kaplan-meier curve comparing 10-year survival outcome for the hotspot group against the neighbourhood in the TCGA validation dataset.  


