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
Run [hotspot_search.py](ERPBC-TDA-hotspot/analysis/hotspot_search.py) to search for a hotspot of patients experiencing relapse before 10 years. We search for a hotspot in Mapper graphs generated from the discovery METABRIC DcT dataset across a predefined range of parameter options. If a hotspot fitting the conditions is found for any combination of Mapper parameters, we investigate the survival outcome for the group of patients contained in the hotspot against the rest of the cohort using log rank tests. If a hotspot group exists that has significantly higher occurence of relapse, we save the settings used to generate the lens function so results can be replicated.

## Build the Mapper graph - Discovery dataset
Run [discovery_mapper_graph.py](ERPBC-TDA-hotspot/analysis/discovery_mapper_graph) to look at the results of the hotspot search. Using the METABRIC discovery DcT dataset, we construct and visualise the Mapper graph and the hotspot using the parameters identified during the algorithm search. 

## Investigate the clinical features of the discovery dataset
Investigate associations between the hotspot group and clinical features (e.g. PAM50 subtype) in the discovery dataset using [discovery_investigating_clinical_features.R](ERPBC-TDA-hotspot/analysis/discovery_investigating_clinical_features.R)

## Perform survival analysis on the discovery hotspot group
[discovery_survival.R](ERPBC-TDA-hotspot/analysis/discovery_survival.R) plots the kaplan-meier curve comparing 10-year relapse-free survival outcome for the hotspot group against the neighbourhood in the METABRIC discovery dataset.  

## Differential expression analysis 



## Perform survival analysis on the validation hotspot group
[validation_survival.R](ERPBC-TDA-hotspot/analysis/validation_survival.R) plots the kaplan-meier curve comparing 10-year survival outcome for the hotspot group against the neighbourhood in the TCGA validation dataset.  






### Code 
The code used for the analysis in the article is available in the 'analysis' folder. This is split between python and R scripts.
1. 'analysis/python/1-match_genes_in_datasets.py': create tumour and normal datasets with matching genes for metabric and tcga using gtex
2. 'analysis/python/2-DSGA.py': run disease specific genomic analysis (https://doi.org/10.1093/bioinformatics/btm033) to transform each breast cancer dataset and highlight the extent that diseased tissue deviates from healthy tissue
3. 'analysis/python/3-METABRIC_hotspot_search.py' : search for a lens that will reveal a hotspot of patients with significantly lower survival in the METABRIC dct dataset
4. 'analysis/python/4-METABRIC_mapper_graph.py' : recreate the Mapper graph identified from the hotspot search in METABRIC
5. 'analysis/R/5-METABRIC_investigating_clinical_features.R' : run chi-square tests on clinical features of METABRIC to test if hotspot differs from neighbourhood
6. 'analysis/R/6-METABRIC_investigate_differentiated_genes.R' : run Kolmogorov-Smirnov tests on METABRIC DcT data to identify which genes differentiate the hotspot group 
7. 'analysis/R/7-METABRIC_kaplan_meier_survival.R' : Plot kaplan-meier curve comparing survival outcome for hotspot against neighbourhood for METABRIC data 
8. 'analysis/python/8-plot-parameter-results.py' : Plot the log-rank p-value for each combination of Mapper interval and overlap parameters that produced a hotspot 
9. 'analysis/R/9-TCGA_distance_to_centroid.R' : Identify the group of tcga patients most similar to the metabric hotspot group centroid according to Kolgorov-Smirnov differentially expressed genes
10. 'analysis/python/10-TCGA_hotspot_search.py' : search for the interval and overlap parameters that reveal the hotspot of patients with significantly lower survival in the TCGA dct dataset, using the same lens function identified from the METABRIC dataset in step 3. 
11. 'analysis/python/11-TCGA_mapper_graph.py' : recreate the Mapper graph identified from the hotspot search in TCGA
12. 'analysis/R/12-TCGA_kaplan_meier_survival.R' : Plot kaplan-meier curve comparing survival outcome for hotspot against neighbourhood for TCGA data 
