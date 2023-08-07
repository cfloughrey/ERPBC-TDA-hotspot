# Detecting hotspots in ER+ breast cancer using TDA
The original datasets used in the analysis are publically available and described in the [main description](README.md). To demonstate the order of the scripts used in the analysis we will use example data.  

## Pre-process files for Disease Specific Genomic Analysis (DSGA)
The DSGA method uses a dataset of tumour samples and a dataset of normal tissue samples as input. Run [1-match_genes_in_datasets.py](ERPBC-TDA-hotspot/analysis/1-match_genes_in_datasets.py) 


#18930 genes, 1429 er+ bc patients
tcga = "tcga/gene_expression.csv" #19957 genes, 790 er+ bc patients
gtex = "gtex/gene_expression.csv" #36043 genes, 169 healthy patients

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
