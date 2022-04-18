## Supplementary material for the article: 
# *Detecting subgroups of poor survival in estrogen receptor-positive breast cancer using topological data analysis*

Here we provide the code for the methodology of our paper "Detecting subgroups of poor survival in estrogen receptor-positive breast cancer using topological data analysis" (Loughrey, C.F., Maguire, S., Dlotko, P., Orr, N., Jurek-Loughrey, A.). For full details of the method see (link to published paper)

This repository also contains ```hot-mapper```, the general implementation of hotspot detection on the TDA tool Mapper for applications in bioinformatics. Mapper builds network graphs from a dataset through a combination of clustering and dimensionality reduction. Hotspot detection identifies anomalous interconnected groups of nodes in the graph. In biomedical terms, these hotspots can represent unusual patients subtypes. 


## Article supplementary material
### Data
The datasets used are found at the following links: 
- **METABRIC**:  raw gene expression data in IDAT format. European Genome-Phenome Archive, Dataset ID EGAD00010000162 (https://ega-archive.org/datasets/)
- **TCGA**:  raw fastq files. National Cancer Institute Genomic Data Commons Data Portal legacy archive, TCGA project access number 16762 (https://portal.gdc.cancer.gov/legacy-archive/search/f)
- **GTEX**:  raw gene read counts. Genotype-Tissue Expression Project, 'GTEx\_Analysis\_2017-06-05\_v8\_RNASeQCv1.1.9\_gene\_reads.gct.gz' (https://gtexportal.org/home/datasets)

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


## hot-mapper
The code to generate a Mapper graph and identify hotspots is available in the 'hotmapper' folder. For an example analysis see ```toy_example.ipynb```.
