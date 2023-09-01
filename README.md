## Supplementary material for the article: 
# *A novel method for subgroup discovery in precision medicine based on topological data analysis*

Here we provide the code for the methodology of our paper "A novel method for subgroup discovery in precision medicine based on topological data analysis" (Loughrey, C.F., Maguire, S., Dlotko, P., Orr, N., Jurek-Loughrey, A.). 

This repository contains ```hot-mapper```, the general implementation of hotspot detection on the TDA tool Mapper for applications in bioinformatics. Mapper builds network graphs from a dataset through a combination of clustering and dimensionality reduction. Hotspot detection identifies anomalous interconnected groups of nodes in the graph. In biomedical terms, these hotspots can represent unusual patients subtypes. 


## Article supplementary material
### Data
The datasets used are found at the following links: 
- **METABRIC**:  raw gene expression data in IDAT format. European Genome-Phenome Archive, Dataset ID EGAD00010000162 (https://ega-archive.org/datasets/)
- **TCGA**:  raw fastq files. National Cancer Institute Genomic Data Commons Data Portal legacy archive, TCGA project access number 16762 (https://portal.gdc.cancer.gov/legacy-archive/search/f)
- **GTEX**:  raw gene read counts. Genotype-Tissue Expression Project, 'GTEx\_Analysis\_2017-06-05\_v8\_RNASeQCv1.1.9\_gene\_reads.gct.gz' (https://gtexportal.org/home/datasets)

### Code 
The code used in the article is available in the 'analysis' folder. The pipeline of analysis can be found in [Detecting hotspots in ER+ breast cancer using TDA](erpos_tda_data_analysis.md)

## hot-mapper
The code to generate a Mapper graph and identify hotspots is available in the 'hotmapper' folder. For an example analysis see [toy demonstration](toy_demonstration.ipynb).
