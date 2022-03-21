## Supplementary material for the article: 
# *Detecting subgroups of poor survival in estrogen receptor-positive breast cancer using topological data analysis*

Here we provide the code for the methodology of our paper "Detecting subgroups of poor survival in estrogen receptor-positive breast cancer using topological data analysis" (Loughrey, C.F., Maguire, S., Dlotko, P., Orr, N., Jurek-Loughrey, A.). For full details of the method see (link to published paper)

This repository also contains ```hot-mapper```, the general implementation of hotspot detection on the TDA tool Mapper for applications in bioinformatics. Mapper builds network graphs from a dataset through a combination of clustering and dimensionality reduction. Hotspot detection identifies anomalous interconnected groups of nodes in the graph. In biomedical terms, these hotspots can represent unusual patients subtypes. 


## Article supplementary material
The code used for the analysis in the article is available in the 'analysis' folder. 

### Data
The datasets used are found at the following links: 
- **METABRIC**:  raw gene expression data in IDAT format. European Genome-Phenome Archive, Dataset ID EGAD00010000162 (https://ega-archive.org/datasets/)
- **TCGA**:  raw fastq files. National Cancer Institute Genomic Data Commons Data Portal legacy archive, TCGA project access number 16762 (https://portal.gdc.cancer.gov/legacy-archive/search/f)
- **GTEX**:  raw gene read counts. Genotype-Tissue Expression Project, 'GTEx\_Analysis\_2017-06-05\_v8\_RNASeQCv1.1.9\_gene\_reads.gct.gz' (https://gtexportal.org/home/datasets)

### Code 


## hot-mapper
The code to generate a Mapper graph and identify hotspots is available in the 'hotmapper' folder. For an example analysis see ```toy_example.ipynb```.
