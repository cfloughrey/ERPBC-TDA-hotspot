## Supplementary material for the article: 
# *Detecting subgroups of poor survival in estrogen receptor-positive breast cancer using topological data analysis*

Here we provide the code for the methodology in our paper: "Detecting subgroups of poor survival in estrogen receptor-positive breast". 
For full details of the method see: (link to published paper)

This repository also contains ```hot-mapper```, the general implementation of hotspot detection on the TDA tool Mapper for applications in bioinformatics. Mapper builds network graphs from a dataset through a combination of clustering and dimensionality reduction. Hotspot detection identifies anomalous interconnected groups of nodes in the graph. In biomedical terms, these hotspots can represent unusual patients subtypes within a larger disease cohort. 


## Article supplementary material
The code used for the analysis in the article is available in the "analysis" folder. 

### datasets 

### folders 


## hot-mapper
The code to generate a Mapper graph and identify hotspots is available in the "hotmapper" folder. 

### How to use? 

The code runs in two steps: 
1. Build a Mapper graph from the dataset 
2. Run hotspot detection on the graph 

For an example analysis see ```toy_example.ipynb```.

