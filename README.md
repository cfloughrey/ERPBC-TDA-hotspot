## Hotspot detection on Mapper graphs

This repository contains ```hot-mapper```, a python implementation of hotspot detection on the TDA tool Mapper for applications in bioinformatics. Mapper builds network graphs from a dataset through a combination of clustering and dimensionality reduction. Hotspot detection identifies anomalous interconnected groups of nodes in the graph. In biomedical terms, these hotspots can represent unusual patients subtypes within a larger disease cohort. 

Here we demonstrate the implementation of the methodology in our paper: "Detecting subgroups of poor survival in estrogen receptor-positive breast". 
For full details of the method see: (link to published paper)


## How to use? 

The code runs in two steps: 
1. Build a Mapper graph from the dataset 
2. Run hotspot detection on the graph 

For an example analysis see ```toy_examples.ipynb```.

