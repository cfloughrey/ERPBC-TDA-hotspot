#-----------------------------------read in files----------------------------------------#
#metabric cohort
X_meta <- "metabric/metabric_dct.csv" #DSGA transformed gene expression data of METABRIC cohort
y_meta <- "metabric/hotspot_labels.csv" #binary labels indiciating hotspot class
genes <- "kolgorov_genes.csv" #names of 112 differentially expressed genes 

#tcga cohort
X_tcga <- "tcga/tcga_dct.csv"#DSGA transformed gene expression data of TCGA cohort



#--------------------------------identify similar profiles---------------------------------#
#for each gene column compute the average value - this is the centroid for that feature 
h <- rownames(y_meta[y_meta$Hotspot == 1,]) #restrict metabric data to hotspot group


#the expression data subset to the hotspot & non-hotspot group by the KS genes
exp.ks.h <- X_meta[h, colnames(X_meta) %in% genes]

#left with a centroid co-ordinate of hotspot cluster 
centroid.h <- colMeans(exp.ks.h)

#subset tcga to 112 genes 
exp.t.ks <- tcga[, colnames(tcga) %in% genes]

#fpr each patient in the TCGA group 
#calculate the euclidean distance to the centroid 
dist.h <- apply(exp.t.ks, 1, function(x) dist(rbind(x, centroid.h), method = "canberra"))

#save file
write.csv(dist.h, "tcga/distance_to_centroid.csv")


