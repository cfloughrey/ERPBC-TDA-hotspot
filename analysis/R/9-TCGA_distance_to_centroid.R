#-----------------------------------read in files----------------------------------------#
#metabric cohort
X_meta <- "metabric/metabric_dct.csv" #DSGA transformed gene expression data of METABRIC cohort
y_meta <- "metabric/hotspot_labels.csv" #binary labels indiciating hotspot class
genes <- "kolgorov_genes.csv" #names of 112 differentially expressed genes 

#tcga cohort
X_tcga <- "tcga/tcga_dct.csv"#DSGA transformed gene expression data of TCGA cohort



#--------------------------------identify similar profiles---------------------------------#
#for each gene column compute the average value - this is the centroid for that feature 
#result is a centroid co-ordinate of the hotspot cluster 
h <- rownames(y_meta[y_meta$Hotspot == 1,]) #restrict metabric data to hotspot group
centroid <- apply(X_meta[h, colnames(X_meta) %in% genes],2,mean) #find the mean of all kologorov-sminorv genes

#fpr each patient in the TCGA group calculate the canberra distance to the centroid 
simlr <- apply(X_tcga[,colnames(X_tcga) %in% genes], 1, function(x) dist(rbind(x, centroid), method = "canberra"))

#reduce the tcga patients to those exceeded the 95th percentile 
#i.e. the closest patients to the centroid
limit <- quantile(simlr, 0.05)
tcga.group <- names(simlr[simlr<limit])

#create a label for the tcga hotspot group
X_tcga$Hotspot <- 0
X_tcga[tcga.group,"Hotspot"] <- 1

#save file 
lab <- X_tcga[,"Hotspot", drop = F]
write.csv(lab, "tcga/hotspot_labels.csv")






