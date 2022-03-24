#------------------------------------read in files---------------------------------------------#
X <- "metabric/metabric_dct.csv" #DSGA transformed gene expression data
y <- "metabric/hotspot_labels.csv" #binary labels indiciating hotspot class


#-------------------------------KOLGOROV-SMIRNOV--------------------------------------------#
#does a gene have the same continuous distribution for both hotspot and non-hotspot 
#split into hotspot and neighbour
h <- rownames(y[y$Hotspot == 1,,drop=FALSE])
n <- rownames(y[y$Hotspot == 0,,drop=FALSE])

#run kolgorov-smirnov test
ks.stat <- apply(X, 2, function(x) ks.test(x[h],x[n])$p.value)

#bonferonni correction
p.new <- round(p.adjust(ks.stat, "BH"), 3)
genes <- p.new[p.new < 0.01]
length(genes)
g.names <- names(genes)

#save results
write.csv(g.names, "output/kolgorov_genes.csv")