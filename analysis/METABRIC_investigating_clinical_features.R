library(openxlsx)


#------------------------------------read in files---------------------------------------------#
y <- "metabric/hotspot_labels.csv" #binary labels indiciating hotspot class
clin <- "metabric/clinical.csv" #all clinical features for metabric cohort


#-------------------------------Chi Square Test--------------------------------------------#
#selct the correct patients 
clin_sub <- clin[rownames(clin) %in% rownames(y),]

#select the relevant clinical features 
#chemothreapy / hormone threapy / menopause / intclust / Claudin Subtype / Three Gene /Histological Subtype 
keep <- c("CHEMOTHERAPY", "INTCLUST", "HER2_SNP6" ,"INFERRED_MENOPAUSAL_STATE","CLAUDIN_SUBTYPE","THREEGENE","HISTOLOGICAL_SUBTYPE")
clin_sub <- clin_sub[,keep]

#divide clinical into hotspot and neighbourhood
clin_sub$HOTSPOT <- as.factor(y$Hotspot)
summary(clin_sub)

#build contigency tables 
cont_list <- apply(clin_sub, 2, function(x) split_categories(x))

#run chi-squared for every categorical column
chisq.stat <- apply(clin_sub, 2, function(x) chisq.test(x, clin_sub$HOTSPOT)$statistic)
chisq.pval <- apply(clin_sub, 2, function(x) chisq.test(x, clin_sub$HOTSPOT)$p.value)

#combine to single dataframe 
chisq.df <- cbind(chisq.stat,chisq.pval)
colnames(chisq.df) <- c("X.squared","P.val")

#add to list
cont_list$CHISQ <- chisq.df
cont_list

#write up excel file of lists 
write.xlsx(cont_list, file = "clinicalfeatures_chisq.xlsx",rowNames=T)


