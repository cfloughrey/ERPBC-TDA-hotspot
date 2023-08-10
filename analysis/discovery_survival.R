library(survminer)
library(survival)


#------------------------------------read in files---------------------------------------------#
rfs = "metabric/10_year_rfs.csv" #relapse free event & time censored to 10 years. 
y <- "metabric/hotspot_labels.csv" #binary labels indiciating hotspot class

#combine survival data and hotspot class for patients 
rfs$Hotspot <- y$Hotspot



#-------------------------------kaplan-meier--------------------------------------------#
#fit curves
fit <- survfit(Surv(Time, Event) ~ Hotspot, data = rfs)

#plot
p <- ggsurvplot(fit, 
                data = rfs, 
                pval = T,
                mark.time = T,
                xlab = "Months", 
                ylab = "10-Year Relapse Free Survival",
                censor = T,  
                break.time.by = 24, 
                risk.table = TRUE,
                conf.int = TRUE,# label curves directly
                legend.labs =  c("Neighbourhood","Hotspot")) # legend instead of direct label)
p

#summary of results
survdiff(Surv(Time, Event) ~ Hotspot, data = rfs)



