library(survminer)
library(survival)


#------------------------------------read in files---------------------------------------------#
surv <- "metabric/10_year_survival.csv" #survival event & time censored to 10 years"
y <- "metabric/hotspot_labels.csv" #binary labels indiciating hotspot class

#combine survival data and hotspot class for patients 
surv$Hotspot <- y$Hotspot



#-------------------------------kaplan-meier--------------------------------------------#
#fit curves
fit <- survfit(Surv(Time, Event) ~ Hotspot, data = surv)

#plot
p <- ggsurvplot(fit, 
                data = surv, 
                pval = T,
                mark.time = T,
                xlab = "Months", 
                ylab = "Overall survival probability",
                censor = T,  
                risk.table = TRUE,
                conf.int = TRUE,# label curves directly
                legend.labs =  c("Neighbourhood","Hotspot")) # legend instead of direct label)
p


#summary of results
survdiff(Surv(Time, Event) ~ Hotspot, data = surv)



