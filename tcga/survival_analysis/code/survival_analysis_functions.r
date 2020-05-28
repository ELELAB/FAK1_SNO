
# Get p-values
p_func <- function(clinical_survival){
  surv_object <- survfit(Surv(time,status)~group, data=clinical_survival)
  # Log-rang test
  sdf <- survdiff(Surv(time,status)~group, data = clinical_survival)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  return(p.val)
}

#function to get patient code
get_patient <- function(list_barcodes){
  IDs <- strsplit(c(as.character(list_barcodes)), "-")
  IDs <- ldply(IDs, rbind)
  IDs$patient <- apply(IDs[,c(1,2,3)],1,paste,collapse = "-" )
  return(IDs$patient)
}

# get survival time: if dead patient time=days_to_death, 
# if it is alive time:days_to_last_followup
get_survival_time <- function(clinical){
  
  for(i in 1:nrow(clinical)){
    if(clinical$vital_status[i]=="dead")
      clinical$time[i] <- clinical$days_to_death[i]
    else
      clinical$time[i] <- clinical$days_to_last_follow_up[i]
  }
  return(clinical)
}



#--------------------------------------------------------------------------------------

### Survival analysis
# method: "cox" or "KM" (Kaplan-Meier)
get_survival_table <- function(dataFilt, clinical, gene, threshDown, threshUp, method){
  
  #get only tumor samples
  dataSmTP <- TCGAquery_SampleTypes(colnames(dataFilt),"TP")
  
  #calculate quantile to be used as threshold to define two groups of samples to compare
  quantile_down <- quantile(dataFilt[gene,dataSmTP],threshDown)
  quantile_up <- quantile(dataFilt[gene,dataSmTP],threshUp)
  samples_down <- colnames(dataFilt[,dataSmTP])[which(dataFilt[gene,dataSmTP]< quantile_down)] 
  samples_up <- colnames(dataFilt[,dataSmTP])[which(dataFilt[gene,dataSmTP]> quantile_up)]
  patient_down <- get_patient(samples_down)
  patient_up <- get_patient(samples_up)
  clinical <- subset(clinical, clinical$bcr_patient_barcode %in% union(patient_up,patient_down),
                     select = c("bcr_patient_barcode","vital_status","days_to_death","days_to_last_follow_up","gender","tumor_stage","age_at_diagnosis","cigarettes_per_day"))
  
  #add survival time in clinical data
  clinical <- get_survival_time(clinical)
  
  # get time in years
  clinical$time <- clinical$time/365
  clinical$group <- as.character(lapply(clinical$bcr_patient_barcode,function(x) get_type_sample(x,patient_up,patient_down)))
  clinical$age_at_diagnosis <- (clinical$age_at_diagnosis)/365
  
  #include age-info and cigarettes-info only if we decide to perform cox regression
  if(method=="cox"){
    clinical <- subset(clinical,!is.na(clinical$age_at_diagnosis))
    median_age <- floor(median(clinical$age_at_diagnosis))
    clinical$age <- as.character(lapply(clinical$age_at_diagnosis,function(x) get_age_range(x,median_age)))
  }
  
  # get vital_status: dead=1, alive=0
  clinical$status <- as.numeric(lapply(clinical$vital_status, FUN=function(x) if(x=="dead") v<-1 else v<-0))
  clinical <- subset(clinical, !is.na(clinical$time))
  
  return(clinical)
}


## Plot survival analysis results
survival_plot <- function(clinical,gene,cancer_type,percentile){
  dir.create("plots", showWarnings = FALSE)
  surv_object <- survfit(Surv(time,status)~group, data=clinical)
  
  # Log-rang test
  sdf <- survdiff(Surv(time,status)~group, data = clinical)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  if(signif(p.val) <= 0.05){
    
    # Make survival plots
    myPlot <- autoplot(surv_object, conf.int = FALSE, ylab = "% surviving",xlab = "time (years)", surv.size = 1.5)+
      ggtitle(paste0(cancer, " GSNOR Logrank p-value=",signif(p.val)))+  
      theme(axis.title=element_text(size=18),axis.text=element_text(size=18),
            legend.text=element_text(size=18),legend.title=element_blank(),
            plot.title = element_text(hjust = 0.5,size=18), plot.subtitle = element_text(size=16),
            panel.background = element_rect(fill = "white",colour = "black",color = "black"))+
      scale_color_manual(values = c("red","blue"))
    ggsave(filename = paste0("plots/survival_plot_",cancer_type,".pdf"),plot = myPlot)
  }
  
}

