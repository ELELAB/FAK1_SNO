#### survival analysis
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(plyr)
library(ggfortify)
library(ggplot2)
library(survminer)


source("code/survival_analysis_functions.r")


# TCGA cancer types
cancer_types <- c("BRCA",
                  "GBM",
                  "LUAD",
                  "UCEC",
                  "KIRC",
                  "HNSC",
                  "THCA",
                  "LUSC",
                  "PRAD",
                  "COAD",
                  "STAD",
                  "BLCA",
                  "LIHC",
                  "KIRP",
                  "ESCA",
                  "READ",
                  "KICH", 
                  "CHOL"
)

# HUGO gene names
geneName <- c("ADH5")


for(cancer in cancer_types){
  
  # Import filered mRNA expression file
  dataFilt <- get(load(paste0("../mRNA_expression_DEA/data/dataFilt/",cancer, "_dataFilt.rda")))
  dataFilt_mat <- as.matrix(dataFilt)
  
  # Sort out the gene of interest
  gene_pos <- which(rownames(dataFilt) %in% geneName)
  if(length(gene_pos)>0){
    
    # Import clinical TCGA data (produced with the TCGAbiolinks function 'GDCquery_clinic')
    clinical <- get(load(paste0("/data/user/shared_projects/gsnor_pan/TCGA/data/clinical/",cancer,"_clinical.rda")))
    
    # Do a survival analysis
    clinical_survival <- get_survival_table(dataFilt_mat,clinical,geneName, 0.25,0.75,"KM")
    
    # Plot survival analysis results
    survival_plot(clinical_survival,geneName,cancer,25)
    
    # Get p-values
    p.val <- p_func(clinical_survival)
  }
}

