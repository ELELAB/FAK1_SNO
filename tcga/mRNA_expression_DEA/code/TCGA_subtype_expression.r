library(SummarizedExperiment)
library(TCGAbiolinks)
#library(GenomicDataCommons)
library(biomaRt)
#library(ggplot2)
library(reshape2)
library(ggpubr)

# Source the revised DEA function
source("/data/tools_scripts_repository/TCGA_DEAfunction_06042020/TCGAanalyze_DEA_revised.R")


# List of cancer types. For example: c("BRCA","GBM","LUAD", ...)
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
                  "KICH"#,
                  #"CHOL"
)

# List of HUGO gene names to sort out for plots
geneName <- c("ADH5")

# Choose tumor sample from project (long name). For example: "Primary solid Tumor"
tumor_sample <- "Primary solid Tumor"

# Short name for tumor sample (must match chosen tumor sample). 
# For example (if tumor sample is "Primary solid Tumor"): "TP" 
tumor_short <- "TP"

# Choose normal sample from project (long name). For example: "Solid Tissue Normal
normal_sample <- "Solid Tissue Normal"

# Short name for normal sample (must match chosen normal sample). 
# For example (if normal sample is "Solid Tissue Normal"): "NT" 
normal_short <- "NT"

# Create an empty data frame to fill with information for every cancer type in the for-loop below
dataDEGs_cancer <- data.frame()
for(cancer in cancer_types){
  
  # Collect information about TCGA cancer subtypes
  tabPanCancer <- as.data.frame(PanCancerAtlas_subtypes())
  
  # Sort out that cancer type from TCGA subtype table
  cancer_subs <- tabPanCancer[which(tabPanCancer$cancer.type == cancer),]
  
  # Get which subtypes belong to that cancer type
  subtypes <- unique(cancer_subs$Subtype_Selected)
  subtypes <- sort(subtypes)
  
  # Subtypes to not include in the analysis
  donts <- c("NotAssigned", "Normal", paste0(cancer,"_LGG.NA"), paste0(cancer, ".NA"))
  donts_len <- which(subtypes %in% donts)
  if(length(donts_len) > 0){
    subtypes <- subtypes[- which(subtypes %in% donts)]
  }
  
  # Create empthy data frames to fill with information for every subtype, one for expression data and one for DEA
  all_dataFilt_sub <- data.frame()
  dataDEGs_all <- data.frame()
  
  # Do the expression and DEA analyses subtype by subtype
  for(subtype in subtypes){
    
    # Collect barcodes/samples for that subtype
    barcodes <- cancer_subs[which(cancer_subs$Subtype_Selected == subtype),]$pan.samplesID
    if(length(barcodes) > 9){
      
      # Collect full barcode names
      project <- paste0("TCGA-", cancer)
      sub.query.exp <- GDCquery(project = project,
                                  data.category = "Transcriptome Profiling",
                                  data.type = "Gene Expression Quantification",
                                  workflow.type = "HTSeq - Counts",
                                  barcode = barcodes,
                                  sample.type = c(tumor_sample),
                                  legacy = FALSE)
      barcodes_long <- getResults(sub.query.exp,cols="cases")
      
      TP <- TCGAquery_SampleTypes(barcodes_long, "TP")
      if(length(TP) > 4){
        
        # Import filtered mRNA expression data file (produced with TCGAbiolinks function 'TCGAanalyze_Filtering')
        dataFi <- get(load(paste0("data/dataFilt/",cancer,"_dataFilt.rda")))
        
        # Get barcodes that belong to normal samples
        NT <- TCGAquery_SampleTypes(colnames(dataFi), "NT")
        
        # Perform a limma + voom transformation
        dataFilt_voom <- limma::voom(dataFi)
        dataFilt <- as.data.frame(dataFilt_voom[[1]])
        
        # Sort out the barcodes that belong to that subtype + all normal samples
        dataFilt_sub <- dataFilt[,which(colnames(dataFilt) %in% c(TP, NT))]
        
        # Sort out the gene of interest
        dataFilt_s_gene <- as.data.frame(t(dataFilt_sub[which(rownames(dataFilt_sub) %in% geneName),]))
        
        # Do not continue if the gene/genes were not found in the filtered data for that subtype
        if(nrow(dataFilt_s_gene) > 0){
          
          # Add subtype and sample type information to 'dataFilt_s_gene', row by row
          dataFilt_s_gene$subtype <- rep("NA", nrow(dataFilt_s_gene))
          dataFilt_s_gene$sample_type <- rep("NA", nrow(dataFilt_s_gene))
          for(rownumb in 1:nrow(dataFilt_s_gene)){
            
            # If that row belong to a normal sample ...
            if(rownames(dataFilt_s_gene)[rownumb] %in% NT){
              
              # Get subtype : Normal tissue, to make sure they are separated from the tumor samples
              # When we later plot these results, we want normal samples to be the first box when we sort the subtypes alphabetically
              dataFilt_s_gene[rownumb,]$subtype <- "AANormal_Tissue"
              dataFilt_s_gene[rownumb,]$sample_type <- "Normal"
            }
            # If that row belong to a tumor sample ...
            if(rownames(dataFilt_s_gene)[rownumb] %in% TP){
              
              # Add subtype name and tumor sample type
              dataFilt_s_gene[rownumb,]$subtype <- subtype
              dataFilt_s_gene[rownumb,]$sample_type <- "Tumor"
            }
          }
          
          # Add information about mRNA expression/the samples for this subtype to the data frame created above
          all_dataFilt_sub <- rbind(all_dataFilt_sub, dataFilt_s_gene)
          
          # Continue to perform a DEA for this subtype
          dataFilt_DEA <- dataFi[,which(colnames(dataFi) %in% c(TP, NT))]
          
          # Divide filtered mRNA expression data file into two files:
          #  one for only tumor samples (for that subtype) and one for only normal samples
          dataT <- dataFilt_DEA[,which(colnames(dataFilt_DEA) %in% TP)]
          dataN <- dataFilt_DEA[,which(colnames(dataFilt_DEA) %in% NT)]
          
          # Use the revised DEA function to perform a DEA
          dataDEGs <- TCGAanalyze_DEA_revised(mat1 = dataN,
                                      mat2 = dataT,
                                      pipeline = "limma", 
                                      batch.factors = "TSS",
                                      Cond1type = "Normal",
                                      Cond2type = "Tumor",
                                      fdr.cut = 0.05 ,
                                      logFC.cut = 0.5,
                                      voom = TRUE)
          save(dataDEGs, file = paste0("data/DEA/",cancer,"_",subtype,"_dataDEGs.rda"))
          
          # Sort out the gene/genes of interest
          dataDEGs_gene <- dataDEGs[which(rownames(dataDEGs) %in% geneName),]
          
          # If the gene/genes were found in the dataDEGs file ...
          if(nrow(dataDEGs_gene) > 0){
            
            # Add subtype information and add it to the data frame created above
            dataDEGs_gene$subtype <- rep(subtype, nrow(dataDEGs_gene))
            dataDEGs_all <- rbind(dataDEGs_all, dataDEGs_gene)
          } 
          
        }
      }
    }
  }
  
  ### Prepare files for plotting
  
  # Remove duplicated rows
  # The data frame below now contain expression values for every subtype for this cancer type
  # This is done in order to plot all subtypes in the same plot
  all_dataFilt_sub <- all_dataFilt_sub[which(duplicated(all_dataFilt_sub) == FALSE),]
  if(length(all_dataFilt_sub) > 0 ){
    
    # Melt data frame to make it fit for ggplot boxplots
    dataFilt_sub_melt <- melt(all_dataFilt_sub)
    
    # Order subtypes alphabetically (and here normal samples will come first)
    dataFilt_sub_melt <- dataFilt_sub_melt[order(dataFilt_sub_melt$subtype),]
    
    # Change the normal samples to "Normal tissue" instead (looks better in the plot)
    dataFilt_sub_melt[which(dataFilt_sub_melt$subtype == "AANormal_Tissue"),]$subtype <- rep("Normal Tissue", length(which(dataFilt_sub_melt$subtype == "AANormal_Tissue")))
    
    # Create a directory where these plots will be saved
    dir.create("plots", showWarnings = FALSE)
    
    # Add information about if there is a significant difference in mRNA expression between subtype and normal samples
    # This information is collected from the DEA
    
    
    # In case the gene of interest were significantly differentially expressed in at least one of the subtypes ...
    if(nrow(dataDEGs_all > 0)){
      
      # Add DEA information to the data frame created above the first for-loop 
      # This data frame will in the end contain DEA results for all cancer types + subtypes for the gene
      dataDEGs_all$cancer_type <- rep(cancer, nrow(dataDEGs_all))
      dataDEGs_cancer <- rbind(dataDEGs_cancer, dataDEGs_all)
      
      # Add significance information 
      sig.dif <- dataDEGs_all[which(dataDEGs_all$subtype %in% dataFilt_sub_melt$subtype),]$subtype
      sig.dif <- sort(sig.dif)
      sig.val <- c()
      for(canc in unique(sig.dif)){
        sig.val <- c(sig.val, sort(dataFilt_sub_melt[which(dataFilt_sub_melt$subtype == canc),]$value)[length(which(dataFilt_sub_melt$subtype == canc))] + 0.5)
      }
      label.df <- data.frame(subtype = sig.dif, value = sig.val, sample_type =rep("Tumor", length(sig.dif)))
      
      textexplain <- c(expression("* DEA logFC" >=  0.5))
      
      # Do a boxplot
      ggboxplot(data = dataFilt_sub_melt, "subtype", "value", color = "subtype", 
                legend.title = "Subtypes") +
        theme(axis.text.x = element_text(size = 12, angle = 90))+
        ggtitle(paste0(cancer," GSNOR expression")) +
        xlab("Subtypes") + 
        ylab("gene expression (logCPM)") +
        labs(caption = textexplain)+
        geom_text(data = label.df, label = "*",size = 5)+
        ggsave(filename = paste0("plots/GSNOR_", cancer, "_subtypes_boxplot.pdf")) 
      
      
      # In case the gene of interest were not significantly differentially expressed in any of the subtypes ...
    }else{
      textexplain <- "No significant differential expression"
      
      # Do a boxplot
      ggboxplot(data = dataFilt_sub_melt, "subtype", "value", color = "subtype", legend.title = "Sample types") +
        theme(axis.text.x = element_text(size = 12, angle = 90))+
        ggtitle(paste0(cancer," GSNOR expression")) +
        xlab("Subtypes") + 
        ylab("gene expression (logCPM)") +
        labs(caption = textexplain)+
        ggsave(filename = paste0("plots/GSNOR_", cancer, "_subtypes_boxplot.pdf")) 
    }
  }  
}

# Save data frame
save(dataDEGs_cancer, file = "DEA_tables/GSNOR_DEGs_subtypes.rda")
