### Work from "/data/user/shared_projects/gsnor_pan/TCGA"

library(ggplot2)
library(reshape2)
library(ggsignif)

source("code/TCGA_gene_expression_functions.r")


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

geneName <- c("ADH5")

### Get metled data.frame showing voom transformed gene expression in tumour and normal samples
## cancer_types = list of TCGA cancer types to get gene expression in, for example "BRCA", "LUAD", etc.
## geneName = list of HUGO gene name/names to get expression data from
## expression_filename = name to give the saved gene expression file. Should be .rda file

expression_filename <- "GSNOR_expression.rda"
expression_file <- get_gene_expression(cancer_types, geneName, expression_filename) 


#expression_file <- get(load(file = "DEA_tables/all_voom_dataFilt.rda"))
#expression_file <- melt(expression_file)


### Do a bootstrapping on cancer vs. normal to find out if there is a significant difference
### in expression between the sample types, and add the results to a separate 'bootstrap_file' data frame.
### Do this on each cancer type separately

## Start by creating 'bootstrap_file', 
bootstrap_file <- data.frame(row.names = cancer_types)

## Create two empty columns in 'bootstrap_file', to which significance stars and their positions in the boxplot can be added to 
bootstrap_file$bootstrap <- rep("", nrow(bootstrap_file))
bootstrap_file$position <- rep(0, nrow(bootstrap_file))

for(cancer in cancer_types){
  
  ## Get the positions tumor/normal samples have for this cancer type, to later and extract their expression values
  cancer.pos <- which(expression_file$cancer_type == cancer)
  tumor.pos <- which(expression_file$sample_type == "Tumor")
  normal.pos <- which(expression_file$sample_type == "Normal")
  tumor.cancer.pos <- cancer.pos[which(cancer.pos %in% tumor.pos)]
  normal.cancer.pos <-  cancer.pos[which(cancer.pos %in% normal.pos)]
  
  ## Do a bootstrap
  ## group1 = numeric vector of expression values from cancer samples
  ## group2 = numeric vector of expression values from normal samples
  ## p.value = p-value significance limit, for example 0.05
  
  group1 <- expression_file$value[tumor.cancer.pos] 
  group2 <- expression_file$value[normal.cancer.pos]
  p.value <- 0.05
  boot_sig <- bootstrap_groups(group1, group2, p.value)
  boot.pos <- which(rownames(bootstrap_file) == cancer)
  if(boot_sig == "significant"){
    bootstrap_file$bootstrap[boot.pos] <- "*"
  }
  
  ### Make positions for significance stars in boxplot
  cancerval <- sort(expression_file[cancer.pos,]$value)[length(expression_file[cancer.pos,]$value)] + 1
  bootstrap_file$position[boot.pos] <- cancerval
}


### Get differential gene expression values
## cancer_types = list of TCGA cancer types to get gene expression in, for example "BRCA", "LUAD", etc.
## geneName = list of HUGO gene name/names to get expression data from
## DEG_filename = name to give the saved gene expression file. Should be .rda file
DEG_filename <- "GSNOR_DEA.rda"
dataDEGs <- get_DEGs(cancer_types, geneName, DEG_filename)

### Do a boxplot showing expression values in tumor and normal samples
## expression_file = output file from 'get_gene_expression'
## boxplot_filename = name to give the saved boxplot file. Should be .pdf or .png file
boxplot_filename <- "GSNOR_expression_boxplot.pdf"
bxp <- boxplot_tumor_vs_normal(expression_file, dataDEGs, boxplot_filename, bootstrap_file)
bxp

### Do a barplot showing differential gene expression values in tumor vs normal
## dataDEGs = output file from 'get_DEGs'
## barplot_filename = name to give the saved barplot file. Should be .pdf file
barplot_filename <- "GSNOR_DEA_barplot.pdf"
brp <- barplot_DEGs(dataDEGs,barplot_filename)
brp