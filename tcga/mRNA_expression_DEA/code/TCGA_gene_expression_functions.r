get_gene_expression <- function(cancer_types, geneName, file_name){
  dir.create("DEA_tables", showWarnings = FALSE)
  dataFilt_all <- data.frame()
  for(cancer in cancer_types){
    dataFilt <- get(load(paste0("data/dataFilt/",cancer,"_dataFilt.rda")))
    dataFilt_voom <- limma::voom(dataFilt)
    dataFilt_voom <- as.data.frame(dataFilt_voom[[1]])
    dataFilt_voom_sub <- dataFilt_voom[which(rownames(dataFilt_voom) %in% geneName),] 
    if(nrow(dataFilt_voom_sub) >0){
      dataFilt_voom_sub_t <- as.data.frame(t(dataFilt_voom_sub))
      dataFilt_voom_sub_t$cancer_type <- rep(cancer,nrow(dataFilt_voom_sub_t))
      dataFilt_voom_sub_t$sample_type <- rep("NA", nrow(dataFilt_voom_sub_t))
      for(rownumb in 1:nrow(dataFilt_voom_sub_t)){
        row_name <- rownames(dataFilt_voom_sub_t)[rownumb]
        subtype <- substr(unlist(stringr::str_split(row_name,"-"))[4], 1, 2)
        if(subtype == "01"){
          dataFilt_voom_sub_t[rownumb,]$sample_type <- "Tumor"
        }else if(subtype == "11"){
          dataFilt_voom_sub_t[rownumb,]$sample_type <- "Normal"
        }
        
      }
      dataFilt_all <- rbind(dataFilt_all, dataFilt_voom_sub_t)
    }
  }
  save(dataFilt_all, file = paste0("DEA_tables/", file_name))
  dfa_melt <- melt(dataFilt_all)
  return(dfa_melt)
}

bootstrap_groups <- function(group1, group2, p.value){
  library(simpleboot)
  library(dplyr)
  
  bootstrapped <- two.boot(group1, group2, mean, 1000)
  stderr <- sort(c(quantile(bootstrapped$t, probs = 1- p.value/2),quantile(bootstrapped$t, probs = p.value/2)))
  
  if(isTRUE(between(0, stderr[1], stderr[2]))){
    cat("There is no significant difference between the two groups", "\n")
    sigval <- "not.significant"
  }else{
    cat("There is a significant difference between the two groups", "\n")
    sigval <- "significant"
  }
  return(sigval)
}

boxplot_tumor_vs_normal <- function(expression_file, dataDEGs, boxplot_filename, bootstrap_file){
  dir.create("plots", showWarnings = FALSE)
  sig.dif <- dataDEGs[which(dataDEGs$cancer_type %in% expression_file$cancer_type),]$cancer_type
  sig.dif <- sort(sig.dif)
  sig.val <- c()
  for(canc in unique(sig.dif)){
    sig.val <- c(sig.val, sort(expression_file[which(expression_file$cancer_type == canc),]$value)[length(which(expression_file$cancer_type == canc))] + 0.5)
  }
  label.df <- data.frame(cancer_type = sig.dif, value = sig.val, sample_type = rep(c("Tumor"), length(sig.val)))
  textexplain <- c(expression("* DEA logFC" >=  0.5))
  ggboxplot(expression_file, "cancer_type", "value", color = "sample_type",
            palette = c("#00AFBB", "#E7B800"), legend.title = "Sample types") +
    theme(axis.text.x = element_text(size = 10, angle = 90))+
    ggtitle("GSNOR expression") +
    xlab("Cancer types") + 
    ylab("gene expression (logCPM)") +
    labs(caption = textexplain)+
    geom_text(data = label.df, label = "*",size = 7)+
    ggsave(filename = paste0("plots/", boxplot_filename))
}

get_DEGs <- function(cancer_types, geneName, file_name){
  dir.create("DEA_tables", showWarnings = FALSE)
  dataDEGs_all <- data.frame()
  for(cancer in cancer_types){
    dataDEGs <- get(load(paste0("/data/user/shared_projects/TCGA_data/DEA/cancer_types/",cancer,"/",cancer,"_dataDEGs.rda")))
    dataDEGs_pos <- dataDEGs[which(dataDEGs$logFC > 0.5),]
    dataDEGs_neg <- dataDEGs[which(dataDEGs$logFC < -0.5),]
    dataDEGs <- rbind(dataDEGs_pos, dataDEGs_neg)
    dataDEGs_gene <- dataDEGs[which(rownames(dataDEGs) %in% geneName),]
    if(nrow(dataDEGs_gene) > 0){
      dataDEGs_gene$cancer_type <- rep(cancer, nrow(dataDEGs_gene))
      dataDEGs_all <- rbind(dataDEGs_all, dataDEGs_gene)
    }
  }
  save(dataDEGs_all, file = paste0("DEA_tables/", file_name))
  return(dataDEGs_all)
}

barplot_DEGs <- function(dataDEGs, barplot_filename){
  dir.create("plots", showWarnings = FALSE)
  dataDEGs$colour <- rep("NA",nrow(dataDEGs))
  for(rownumb in 1:nrow(dataDEGs)){
    if(dataDEGs[rownumb,]$logFC > 0){
      dataDEGs[rownumb,]$colour <- "red"
    }else{
      dataDEGs[rownumb,]$colour <- "blue"
    }
  }
  pdf(paste0("plots/", barplot_filename))
  bp <- barplot(dataDEGs$logFC,
                names.arg = dataDEGs$cancer_type,
                cex.axis=0.9,
                cex.names=0.9,
                ylim = c(-1.5, 1.5),
                main = paste("GSNOR differential expression"),
                ylab = "logFC",
                axes = TRUE, 
                col = dataDEGs$colour,
                las = 2)
  dev.off()
  return(bp)
}