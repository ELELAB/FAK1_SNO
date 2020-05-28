library(SummarizedExperiment)
library(TCGAbiolinks)


source("code/TCGA_mutations_functions.r")

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
# Which pipeline to use to find muatations. Possible pipelines are: muse, varscan2, somaticsniper and mutect2
pipeline <- c("mutect2")

# Pick genes of interest to look for mutations in
geneNames <- c("ADH5", "PTK2")

# Create an empty data frame to fill with mutation information for every cancer type
all_muts_pipe <- data.frame()

for(cancer in cancer_types){
  # Find mutations in that cancer type using the selected pipeline
  query_maf <- mutect2_tcga_maf(cancer, pipeline)
  
  # See if any mutations are found in the genes of interest
  query_gen <- query_maf[which(query_maf$Hugo_Symbol %in% geneNames),]
  
  # Add information to table about cancer type and pipeline
  query_gen$Cancer <- rep(cancer, nrow(query_gen))
  query_gen$Pipeline <- rep(pipeline, nrow(query_gen))
  
  # Save to the data frame created above
  all_muts_pipe <- rbind(all_muts_pipe, query_gen)
  
}

# Save results
possposs <- c(1,5,6,7,8,9,35,36,52,53,54,55,56,57,61, 69,73,74,94, 121, 122)

all_muts_pipe_good <- all_muts_pipe[,possposs]
all_muts_pipe_good <- all_muts_pipe_good[,c(21,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]
write.csv(all_muts_pipe_good, file = "tables/gsnor_fak1_mutations.csv", quote = F)

# Save separate files for every gene
for(gene in unique(all_muts_pipe$Hugo_Symbol)){
  all_muts_pipe_good_gene <- all_muts_pipe_good[which(all_muts_pipe_good$Hugo_Symbol == gene),]
  write.csv(all_muts_pipe_good_gene, file = paste0("tables/", gene, "_mutations.csv"), quote = F)
}