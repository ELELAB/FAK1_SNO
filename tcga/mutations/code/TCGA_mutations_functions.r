
mutect2_tcga_maf <- function(cancer, pipe){
  dir.create("data", showWarnings = FALSE)
  dir.create("tables", showWarnings = FALSE)
  query.maf <- GDCquery_Maf(cancer, pipelines = pipe, directory = "data/GDCdata")
  #save(query.maf, file = paste0("tables/",cancer, "_somatic_variations.rda"))
  return(query.maf)
}
