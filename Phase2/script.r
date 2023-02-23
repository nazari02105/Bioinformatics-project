# This is the code for the Phase 2 of the Bioinformatics project
# at Sharif University of Technology in Fall Semester of 2022
# By Ali Nazari, Parham Bateni, Masih Najafi
setRepositories() # include number 2 repositories
install.packages(c("BiocManager","GEOquery", "limma", "pheatmap", "ggplot2", "gplots", "reshape2", "plyr", "umap"))
BiocManager::install('GEOquery', force = TRUE)
library(GEOquery)
library(pheatmap)

Sys.setenv(VROOM_CONNECTION_SIZE=500072)
setwd("D://Codes/R RStudio/Bio_Phase2")
options('download.file.method.GEOquery' = 'libcurl')
gset<-getGEO(GEO='GSE48558',GSEMatrix = TRUE, AnnotGPL=TRUE)
gset_copy = gset
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
ex <- exprs(gset)
normal_cols = which(sml == "0")
cancer_cols = which(sml == "1")
main_normal = gset[ ,normal_cols]
main_cancer = gset[ ,cancer_cols]
main_normal = exprs(main_normal)
main_cancer = exprs(main_cancer)
all_amls <- main_cancer
all_normals <- main_normal
CD34_list <- c("GSM1180849", "GSM1180853", "GSM1180857")

aml_group <- all_amls
cd34_group <- all_normals[, CD34_list]
pheatmap(cor(aml_group, cd34_group), main = "aml VS. cd34")

dim(aml_group)
dim(cd34_group)
row_names = rownames(aml_group)
genes_with_p_value_less_than_0.05 = c()
for(i in 1:32321) {
  p_value = t.test(aml_group[row_names[i],], cd34_group[row_names[i],])$p.value
  if(p_value < 0.05){
    genes_with_p_value_less_than_0.05 = c(genes_with_p_value_less_than_0.05, row_names[i])
  }
}
print(genes_with_p_value_less_than_0.05)
length(genes_with_p_value_less_than_0.05)

more_in_aml = c()
more_in_cd34 = c()
for(i in 1:7342) {
  aml_mean = mean(aml_group[genes_with_p_value_less_than_0.05[i],])
  cd34_mean = mean(cd34_group[genes_with_p_value_less_than_0.05[i],])
  if (aml_mean > cd34_mean){
    more_in_aml = c(more_in_aml, genes_with_p_value_less_than_0.05[i])
  }
  if (aml_mean <= cd34_mean){
    more_in_cd34 = c(more_in_cd34, genes_with_p_value_less_than_0.05[i])
  }
}
print(more_in_aml)
length(more_in_aml)
print(more_in_cd34)
length(more_in_cd34)

unique_genes <- unique(genes_with_p_value_less_than_0.05)
write.table(unique_genes, file = "./AML_CD34_Up.txt", quote = F, row.names = F, col.names = F)
