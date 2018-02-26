library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)

install_url("https://github.com/satijalab/seurat/releases/download/v2.0.0/Seurat_2.0.0_R3.4.tgz",
            binary = TRUE)

# read the 10X data
pbmc.data <- Read10X(data.dir = "datasets/filtered_gene_bc_matrices/hg19/")
class(pbmc.data)
dim(pbmc.data)
pbmc.data[1:6, 1:6]


# summary of total expression per single cell
summary(colSums(pbmc.data))

# check how many genes have at least one transcript in each cell
at_least_one <- apply(pbmc.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")


hist(colSums(pbmc.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")


# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(pbmc.data, 1, function(x) sum(x>0))
table(tmp>=3)

