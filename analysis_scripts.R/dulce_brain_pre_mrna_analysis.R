library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)

# read the 10X data
dat <- Read10X(data.dir = "datasets/filtered_gene_bc_matrices/hg19.premrna/")
class(dat)
dim(dat)
dat[1:6, 1:6]


# summary of total expression per single cell
summary(colSums(dat))

# check how many genes have at least one transcript in each cell
at_least_one <- apply(dat, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")


hist(colSums(dat),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")


# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(dat, 1, function(x) sum(x>0))
table(tmp>=3)

# all cells have at least 200 detected genes
keep <- tmp >= 3
tmp <- dat[keep,]
at_least_one <- apply(tmp, 2, function(x) sum(x>0))
summary(at_least_one)

s_dat <- CreateSeuratObject(raw.data = dat,
                           min.cells = 3,
                           min.genes = 200,
                           project = "hg19_pre_mrna")

s_dat

# Count mito genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = s_dat@data), value = TRUE)
length(mito.genes)
percent.mito <- Matrix::colSums(s_dat@raw.data[mito.genes, ]) /
  Matrix::colSums(s_dat@raw.data)

# Histogram of mitochondrial percentage
hist(percent.mito)

# Add the mito data to the seurat object
head(s_dat@meta.data) 

s_dat <- AddMetaData(object = s_dat,
                    metadata = percent.mito,
                    col.name = "percent.mito")
head(s_dat@meta.data) 


VlnPlot(object = s_dat,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = s_dat, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = s_dat, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')


# Filter and normalize

# manual check; I already know all cells have >200 genes
table(s_dat@meta.data$percent.mito < 0.05 & s_dat@meta.data$nGene<2500)

# perform the filtering using FilterCells()
s_dat <- FilterCells(object = s_dat,
                    subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf),
                    high.thresholds = c(6000, 0.05))

hist(colSums(s_dat@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")


s_dat <- NormalizeData(object = s_dat,
                       normalization.method = "LogNormalize",
                       scale.factor = 1e4)

hist(colSums(s_dat@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")

# Find variable genes
# the variable genes slot is empty before the analysis
s_dat@var.genes


s_dat <- FindVariableGenes(object = s_dat,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR, do.text=FALSE)



