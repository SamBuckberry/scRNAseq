library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(reshape2)
library(stringr)

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

# vector of variable genes
head(s_dat@var.genes)

length(s_dat@var.genes)

# mean and variance of genes are stored pbmc@hvg.info
head(s_dat@hvg.info)



# slot is empty before running ScaleData()
s_dat@scale.data


# build linear model using nUMI and percent.mito
s_dat <- ScaleData(object = s_dat,
                  vars.to.regress = c("nUMI", "percent.mito"))

class(s_dat@scale.data)
s_dat@scale.data[1:6, 1:6]


s_dat <- RunPCA(object = s_dat,
               pc.genes = s_dat@var.genes,
               do.print = TRUE,
               pcs.print = 1:5,
               genes.print = 5)

PrintPCAParams(s_dat)

PrintPCA(object = s_dat, pcs.print = 1:2,
         genes.print = 5, use.full = FALSE)

# visualise top genes associated with principal components
VizPCA(object = s_dat, pcs.use = 1:2)

PCAPlot(object = s_dat, dim.1 = 1, dim.2 = 2)


# the results of the projected PCA can be explored by setting use.full=TRUE in the functions above
s_dat <- ProjectPCA(object = s_dat, do.print = FALSE)

# 500 cells
PCHeatmap(object = s_dat,
          pc.use = 1,
          cells.use = 500,
          do.balanced = TRUE,
          label.columns = FALSE)

# All cells
PCHeatmap(object = s_dat,
          pc.use = 1,
          do.balanced = TRUE,
          label.columns = FALSE)

# Determine number of PC's to use

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
system.time(
  s_dat <- JackStraw(object = s_dat,
                    num.replicate = 100,
                    do.print = FALSE)
)

JackStrawPlot(object = s_dat, PCs = 1:15)

PCElbowPlot(object = s_dat)


s_dat <- FindClusters(object = s_dat,
                     reduction.type = "pca",
                     dims.use = 1:8,
                     resolution = 0.6,
                     print.output = 0,
                     save.SNN = TRUE)

# use PrintFindClustersParams() to print summary
# of parameters used to FindClusters()
PrintFindClustersParams(object = s_dat)


s_dat <- RunTSNE(object = s_dat,
                dims.use = 1:8,
                do.fast = TRUE)

TSNEPlot(object = s_dat, do.label = TRUE)


### Find markers of clusters

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = s_dat,
                                ident.1 = 1,
                                min.pct = 0.25)

head(cluster1.markers)

# find all markers distinguishing cluster 5 from clusters 2 and 4
cluster5.markers <- FindMarkers(object = s_dat,
                                ident.1 = 5,
                                ident.2 = c(2,4),
                                min.pct = 0.25)
head(cluster5.markers)


# find markers for every cluster compared to all remaining cells, report only the positive ones
s_dat.markers <- FindAllMarkers(object = s_dat,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               thresh.use = 0.25)

head(s_dat.markers)


### Or the top 2 genes for each cluster

s_dat.markers %>% group_by(cluster) %>% top_n(2, avg_diff)

## Stats tests
levels(s_dat@ident)
table(s_dat@ident)

my_bimod <- FindMarkers(object = s_dat,
                        ident.1 = 1,
                        thresh.use = 0.25,
                        test.use = "bimod",
                        only.pos = TRUE)

my_roc <- FindMarkers(object = s_dat,
                      ident.1 = 1,
                      thresh.use = 0.25,
                      test.use = "roc",
                      only.pos = TRUE)

my_t <- FindMarkers(object = s_dat,
                    ident.1 = 1,
                    thresh.use = 0.25,
                    test.use = "t",
                    only.pos = TRUE)

my_tobit <- FindMarkers(object = s_dat,
                        ident.1 = 1,
                        thresh.use = 0.25,
                        test.use = "tobit",
                        only.pos = TRUE)


# identical set of genes
dim(my_bimod)
dim(my_roc)
dim(my_t)
dim(my_tobit)


my_gene <- row.names(my_bimod)
a <- 1:length(my_gene)
b <- match(my_gene, row.names(my_roc))
c <- match(my_gene, row.names(my_t))
d <- match(my_gene, row.names(my_tobit))


# bimod vs. bimod
cor(a, a, method = "spearman")

# bimod vs. roc
cor(a, b, method = "spearman")

# bimod vs. t
cor(a, c, method = "spearman")

# bimod vs. tobit
cor(a, d, method = "spearman")

par(mfrow=c(2,2))
barplot(a, main = 'bimod')
barplot(b, main = 'roc')
barplot(c, main = 't')
barplot(d, main = 'tobit')



# Visualise some markers

VlnPlot(object = s_dat, features.plot = c("PLP1", "ST18"))


# And with raw UMI counts

# you can plot raw UMI counts as well
VlnPlot(object = s_dat,
        features.plot = c("FGF13", "RELN"),
        use.raw = TRUE,
        y.log = TRUE)

## And the tSNE with marker gene expression
FeaturePlot(object = s_dat,
            features.plot = c("PLP1", "ST18", "FGF13", "RELN", "LHFPL3", "PCDH15", "NXPH1", "GRIK1", "CBLN2"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

FeaturePlot(object = s_dat,
            features.plot = head(row.names(my_tobit), 9),
            cols.use = c("grey", "blue"))



## Heatmaps
head(s_dat.markers)


s_dat.markers %>%
  group_by(cluster) %>%
  top_n(10, avg_diff) -> top10

head(top10)


DoHeatmap(object = s_dat,
          genes.use = top10$gene,
          order.by.ident = TRUE,
          slim.col.label = TRUE,
          remove.key = TRUE)


###### Make barplot of marker genes

## Get top 2 markers
top_markers <- s_dat.markers %>% group_by(cluster) %>% top_n(1, avg_diff)

markers <- top_markers$gene
object <- s_dat


dim(object@raw.data)
dim(object@data)

barplot_markers <- function(object, markers, y_text_size=8){

  
  # Check if there are replicate marker ids
  stopifnot(all(unique(markers) == markers))
  
  # Get the ID's of the cells to include after filtering
  cell_include <- colnames(object@data)
  
  ## Get the raw UMI data
  dat <- object@raw.data[rownames(object@raw.data) %in% markers, ] %>%
    as.matrix() %>% t() %>% data.frame()
  dat <- dat[rownames(dat) %in% cell_include, ]
  
  stopifnot(all(rownames(dat) == names(object@ident)))
  
  # Add the cluster information
  stopifnot(length(object@ident) == nrow(dat))
  dat$cluster <- as.numeric(object@ident)
  
  # Reorder by cluster
  dat <- dat[order(dat$cluster), ]
  dat$cell <- dat$cell <- 1:nrow(dat)

  dat <- melt(dat, id.vars = c("cluster", "cell"))
  
  #dat$variable <- str_replace(dat$variable, pattern = "-", replacement = "_")
  #marker_correct <- str_replace(markers, pattern = "-", replacement = "_")
  
  dat$variable <- factor(x = dat$variable, levels=markers)

  ggplot(data = dat, aes(x = cell, y = value+0.1, fill=factor(cluster))) +
    geom_col() +
    facet_grid(variable~., scales = "free_y") +
    ylab(label = "UMI count") +
    labs(x = "", fill = "Cluster") +
    theme(strip.text.y = element_text(angle = 0, hjust = 0),
          strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = y_text_size))
  
  
}

barplot_markers(object = s_dat, markers =c("NXPH1", "PLP1", "TSHZ2"))


library(pheatmap)
heatmap_cells <- function(object){

  # Get variable gene IDs  
  var_genes <- object@var.genes
  
  # Get the ID's of the cells to include after filtering
  cell_include <- colnames(object@data)

  # Get UMI data for variable genes
  dat <- object@raw.data[rownames(object@raw.data) %in% var_genes, ] %>%
    as.matrix() %>% t() %>% data.frame()
  dat <- dat[rownames(dat) %in% cell_include, ]
  
  # Add the cluster information
  stopifnot(all(names(object@ident) == rownames(dat)))
  cluster <- as.numeric(object@ident)
  ord <- order(cluster)
  
  # Reorder by cluster
  dat <- dat[ord, ]
  
  pheatmap(as.matrix(dat), cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE)
  
  
  
}






qplot(test$cell, test$value, colour=factor(test$cluster),
      geom = 'bar', stat = 'identity')

  facet_grid(variable~., space = "free_y")
gg_bar
