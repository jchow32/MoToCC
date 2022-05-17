library(Seurat)  # 3.2.1
library(dplyr)
library(patchwork)
library(Matrix)
library(data.table)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(gridGraphics)
library(RColorBrewer)

make_2d_tsne = function(input_file, output_name) {
  selected = read.table(input_file, header = TRUE, sep = "\t", row.names="cell")
  # these are the selected cells 
  to_highlight = rownames(selected)
  
  g = DimPlot(in_file, reduction = "tsne_2D", cells.highlight=to_highlight) + 
    ggtitle(output_name) + NoLegend()
  return(g)
}

make_3d_tsne = function(input_file, output_name) {
  selected = read.table(input_file, header = TRUE, sep = "\t", row.names="cell")
  # these are the selected cells 
  to_highlight = rownames(selected)
  
  g = DimPlot(in_file, reduction = "tsne_3D", cells.highlight=to_highlight) + 
    ggtitle(output_name) + NoLegend()
  return(g)
}

make_umap = function(input_file, output_name) {
  selected = read.table(input_file, header = TRUE, sep = "\t", row.names="cell")
  # these are the selected cells 
  to_highlight = rownames(selected)
  
  g = DimPlot(in_file, reduction = "umap", cells.highlight=to_highlight) + 
    ggtitle(output_name) + NoLegend()
  return(g)
}

# To save normalized and scaled expression data for all genes
in_file = readRDS("local.rds")
in_file <- NormalizeData(in_file)
in_file <- ScaleData(in_file)  # mean 0, STD 1
write.table(t(as.matrix(GetAssayData(object = in_file, slot = "scale.data"))), sep=",", file='norm_scale_exp.csv', quote=FALSE)

# To save cell-cell similarity and KNN graph
in_file = readRDS("local.rds")
in_file <- NormalizeData(in_file)
in_file <- FindVariableFeatures(in_file) 
in_file <- ScaleData(in_file)  # mean 0, STD 1

in_file <- RunPCA(in_file, features = VariableFeatures(object = in_file), npcs = 20)
in_file <- FindNeighbors(in_file, dims = 1:20, compute.SNN = TRUE)
write.table(
  as.matrix(in_file$RNA_nn), sep=",", 
  file="KNN_graph.csv", quote=FALSE)
write.table(
  as.matrix(in_file$RNA_snn), sep=",", 
  file="SNN_graph.csv", quote=FALSE)

in_file <- FindClusters(in_file)
in_file <- RunTSNE(in_file, dims=1:20, reduction.name="tsne_2D")
in_file <- RunTSNE(in_file, dims=1:20, dim.embed=3, reduction.name="tsne_3D")
#  in_file <- RunUMAP(in_file, dims=1:20)  # already calculated for this data

exit()

######
# To print the 2D t-SNE with the paper's identity labels
original = DimPlot(in_file, reduction="tsne_2D")
file_name_or = "original_2DtSNE.png"
png(file_name_or)
print(original)
dev.off()

# Examples for how to make 2D and 3D t-SNE and UMAP plots
# the solution file from C-TULPSE is named ${output_string}_${k}_LPonLSCC.txt

tsne_2_250 = make_2d_tsne('SLC2A2_250_LPonLSCC.txt', '250', input_data)
tsne_3_250 = make_3d_tsne('SLC2A2_250_LPonLSCC.txt', '250', input_data)
umap_250 = make_umap('SLC2A2_250_LPonLSCC.txt', '250', input_data)
pdf("Example_SLC2A2_250.pdf")
grid.arrange(tsne_2_250, tsne_3_250, umap_250, nrow=1)
dev.off()
