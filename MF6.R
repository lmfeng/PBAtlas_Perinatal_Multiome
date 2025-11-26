library(gplots)
library(ComplexHeatmap)
library(circlize) 
library(WGCNA)
library(limma)
library(RColorBrewer)
set.seed(100)

disease <- read.table("/PBA/disease_science_lmf.txt", row.names = 1, header=TRUE, sep = '\t')
col_names <- colnames(disease)
gene_lists <- list()

for (i in 1:ncol(disease)) {
  genes <- rownames(disease)[disease[, i] == 1]
  gene_lists[[col_names[i]]] <- genes
}

disease2 <- read.table("/PBA/disease_biorxiv_2842_only.txt", row.names = 1, header=TRUE, sep = '\t')
col_names <- colnames(disease2)
gene_lists2 <- list()

for (i in 1:ncol(disease2)) {
  genes <- rownames(disease2)[disease2[, i] == 1]
  gene_lists2[[col_names[i]]] <- genes
}

disease_list <- c(gene_lists, gene_lists2)

nddList = disease_list

gexpr <- NULL
meta <- NULL
load("/PBA/PBAtlas.geneExpr.norm.RData")
gexpr.orig = gexpr
meta.orig = meta
regionRef = c("CERE","STR","HIP","THA","PFC","TMP","OCC")
ageRef = c("E76","E85","E94","E104","E109","P0","P3","P30")

##
idx.use = which(meta.orig$Age %in% c("E76","E85","E94"))
gexpr.lf = gexpr.orig[, idx.use]
calcor.lf = get_ndd_cor(gexpr.lf, nddList)
##
idx.use = which(meta.orig$Age %in% c("E104","E109","P0","P3"))
gexpr.peri = gexpr.orig[, idx.use]
calcor.peri = get_ndd_cor(gexpr.peri, nddList)


rpkm <- NULL
meta <- NULL
load("/PBA/BrainSpan.SestanScience2018.RData")
gsym = substring(rownames(rpkm), 17)
rpkm = rpkm[!duplicated(gsym), ]
rownames(rpkm) = gsym[!duplicated(gsym)]
meta  = meta 

###
regionRef.use = c("OFC","DFC","VFC","MFC","ITC","OC","HIP","V1C","STR","MD","CBC")
idx.use = which(meta$Regioncode %in% regionRef.use)
meta.bs = meta[idx.use, ]
gexpr.bs = rpkm[, idx.use]

####normaliztion
gexpr.bs <- log2(as.matrix(gexpr.bs) + 1)
gexpr.bs <- normalizeBetweenArrays(gexpr.bs, method = "quantile")

###calculae correlation
##
idx.use = which(meta.bs$Period %in% c(2,3,4))
gexpr.ef = gexpr.bs[, idx.use]
calcor.ef = get_ndd_cor(gexpr.ef, nddList)

##
idx.use = which(meta.bs$Period %in% c(4,5,6))
gexpr.mf = gexpr.bs[, idx.use]
calcor.mf = get_ndd_cor(gexpr.mf, nddList)

##
idx.use = which(meta.bs$Period %in% c(11,12,13))
gexpr.adult = gexpr.bs[, idx.use]
calcor.adult = get_ndd_cor(gexpr.adult, nddList)

col_fun <- viridis::inferno(256)
##
column_anno <- HeatmapAnnotation(
  df = data.frame(Stage = rep("Early-fetal", ncol(calcor.ef))),
  col = list(Stage = c("Early-fetal" = "#FEE0D2")),
  show_annotation_name = FALSE
)
ht1 <- Heatmap(calcor.ef, name = "Early-fetal", col = col_fun,column_names_rot = -60,
               column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 12),top_annotation = column_anno)


row_order <- row_order(ht1); column_order <- column_order(ht1)

column_anno <- HeatmapAnnotation(
  df = data.frame(Stage = rep("Mid-fetal", ncol(calcor.ef))),
  col = list(Stage = c("Mid-fetal" = "#FC9272")),
  show_annotation_name = FALSE
)
ht2 <- Heatmap(calcor.mf, name = "Mid-fetal", col = col_fun, cluster_rows = FALSE, 
               cluster_columns = FALSE, row_order = row_order, 
               column_order = column_order,show_column_names = T,column_names_rot = -60,
               row_names_gp = gpar(fontsize = 12),
               top_annotation = column_anno)

column_anno <- HeatmapAnnotation(
  df = data.frame(Stage = rep("Later-fetal", ncol(calcor.ef))),
  col = list(Stage = c("Later-fetal" = "#FB6A4A")),
  show_annotation_name = FALSE
)
ht3 <- Heatmap(calcor.lf, name = "Later-fetal", col = col_fun, cluster_rows = FALSE,
               cluster_columns = FALSE, row_order = row_order, 
               column_order = column_order,show_column_names = T,column_names_rot = -60,
               row_names_gp = gpar(fontsize = 12),
               top_annotation = column_anno)

column_anno <- HeatmapAnnotation(
  df = data.frame(Stage = rep("Perintal", ncol(calcor.ef))),
  col = list(Stage = c("Perintal" = "#CB181D")),
  show_annotation_name = FALSE
)
ht4 <- Heatmap(calcor.peri, name = "Perintal", col = col_fun, cluster_rows = FALSE, 
               cluster_columns = FALSE, row_order = row_order, 
               column_order = column_order,show_column_names = T,column_names_rot = -60,
               row_names_gp = gpar(fontsize = 12),
               top_annotation = column_anno)

column_anno <- HeatmapAnnotation(
  df = data.frame(Stage = rep("Adult", ncol(calcor.ef))),
  col = list(Stage = c("Adult" = "#67000D")),
  show_annotation_name = FALSE
)
ht5 <- Heatmap(calcor.adult, name = "Adult", col = col_fun, cluster_rows = FALSE, 
               cluster_columns = FALSE, row_order = row_order, column_order = column_order,
               show_column_names = T,column_names_rot = -60,
               row_names_gp = gpar(fontsize = 12),
               top_annotation = column_anno)

par(omi = c(0.1, 0.1, 0.1, 0.1))
ht1 + ht2 + ht3 + ht4 + ht5
