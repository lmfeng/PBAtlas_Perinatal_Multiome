library(dplyr)
library(tidyverse)
library(readr)
library(reshape2)
library(cluster)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(circlize)

source("/PBA/Function.R")

#####————3B and 3D————#####
gexpr <- read.table("/PBA/geneRPKM_PBA.txt", header = TRUE, sep = "\t", row.names = 1)
gexpr <- t(gexpr)
meta <-  read.table("/PBA/original_data/PBAtlas/samples_all182.txt", header = TRUE, row.names = 1)
reg_ref <- c("PFC","TMP","OCC","STR","THA","HIP","CERE")
age_ref <- c("E76","E85","E94","E104","E109","P0","P3","P30")

hiexpindex = which(apply(gexpr,1,function (x) {y=length(which(x >=1))}) >= 5)
gexpr = gexpr[hiexpindex,]

sampList = colnames(gexpr)
geneList = as.character(rownames(gexpr))
bakeup = gexpr

gexpr = as.numeric(gexpr)
dim(gexpr) = dim(bakeup)
myAve = apply(gexpr,1, mean);
myStd = apply(gexpr,1, sd)

recIndex = c();
for(i in 1:nrow(gexpr)){
  eachExp = as.numeric(gexpr[i,]);
  eachIndex = which(eachExp > (max(eachExp)*0.75));
  eachNumber = length(eachIndex);
  if(eachNumber > 2){
    eachSample = sampList[eachIndex];
    reg = meta[eachSample,"Region"]
    age = meta[eachSample,"Age"]
    rec = data.frame(reg = reg,
                     age = age)
    x = table(as.character(rec$reg))
    y = table(as.character(rec$age))
    if(length(x) > 1 && length(y) < nrow(rec))	recIndex = c(recIndex,i)
  }
}

recIndex = intersect(which(myStd >1),recIndex)
gexpr = gexpr[recIndex,];
geneList = geneList[recIndex];
length(geneList)
# 7891

##final data
exprMat <- bakeup[geneList,]
exprMat = log2(exprMat+1)
meta <- meta[,c(1,8)]

###Region RNA
heatmap_list <- list()
colors <- rev(brewer.pal(5, "RdYlBu"))
expected_order <- c("E76", "E85", "E94", "E104", "E109", "P0", "P3", "P30")

heatmap_list <- list()
for (region in reg_ref) {
  age_samples <- lapply(expected_order, function(age) {
    colnames(exprMat)[meta$Age == age & meta$Region == region]
  })
  names(age_samples) <- expected_order
  
  age_avg <- lapply(age_samples, function(samples) {
    if (length(samples) > 1) {
      rowMeans(exprMat[, samples, drop = FALSE]) 
    }else {
      exprMat[, samples, drop = FALSE]
    }
  })
  age_avg <- age_avg[!sapply(age_avg, is.null)]
  combined_df <- do.call(cbind, age_avg)
  cor_matrix <- cor(combined_df, method = "s")
  rownames(cor_matrix) <- age_ref;colnames(cor_matrix) <- age_ref
  
  min_val <- min(cor_matrix)
  col_fun <- colorRamp2(seq(min_val, 1, length.out = 5), colors)
  
  heatmap_list[[region]] <- Heatmap(
    cor_matrix,
    row_order = age_ref, 
    column_order = age_ref,
    col = col_fun,
    name = "correlation",
    column_title = paste(region),
    show_column_names = T,
    show_row_names = T,
  )
}
pdf("/PBA/MF3_SRDA.hm.g.pdf",height = 3.2,width = 19)
print(heatmap_list[[1]]+heatmap_list[[2]]+heatmap_list[[3]]+heatmap_list[[4]]+heatmap_list[[5]]+heatmap_list[[6]]+heatmap_list[[7]])
dev.off()

###Region Protein
load("/home/wangzm/Project/PBA/original_data/PBAtlas/PBA_proteins/new_20241227/ProteinBatchData.RData")
expr = batch_data
meta = pmeta

hiexpindex = which(apply(expr,1,function (x) {y=length(which(x >=0))}) >= 3)
exprMat = expr[hiexpindex,]

colors <- rev(brewer.pal(5, "RdYlBu"))
expected_order <- c("E76", "E85", "E94", "E104", "E109", "P0", "P3", "P30")

heatmap_list <- list()
for (region in reg_ref) {
  age_samples <- lapply(expected_order, function(age) {
    colnames(exprMat)[meta$Age == age & meta$Region == region]
  })
  names(age_samples) <- expected_order
  
  age_avg <- lapply(age_samples, function(samples) {
    if (length(samples) > 1) {
      rowMeans(exprMat[, samples, drop = FALSE]) 
    }else {
      exprMat[, samples, drop = FALSE]
    }
  })
  age_avg <- age_avg[!sapply(age_avg, is.null)]
  combined_df <- do.call(cbind, age_avg)
  cor_matrix <- cor(combined_df, method = "spearman")
  rownames(cor_matrix) <- age_ref;colnames(cor_matrix) <- age_ref
  
  min_val <- min(cor_matrix)
  col_fun <- colorRamp2(
    breaks = seq(min_val, 1, length.out = 5),
    colors = colors
  )
  
  heatmap_list[[region]] <- Heatmap(
    cor_matrix,
    row_order = age_ref, 
    column_order = age_ref,
    col = col_fun,
    name = "correlation",
    column_title = paste(region),
    show_column_names = T,
    show_row_names = T,
  )
}

pdf("/PBA/MF3_SRDA.hm.p.pdf",height = 3.2,width = 19)
heatmap_list[[1]]+heatmap_list[[2]]+heatmap_list[[3]]+heatmap_list[[4]]+heatmap_list[[5]]+heatmap_list[[6]]+heatmap_list[[7]]
dev.off()

#####————3E and 3F————#####
regions <- c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE")
file_names <- c("E85vs.E76", "E85vs.E94", "E104vs.E94", "E104vs.E109","P0vs.E109", "P0vs.P3","P3vs.P30")
dir <- "/PBA/SameRegion_across_DiffAge/"
all_reg_cols <- c(
  "PFC" = "#d71345",
  "TMP" = "#f26522",
  "OCC" = "#f7acbc",
  "STR" = "#6950a1",
  "THA" = "#009ad6",
  "HIP" = "#3CB371",
  "CERE" = "#8B4513"
)

region_dfs <- lapply(regions, function(region){
  ProcessDEanalysisResult_sRdA(region = region, data_type = "Gene", pv_cut = 0.05, fc_cut = log2(1.25))
})

pp <- ggplot(combined_df, aes(x = Age_range, y = RowCount, color = Region, group = Region)) +
  geom_smooth(method = "loess", se = F, linewidth = 2, span = 0.7) +
  geom_point(size = 0) +
  scale_color_manual(values = all_reg_cols) +
  labs(title = "", 
       x = "Age", 
       y = "Count") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "grey75"), 
        panel.grid.minor = element_line(color = "grey85"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  geom_vline(xintercept = "E94vsE104", linetype = "twodash", color = "#363636", linewidth = 1.0)

regions <- c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE")
file_names <- c("E76vsE85", "E85vsE94", "E94vsE104", "E104vsE109","E109vsP0", "P0vsP3","P3vsP30")
dir <- "/PBA/SameRegAcrossDiffAge/"
all_reg_cols <- c(
  "PFC" = "#d71345",
  "TMP" = "#f26522",
  "OCC" = "#f7acbc",
  "STR" = "#6950a1",
  "THA" = "#009ad6",
  "HIP" = "#3CB371",
  "CERE" = "#8B4513"
)
region_dfs <- lapply(regions, function(region){
  ProcessDEanalysisResult_sRdA(region = region, data_type = "Protein", pv_cut = 0.05, fc_cut = log2(1.25))
})
combined_df <- bind_rows(region_dfs)
combined_df$File <- factor(combined_df$File,levels = c("E76vsE85","E85vsE94","E94vsE104","E104vsE109","E109vsP0","P0vsP3","P3vsP30"))
pp <- ggplot(combined_df, aes(x = File, y = RowCount, color = Region, group = Region)) +
  geom_smooth(method = "loess", se = F, linewidth = 2, span = 0.7) +
  geom_point(size = 0) +
  scale_color_manual(values = all_reg_cols) +
  labs(title = "DEP(SameRegion across DiffAge)",
       x = "Age",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "grey85"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  geom_vline(xintercept = "E94vsE104", linetype = "twodash", color = "#363636", linewidth = 1.0)

#####————3G and 3H————#####
age_group <- c("E76vsE85","E85vsE94","E94vsE104","E104vsE109","E109vsP0","P0vsP3","P3vsP30")
regions <- c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE")
file_names <- c("E85vs.E76", "E85vs.E94", "E104vs.E94", "E104vs.E109","P0vs.E109", "P0vs.P3","P3vs.P30")
dir <- "/PBA/SameRegion_across_DiffAge/"

region_dfs <- lapply(regions, function(region){
  suppressMessages(ProcessDEanalysisResult_sRdA_TRG(region,"gene",pv_cut = 0.05,fc_cut = log2(1.25),age_group))
}
)

merged_df <- do.call(rbind, region_dfs)
merged_df$Region <- factor(merged_df$Region, levels = regions)
merged_df$Age_group <- factor(merged_df$Age_group, levels = age_group)

spectral_palette <- brewer.pal(9, "Reds")
color_gradient <- colorRampPalette(spectral_palette)(100)

pp <- ggplot(merged_df, aes(x = Age_group, y = Region, size = Allcounts, fill = TRGpct)) +
  geom_point(shape = 21, color = "black", stroke = 0.8) +
  scale_size_continuous(range = c(1, 12)) +
  scale_fill_gradientn(colors = color_gradient) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "",
       x = "Age",
       y = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

age_group <- c("E76vsE85","E85vsE94","E94vsE104","E104vsE109","E109vsP0","P0vsP3","P3vsP30")
regions <- c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE")
file_names <- c("E76vsE85", "E85vsE94", "E94vsE104", "E104vsE109","E109vsP0", "P0vsP3","P3vsP30")
dir <- "/PBA/SameRegAcrossDiffAge/"

region_dfs <- lapply(regions, function(region){
  suppressMessages(ProcessDEanalysisResult_sRdA_TRG(region,"protein",pv_cut = 0.05,fc_cut = log2(1.25),age_group))
}
)

merged_df <- do.call(rbind, region_dfs)
merged_df$Region <- factor(merged_df$Region, levels = regions)
merged_df$Age_group <- factor(merged_df$Age_group, levels = age_group)

spectral_palette <- brewer.pal(9, "Reds")
color_gradient <- colorRampPalette(spectral_palette)(100)

pp <- ggplot(merged_df, aes(x = Age_group, y = Region, size = Allcounts, fill = TRGpct)) +
  geom_point(shape = 21, color = "black", stroke = 0.8) +
  scale_size_continuous(range = c(1, 12)) +
  scale_fill_gradientn(colors = color_gradient) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "",
       x = "Age",
       y = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#####————3K————#####
base_path <- "/PBA/SameRegion_across_DiffAge/"
reg_ref <- c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE")
filename <- "E104vs.E94.csv"

base_path <- "/PBA/SameRegion_across_DiffAge/"
reg_ref <- c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE")
filename <- "E104vs.E94.csv"

Deg_ls_E104 <- list(
  PFC = data.frame(),
  TMP = data.frame(),
  OCC = data.frame(),
  STR = data.frame(),
  THA = data.frame(),
  HIP = data.frame(),
  CERE = data.frame()
)

for (i in 1:length(reg_ref)) {
  reg <- reg_ref[i]
  filepath <- paste0(base_path, reg, "/", reg, "_", filename)
  df <- suppressMessages(read_csv(filepath))
  df <- subset(df, padj < 0.05 & abs(log2FoldChange) > log2(1.25))
  df <- df[df$log2FoldChange < 0, , drop = FALSE]
  df <- df[!grepl("^ENS", df$...1), , drop = FALSE]
  df$log2FoldChange <- -df$log2FoldChange
  df <- df[order(df$log2FoldChange, decreasing = TRUE), ]
  Deg_ls_E104[[reg]] <- df
}

Deg_ls_E94 <- list(
  PFC = data.frame(),
  TMP = data.frame(),
  OCC = data.frame(),
  STR = data.frame(),
  THA = data.frame(),
  HIP = data.frame(),
  CERE = data.frame()
)
for (i in 1:length(reg_ref)) {
  reg <- reg_ref[i]
  filepath <- paste0(base_path, reg, "/", reg, "_", filename)
  df <- suppressMessages(read_csv(filepath))
  df <- subset(df, padj < 0.05 & abs(log2FoldChange) > log2(1.25))
  df <- df[df$log2FoldChange > 0, , drop = FALSE]
  df <- df[!grepl("^ENS", df$...1), , drop = FALSE]
  df <- df[order(df$log2FoldChange, decreasing = TRUE), ]
  Deg_ls_E94[[reg]] <- df
}

Deg_ls_E104_sep <- list(
  PFC = Deg_ls_E104[["PFC"]]$...1[1:50],
  TMP = Deg_ls_E104[["TMP"]]$...1[1:50],
  OCC = Deg_ls_E104[["OCC"]]$...1[1:50],
  STR = Deg_ls_E104[["STR"]]$...1[1:50],
  THA = Deg_ls_E104[["THA"]]$...1[1:50],
  HIP = Deg_ls_E104[["HIP"]]$...1[1:50],
  CERE = Deg_ls_E104[["CERE"]]$...1[1:50]
)

Deg_ls_E94_sep <- list(
  PFC = Deg_ls_E94[["PFC"]]$...1[1:50],
  TMP = Deg_ls_E94[["TMP"]]$...1[1:50],
  OCC = Deg_ls_E94[["OCC"]]$...1[1:50],
  STR = Deg_ls_E94[["STR"]]$...1[1:50],
  THA = Deg_ls_E94[["THA"]]$...1[1:50],
  HIP = Deg_ls_E94[["HIP"]]$...1[1:50],
  CERE = Deg_ls_E94[["CERE"]]$...1[1:50]
)

top_genes_intersect <- NULL
for (reg in names(Deg_ls_E104)) {
  df <- Deg_ls_E104[[reg]]
  if (is.data.frame(df) && nrow(df) > 0) {
    top_genes_reg <- df$...1[1:min(200, nrow(df))]
    
    if (is.null(top_genes_intersect)) {
      top_genes_intersect <- top_genes_reg
    } else {
      top_genes_intersect <- intersect(top_genes_intersect, top_genes_reg)
    }
  }
}

print(top_genes_intersect)
# "MBP","SLC14A1","FA2H","TF","SERPINF2","ATP1A4"
wb_genes <- c("MBP","SLC14A1")

ctx_genes <- setdiff(
  setdiff(intersect(intersect(Deg_ls_E104$PFC$...1[1:floor(nrow(Deg_ls_E104$PFC)/4)],Deg_ls_E104$TMP$...1[1:floor(nrow(Deg_ls_E104$TMP)/4)]),Deg_ls_E104$OCC$...1[1:floor(nrow(Deg_ls_E104$OCC)/4)]),
          union(union(Deg_ls_E104$STR$...1,Deg_ls_E104$THA$...1),Deg_ls_E104$HIP$...1)),
  Deg_ls_E104$CERE$...1
)

scx_genes <- setdiff(
  setdiff(intersect(intersect(Deg_ls_E104$STR$...1[1:floor(nrow(Deg_ls_E104$PFC)/4)],Deg_ls_E104$THA$...1[1:floor(nrow(Deg_ls_E104$PFC)/4)]),Deg_ls_E104$HIP$...1[1:floor(nrow(Deg_ls_E104$PFC)/4)]),
          union(union(Deg_ls_E104$PFC$...1,Deg_ls_E104$TMP$...1),Deg_ls_E104$OCC$...1)),
  Deg_ls_E104$CERE$...1
)

cere_genes <- setdiff(Deg_ls_E104$CERE$...1[1:floor(nrow(Deg_ls_E104$CERE)/70)],
                      union(Deg_ls_E104$PFC$...1,
                            union(Deg_ls_E104$TMP$...1,
                                  union(Deg_ls_E104$OCC$...1,
                                        union(Deg_ls_E104$STR$...1,
                                              union(Deg_ls_E104$THA$...1,Deg_ls_E104$HIP$...1)
                                        )
                                  )
                            )
                      )
)
