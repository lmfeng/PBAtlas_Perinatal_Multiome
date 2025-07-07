source("/PBA/Function.R")
library(tidyverse)
library(reshape2)
#####————Development velocity————#####
###Gene
gexpr <- read.table("/PBA/geneCount_PBA.txt", header = TRUE, sep = "\t", row.names = 1)
gexpr <- t(gexpr)
meta <-  read.table("/PBA/sampleInfo_PBA.txt", header = TRUE, row.names = 1)
age_ref <- c("E76","E85","E94","E104","E109","P0","P3","P30")
reg_ref <- c("TMP","OCC","STR","THA","HIP","CERE")

gexpr_PFC <- gexpr[,rownames(meta[meta$Region == "PFC",])]
meta_PFC <- meta[meta$Region == "PFC",]

gexpr_PFC_mean <- data.frame(matrix(NA,nrow = nrow(gexpr_PFC)))
for (i in 1:length(age_ref)){
  age_idx <- age_ref[i]
  sample_idx <- rownames(meta_PFC[meta_PFC$Age == age_idx,])
  df_tmp <- gexpr_PFC[,sample_idx]
  row_means <- apply(df_tmp, 1, mean, na.rm = TRUE)
  gexpr_PFC_mean <- cbind(gexpr_PFC_mean,row_means)
}
gexpr_PFC_mean <- gexpr_PFC_mean[,-1]
rownames(gexpr_PFC_mean) <- rownames(gexpr_PFC)
colnames(gexpr_PFC_mean) <- age_ref

cor_PFC <- cor(gexpr_PFC_mean)
cor_PFC_slc <- cor_PFC[8,,drop = F]
max_cor <- max(cor_PFC_slc)
min_cor <- min(cor_PFC_slc)

cor_PFC_slc_scale <- (cor_PFC_slc-min_cor)/(max_cor-min_cor)

PFC_P30 <- gexpr_PFC_mean[,8,drop = F]
df_res <- data.frame(matrix(NA,ncol = 8,nrow = 6))
colnames(df_res) <- age_ref
rownames(df_res) <- reg_ref

for (i in 1:length(reg_ref)){
  reg_idx <- reg_ref[i]
  for (j in 1:length(age_ref)){
    age_idx <- age_ref[j]
    sample_idx <- rownames(meta[meta$Age == age_idx & meta$Region == reg_idx,])
    gexpr_tmp <- gexpr[,sample_idx,drop = F]
    if(ncol(gexpr_tmp) > 1){
      gexpr_tmp_mean <- as.data.frame(apply(gexpr_tmp, 1, mean, na.rm = TRUE))
    }else{
      gexpr_tmp_mean <- gexpr_tmp
    }
    colnames(gexpr_tmp_mean) <- age_idx
    cor_idx <- cor(PFC_P30,gexpr_tmp_mean)
    df_res[i,j] <- cor_idx
  }
  cat(reg_idx,"\n")
}

df_res_scale <- (df_res - min_cor)/(max_cor-min_cor)
rownames(cor_PFC_slc_scale) <- "PFC"
df_res <- rbind(cor_PFC_slc_scale,df_res_scale)

df_long <- df_res %>%
  rownames_to_column("Region") %>%
  pivot_longer(cols = -Region, 
               names_to = "Real_Age", 
               values_to = "Fit_Age")
df_long <- as.data.frame(df_long)
df_long$Region <- factor(df_long$Region, levels = c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE"))
df_long$Real_Age_lab <- df_long$Real_Age
df_long$Real_Age[df_long$Real_Age == "E76"] <- df_long$Fit_Age[df_long$Region == "PFC"][1]
df_long$Real_Age[df_long$Real_Age == "E85"] <- df_long$Fit_Age[df_long$Region == "PFC"][2]
df_long$Real_Age[df_long$Real_Age == "E94"] <- df_long$Fit_Age[df_long$Region == "PFC"][3]
df_long$Real_Age[df_long$Real_Age == "E104"] <- df_long$Fit_Age[df_long$Region == "PFC"][4]
df_long$Real_Age[df_long$Real_Age == "E109"] <- df_long$Fit_Age[df_long$Region == "PFC"][5]
df_long$Real_Age[df_long$Real_Age == "P0"] <- df_long$Fit_Age[df_long$Region == "PFC"][6]
df_long$Real_Age[df_long$Real_Age == "P3"] <- df_long$Fit_Age[df_long$Region == "PFC"][7]
df_long$Real_Age[df_long$Real_Age == "P30"] <- df_long$Fit_Age[df_long$Region == "PFC"][8]

df_long$Real_day[df_long$Real_Age_lab == "E76"] <- 76
df_long$Real_day[df_long$Real_Age_lab == "E85"] <- 85
df_long$Real_day[df_long$Real_Age_lab == "E94"] <- 94
df_long$Real_day[df_long$Real_Age_lab == "E104"] <- 104
df_long$Real_day[df_long$Real_Age_lab == "E109"] <- 109
df_long$Real_day[df_long$Real_Age_lab == "P0"] <- 115
df_long$Real_day[df_long$Real_Age_lab == "P3"] <- 118
df_long$Real_day[df_long$Real_Age_lab == "P30"] <- 145

df_long$Real_Age <- as.numeric(df_long$Real_Age)
df_long$Fit_Age <- as.numeric(df_long$Fit_Age)

reg_ref_all <- c("PFC", reg_ref)
res <- NULL
res <- data.frame(Slope = numeric(0), Intercept = numeric(0), Velocity = numeric(0))
for (reg in reg_ref_all) {
  df_region <- df_long[df_long$Region == reg, ]
  x <- df_region$Real_Age
  y <- df_region$Fit_Age
  
  model <- lm(y ~ x)
  
  slope <- coefficients(model)[2]
  intercept <- coefficients(model)[1]
  velocity <- 1 / slope
  
  res <- rbind(res, data.frame(Slope = slope, Intercept = intercept, Velocity = velocity))
}

rownames(res) <- reg_ref_all
res <- res[order(res$Velocity,decreasing = T),]
res

###Protein
load("/PBA/Proteins_originalData.RData")
expr = pexpr
meta = pmeta
age_ref <- c("E76","E85","E94","E104","E109","P0","P3","P30")
reg_ref <- c("TMP","OCC","STR","THA","HIP","CERE")

expr_PFC <- expr[,rownames(meta[meta$Region == "PFC",])]
meta_PFC <- meta[meta$Region == "PFC",]
expr_PFC_mean <- data.frame(matrix(NA,nrow = nrow(expr_PFC)))
for (i in 1:length(age_ref)){
  age_idx <- age_ref[i]
  sample_idx <- rownames(meta_PFC[meta_PFC$Age == age_idx,])
  df_tmp <- expr_PFC[,sample_idx]
  row_means <- apply(df_tmp, 1, mean, na.rm = TRUE)
  expr_PFC_mean <- cbind(expr_PFC_mean,row_means)
}
expr_PFC_mean <- expr_PFC_mean[,-1]
rownames(expr_PFC_mean) <- rownames(expr_PFC)
colnames(expr_PFC_mean) <- age_ref

cor_PFC <- cor(expr_PFC_mean)
cor_PFC_slc <- cor_PFC[8,,drop = F]
max_cor <- max(cor_PFC_slc)
min_cor <- min(cor_PFC_slc)

cor_PFC_slc_scale <- (cor_PFC_slc-min_cor)/(max_cor-min_cor)

PFC_P30 <- expr_PFC_mean[,8,drop = F]
df_res <- data.frame(matrix(NA,ncol = 8,nrow = 6))
colnames(df_res) <- age_ref
rownames(df_res) <- reg_ref

for (i in 1:length(reg_ref)){
  reg_idx <- reg_ref[i]
  for (j in 1:length(age_ref)){
    age_idx <- age_ref[j]
    sample_idx <- rownames(meta[meta$Age == age_idx & meta$Region == reg_idx,])
    expr_tmp <- expr[,sample_idx,drop = F]
    if(ncol(expr_tmp) > 1){
      expr_tmp_mean <- as.data.frame(apply(expr_tmp, 1, mean, na.rm = TRUE))
    }else{
      expr_tmp_mean <- expr_tmp
    }
    colnames(expr_tmp_mean) <- age_idx
    cor_idx <- cor(PFC_P30,expr_tmp_mean)
    df_res[i,j] <- cor_idx
  }
  cat(reg_idx,"\n")
}

df_res_scale <- (df_res - min_cor)/(max_cor-min_cor)
rownames(cor_PFC_slc_scale) <- "PFC"
df_res <- rbind(cor_PFC_slc_scale,df_res_scale)

df_long <- df_res %>%
  rownames_to_column("Region") %>%
  pivot_longer(cols = -Region, 
               names_to = "Real_Age", 
               values_to = "Fit_Age")
df_long <- as.data.frame(df_long)
df_long$Region <- factor(df_long$Region, levels = c("PFC", "TMP", "OCC", "STR", "THA", "HIP", "CERE"))
df_long$Real_Age_lab <- df_long$Real_Age
df_long$Real_Age[df_long$Real_Age == "E76"] <- df_long$Fit_Age[df_long$Region == "PFC"][1]
df_long$Real_Age[df_long$Real_Age == "E85"] <- df_long$Fit_Age[df_long$Region == "PFC"][2]
df_long$Real_Age[df_long$Real_Age == "E94"] <- df_long$Fit_Age[df_long$Region == "PFC"][3]
df_long$Real_Age[df_long$Real_Age == "E104"] <- df_long$Fit_Age[df_long$Region == "PFC"][4]
df_long$Real_Age[df_long$Real_Age == "E109"] <- df_long$Fit_Age[df_long$Region == "PFC"][5]
df_long$Real_Age[df_long$Real_Age == "P0"] <- df_long$Fit_Age[df_long$Region == "PFC"][6]
df_long$Real_Age[df_long$Real_Age == "P3"] <- df_long$Fit_Age[df_long$Region == "PFC"][7]
df_long$Real_Age[df_long$Real_Age == "P30"] <- df_long$Fit_Age[df_long$Region == "PFC"][8]

df_long$Real_day[df_long$Real_Age_lab == "E76"] <- 76
df_long$Real_day[df_long$Real_Age_lab == "E85"] <- 85
df_long$Real_day[df_long$Real_Age_lab == "E94"] <- 94
df_long$Real_day[df_long$Real_Age_lab == "E104"] <- 104
df_long$Real_day[df_long$Real_Age_lab == "E109"] <- 109
df_long$Real_day[df_long$Real_Age_lab == "P0"] <- 115
df_long$Real_day[df_long$Real_Age_lab == "P3"] <- 118
df_long$Real_day[df_long$Real_Age_lab == "P30"] <- 145

df_long$Real_Age <- as.numeric(df_long$Real_Age)
df_long$Fit_Age <- as.numeric(df_long$Fit_Age)

reg_ref_all <- c("PFC", reg_ref)
res <- NULL
for (reg in reg_ref_all) {
  df_region <- df_long[df_long$Region == reg, ]
  x <- df_region$Real_Age
  y <- df_region$Fit_Age
  
  model <- lm(y ~ x)
  
  slope <- coefficients(model)[2]
  intercept <- coefficients(model)[1] 
  velocity <- 1 / slope
  
  res <- rbind(res, data.frame(Slope = slope, Intercept = intercept, Velocity = velocity))
}

rownames(res) <- reg_ref_all
res <- res[order(res$Velocity,decreasing = T),]
res

#####————DVG and DVP————#####
###DVG
results_list <- list()
for (region in reg_ref) {

  gexpr_region <- gexpr[, rownames(meta[meta$Region == region, ])]
  meta_region <- meta[meta$Region == region, ]

  gexpr_region_mean <- data.frame(matrix(NA, nrow = nrow(gexpr_region)))
  for (i in 1:length(age_ref)) {
    age_idx <- age_ref[i]
    sample_idx <- rownames(meta_region[meta_region$Age == age_idx, ])
    df_tmp <- gexpr_region[, sample_idx,drop = F]
    if (ncol(df_tmp) > 1){
      row_means <- as.data.frame(apply(df_tmp, 1, mean, na.rm = TRUE))
    }else{
      row_means <- df_tmp
    }
    gexpr_region_mean <- cbind(gexpr_region_mean, row_means)
  }
  gexpr_region_mean <- gexpr_region_mean[, -1]
  
  results_region <- data.frame(GeneName = character(), correlation = numeric(), Region = character(), stringsAsFactors = FALSE)

  for (gene in rownames(gexpr_region_mean)) {
    test <- gexpr_region_mean[gene, ]
    
    test_long <- test %>%
      rownames_to_column("GeneName") %>%
      pivot_longer(cols = -GeneName, names_to = "Real_Age", values_to = "y")
    test_long <- as.data.frame(test_long)
    
    test_long$Real_Age <- factor(test_long$Real_Age, levels = age_ref)
    test_long$x <- df_long$Real_Age[1:8]
    
    correlation_value <- cor(test_long$y, df_long$Fit_Age[df_long$Region == region])
    results_region <- rbind(results_region, data.frame(GeneName = gene, correlation = correlation_value, Region = region))
  }
  results_region <- results_region %>%
    filter(!is.na(correlation))
  
  results_list[[region]] <- results_region
  cat(region,"\n")
}

PFC_dvg_slc <- results_list$PFC %>%
  dplyr::filter(correlation > 0.7)
TMP_dvg_slc <- results_list$TMP %>%
  filter(correlation > 0.7)
OCC_dvg_slc <- results_list$OCC %>%
  filter(correlation > 0.7)
STR_dvg_slc <- results_list$STR %>%
  filter(correlation > 0.7)
THA_dvg_slc <- results_list$THA %>%
  filter(correlation > 0.7)
HIP_dvg_slc <- results_list$HIP %>%
  filter(correlation > 0.7)
CERE_dvg_slc <- results_list$CERE %>%
  filter(correlation > 0.7)

###DVP
results_list <- list()
for (region in reg_ref) {
  expr_region <- expr[, rownames(meta[meta$Region == region, ])]
  meta_region <- meta[meta$Region == region, ]
  expr_region_mean <- data.frame(matrix(NA, nrow = nrow(expr_region)))
  for (i in 1:length(age_ref)) {
    age_idx <- age_ref[i]
    sample_idx <- rownames(meta_region[meta_region$Age == age_idx, ])
    df_tmp <- expr_region[, sample_idx,drop = F]
    if (ncol(df_tmp) > 1){
      row_means <- as.data.frame(apply(df_tmp, 1, mean, na.rm = TRUE))
    }else{
      row_means <- df_tmp
    }
    expr_region_mean <- cbind(expr_region_mean, row_means)
  }
  expr_region_mean <- expr_region_mean[, -1]

  results_region <- data.frame(GeneName = character(), correlation = numeric(), Region = character(), stringsAsFactors = FALSE)

  for (gene in rownames(expr_region_mean)) {
    test <- expr_region_mean[gene, ]
    
    test_long <- test %>%
      rownames_to_column("GeneName") %>%
      pivot_longer(cols = -GeneName, names_to = "Real_Age", values_to = "y")
    test_long <- as.data.frame(test_long)
    
    test_long$Real_Age <- factor(test_long$Real_Age, levels = age_ref)
    test_long$x <- df_long$Real_Age[1:8]
    

    correlation_value <- cor(test_long$y, df_long$Fit_Age[df_long$Region == region])
    results_region <- rbind(results_region, data.frame(GeneName = gene, correlation = correlation_value, Region = region))
  }
  results_region <- results_region %>%
    filter(!is.na(correlation))
  
  results_list[[region]] <- results_region
  cat(region,"\n")
}

PFC_dvp_slc <- results_list$PFC %>%
  dplyr::filter(correlation > 0.7)
TMP_dvp_slc <- results_list$TMP %>%
  filter(correlation > 0.7)
OCC_dvp_slc <- results_list$OCC %>%
  filter(correlation > 0.7)
STR_dvp_slc <- results_list$STR %>%
  filter(correlation > 0.7)
THA_dvp_slc <- results_list$THA %>%
  filter(correlation > 0.7)
HIP_dvp_slc <- results_list$HIP %>%
  filter(correlation > 0.7)
CERE_dvp_slc <- results_list$CERE %>%
  filter(correlation > 0.7)

#####————DVG classification————#####
pfc_spec <- read.table("/PBA/PFC/pfc_hum.spec.xls", sep = "\t", header = TRUE)
tmp_spec <- read.table("/PBA/TMP/tmp_hum.spec.xls", sep = "\t", header = TRUE)
occ_spec <- read.table("/PBA/OCC/occ_hum.spec.xls", sep = "\t", header = TRUE)
str_spec <- read.table("/PBA/STR/str_hum.spec.xls", sep = "\t", header = TRUE)
tha_spec <- read.table("/PBA/THA/tha_hum.spec.xls", sep = "\t", header = TRUE)
str_spec <- read.table("/PBA/STR/str_hum.spec.xls", sep = "\t", header = TRUE)
hip_spec <- read.table("/PBA/HIP/hip_hum.spec.xls", sep = "\t", header = TRUE)
cere_spec <- read.table("/PBA/CERE/cere_hum.spec.xls", sep = "\t", header = TRUE)

dvg_list <- list(
  PFC = PFC_dvg_slc$GeneName,
  TMP = TMP_dvg_slc$GeneName,
  OCC = OCC_dvg_slc$GeneName,
  STR = STR_dvg_slc$GeneName,
  THA = THA_dvg_slc$GeneName,
  HIP = HIP_dvg_slc$GeneName,
  CERE = CERE_dvg_slc$GeneName
)

spec_list <- list(
  PFC = pfc_spec,
  TMP = tmp_spec,
  OCC = occ_spec,
  STR = str_spec,
  THA = tha_spec,
  HIP = hip_spec,
  CERE = cere_spec
)

results <- lapply(names(dvg_list), function(region) {
  dvg_genes <- dvg_list[[region]]
  spec_df <- spec_list[[region]]
  
  hcdvg <- get_ct_related_hcdvg(region, dvg_genes)
  classifications <- sapply(dvg_genes, function(g) classify_ct_related_dvg(g, region, hcdvg, spec_df))
  
  all_categories <- c("High-confidence cell-type related",
                      "Low-confidence cell-type related",
                      "Non-cell-type related",
                      "Unknown")
  counts <- table(factor(classifications, levels = all_categories))
  counts[is.na(counts)] <- 0
  
  data.frame(
    Region = region,
    Category = names(counts),
    Count = as.numeric(counts),
    Percentage = as.numeric(counts / length(dvg_genes) * 100))
}) %>% bind_rows()


results$Region <- factor(results$Region,levels = c("PFC","TMP","OCC","STR","THA","HIP","CERE"))
results$Category <- factor(results$Category,levels = c("Low-confidence cell-type related","High-confidence cell-type related","Non-cell-type related","Unknown"))
pdf("/PBA/dvg_composition.bar.pdf",width = 10,height = 7)
ggplot(results, aes(x = Region, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF")) +
  labs(x = "Brain Region", y = "Percentage (%)", fill = "DVG Category") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(nrow = 2))
dev.off()
