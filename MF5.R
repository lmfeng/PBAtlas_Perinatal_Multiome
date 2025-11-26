library(tidyverse)
library(WGCNA)
library(fastcluster)
library(clusterProfiler)
library(openxlsx)
library(dlookr) 

#####—————integrate gene and protein—————#####
set.seed(100)

age_ref <- c("E76","E85","E94","E104","E109","P0","P3","P30")
reg_ref <- c("PFC","TMP","OCC","STR","THA","HIP","CERE")

load("/home/wangzm/Project/PBA/basic_analysis/MF5_related/PBAtlas.geneExpr.norm.RData")
gexpr.g = gexpr
meta.g = meta

load("/home/wangzm/Project/PBA/basic_analysis/MF5_related/PBAtlas.proteinExpr.norm.RData")
gexpr.p = combat_edata1
meta.p = meta

sampList = intersect(colnames(gexpr.g), colnames(gexpr.p))
x01 = match(sampList, colnames(gexpr.g))
x02 = match(sampList, colnames(gexpr.p))
gexpr2.g = gexpr.g[,x01]
gexpr2.p = gexpr.p[,x02]
meta2.g = meta.g[x01, ]
meta2.p = meta.p[x02, ]

genes.pc = intersect(rownames(gexpr2.g), rownames(gexpr2.p))
gexpr2.g = gexpr2.g[genes.pc, ]
gexpr2.p = gexpr2.p[genes.pc, ]
gexpr2.p <- as.matrix(gexpr2.p)
gexpr2.g

bc_trans<- function(x){
  x = 2^x
  bc_res = boxcox(x ~ 1, plotit = F)
  best_lambda = bc_res$x[which.max(bc_res$y)]
  y = (x^best_lambda - 1)/best_lambda
  z = sign(y) * log2(abs(y))
  return(z)
}
gexpr2.p = apply(gexpr2.p, 2, bc_trans)

for(i in 1:nrow(gexpr2.p)){
  eachline = gexpr2.p[i,]
  df = data.frame(sample = names(eachline), value = eachline)
  xx= imputate_outlier(df, value, method = "capping", cap_ntiles = c(0.01, 0.95))
  outlier_num = length(attr(xx, "outliers"))
  ##
  if(outlier_num > 0 ){
    gexpr2.p[i,] = as.numeric(xx)
  }  
}

gexpr2.g = base::scale(gexpr2.g, center = T, scale = T)
gexpr2.p = base::scale(gexpr2.p, center = T, scale = T)
boxplot(gexpr2.g, outline = F)
boxplot(gexpr2.p, outline = F)

##merge data
rownames(gexpr2.g) = paste0("g.",rownames(gexpr2.g))
rownames(gexpr2.p) = paste0("p.",rownames(gexpr2.p))
gexpr <- NULL
gexpr = rbind(gexpr2.g, gexpr2.p)
geneList = rownames(gexpr)
meta = cbind(meta2.g, meta2.p)

ageList = meta$Age
dayList = ageList
dayList = rep(0,length(ageList))
dayList[which(ageList == "E76")] = 76
dayList[which(ageList == "E85")] = 85
dayList[which(ageList == "E94")] = 94
dayList[which(ageList == "E104")] = 104
dayList[which(ageList == "E109")] = 109
dayList[which(ageList == "P0")] = 115
dayList[which(ageList == "P3")] = 118
dayList[which(ageList == "P30")] = 145
meta$Day = dayList

save(gexpr, meta, file="/PBA/pg.WGCNA.RData")

#####—————WGCNA—————#####
set.seed(100)
load("/PBA/pg.WGCNA.RData")

recGenexp = t(gexpr)
powers = c(c(1:15), seq(from = 15, to=30, by=2))
sft = pickSoftThreshold(recGenexp, powerVector = powers, verbose = 3)

recGenexp = t(gexpr)
meta = meta

geneNet = blockwiseModules(recGenexp, power=6, minModuleSize=10,
                           mergeCutHeight=0.2,networkType = "signed",
                           TOMType = "signed",numericLabels=T,pamRespectsDendro=F,
                           saveTOMs=F,verbose=0)
mergedColors = labels2colors(geneNet$colors)
plotDendroAndColors(geneNet$dendrograms[[1]],mergedColors[geneNet$blockGenes[[1]]],
                    "Module colors", addTextGuide=T, dendroLabels = F, hang = 0.03,addGuide = TRUE,
                    guideHang = 0.05,marAll = c(1, 5, 3,1))

save(gexpr, geneNet, meta, file="/PBA/pg.WGCNA.gene-net-signed.RData")

load("/PBA/pg.WGCNA.gene-net-signed.RData")

eigenexp = geneNet$MEs
moduleList = colnames(eigenexp)
gexpr = t(gexpr)
geneList = colnames(gexpr)
colorList = geneNet$colors

outFile = "/PBA/pg.WGCNA.geneModule.reassigned.xls"
if(file.exists(outFile)) file.remove(outFile)

res <- c()
for(i in 0:max(colorList)){
  myindex = which(colorList == i)
  clusterGene = as.character(geneList[myindex])
  cat(i, "\n")
  
  for(j in 1:length(clusterGene)){
    geneName = as.character(clusterGene[j])
    myindex2 = which(geneList == geneName)
    oneGenexp = as.numeric(gexpr[,myindex2])
    oneGenecor = cor(oneGenexp, eigenexp,method="p",use = 'pairwise.complete.obs')

    moduleName.old = paste("ME",i,sep="")
    myindex3 = which(moduleList == moduleName.old)
    
    moduleCor.old = oneGenecor[myindex3]
    moduleCor.new = max(oneGenecor)

    output = c(geneName, moduleName.old)
    if((moduleCor.new -moduleCor.old) >0.3){
      maxIndex = which(oneGenecor == max(oneGenecor))
      output = c(geneName, moduleList[maxIndex])
    }	
    
    ##recording	
    if(length(res) == 0){
      res = output
    }
    else{
      res = rbind(res,output)
    }
  }
}	

for(i in 0:max(colorList)){
  moduleName = paste("ME",i,sep="")
  myindex = which(as.character(res[,2]) == moduleName)
  clusterGene = as.character(res[myindex,1])
  moduleSize = length(clusterGene)
  geneNum = length(grep("^g", clusterGene))
  proteinNum = length(grep("^p", clusterGene))
  cat(moduleName,moduleSize, geneNum, proteinNum, "\n",sep="\t")
  cat(moduleName,clusterGene,"\n",file = outFile, sep="\t",append=T)
}

annotation_colors <- list(
  Age = c(
    "E76" = "#f3993a",
    "E85" = "#ffee6f",
    "E94" = "#add5a2",
    "E104" = "#7a7b78",
    "E109" = "#9933cc",
    "P0" = "#0066ff",
    "P3" = "#33cccc",
    "P30" = "#ff66cc"),
  Region=c(
    "PFC" = "#d71345",
    "TMP" = "#f26522",
    "OCC" = "#f7acbc",
    "STR" = "#6950a1",
    "THA" = "#009ad6",
    "HIP" = "#3CB371",
    "CERE" = "#8B4513"
  )
)
set.seed(100)

ageRef = c("E76","E85","E94", "E104", "E109","P0","P3","P30")	
dayRef = c(76, 85, 94, 104, 109, 115, 118, 145)	
regionRef = c("PFC","TMP","OCC","STR","THA","HIP","CERE")
dayList = meta$Day
ageList = meta$Age

gexpr= t(gexpr)
sampList = rownames(gexpr)
geneList = colnames(gexpr)

modData = readLines("/PBA/pg.WGCNA.geneModule.reassigned.xls")
resEigen <- c()
smEigen <- c()
recPlts <- list()

##
for(i in 1:length(modData)){
  rec = unlist(strsplit(as.character(modData[i]),split="\t",fixed=T))
  modName = rec[1]
  clustGene = rec[2:length(rec)]

  moduleSize = length(clustGene)
  geneNum = length(grep("^g", clustGene))
  proteinNum = length(grep("^p", clustGene))
  cat(modName,moduleSize, geneNum, proteinNum, "\n",sep="\t")
  eachName = paste(modName,moduleSize, geneNum, proteinNum, sep = "_")

  idx = match(clustGene,geneList)
  moduleGenexp = gexpr[,idx]

  myEigen = moduleEigengenes(moduleGenexp,colors=rep("red",ncol(moduleGenexp)))
  geneigen = unlist(myEigen$eigengenes)

  if(length(resEigen) == 0){
    resEigen = geneigen
  }	else{
    resEigen = rbind(resEigen,geneigen)
  }	

  res <- NULL
  res = data.frame(eigen.val = as.numeric(geneigen), Age = ageList, Day = dayList,
                   Region = meta$Region, Sex = meta$Sex)
  res$Region = factor(res$Region, regionRef)

  p0 <- ggplot(data=res, aes(x = Day, y = eigen.val)) +
    labs(title = eachName, x = "", y = "Eigengene") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 16),
          panel.background = element_rect(fill = "white", color = "white"),
          panel.grid.major = element_line(color = "grey80"), 
          panel.grid.minor = element_line(color = "grey90"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
    geom_point(aes(x = Day, y = eigen.val, colour = Region)) +
    scale_x_continuous(breaks = dayRef, labels = ageRef) +
    scale_colour_manual(name = "Region",values = annotation_colors$Region) +
    geom_smooth(level = 0, span = 0.75, aes( group = Region,colour = Region,fill = Region,alpha=0.5),size=1.5) +
    scale_fill_manual(name = "Region",values = annotation_colors$Region) +
    geom_vline(xintercept = 115,colour="grey20",linewidth =1,linetype = "dashed")

  recPlts[[modName]] = p0

  smooth_values <- ggplot_build(p0)$data[[1]]
  if(length(smEigen) == 0){
    smEigen = smooth_values$y
  }	else{
    smEigen = rbind(smEigen,smooth_values$y)
  }	
}

colnames(resEigen) = sampList
rownames(resEigen) = paste("ME",0:(length(modData)-1),sep="")

colnames(smEigen) = sampList
rownames(smEigen) = paste("ME",0:(length(modData)-1),sep="")

save(gexpr, meta, geneNet, resEigen, smEigen,recPlts, file="/PBA/pg.WGCNA.eigengenes.RData")

#####—————Traits—————#####
load(file="/PBA/pg.WGCNA.eigengenes.RData")
modList = rownames(resEigen)
dayList = meta$Day
regionRef = c("PFC","TMP","OCC", "STR","THA", "HIP", "CERE")
modData = readLines("/PBA/pg.WGCNA.geneModule.reassigned.xls")

pctMat <- NULL
for(i in 1:length(modData)){
  rec = unlist(strsplit(as.character(modData[i]),split="\t",fixed=T))
  modName = rec[1]
  clustGene = rec[2:length(rec)]
  moduleSize = length(clustGene)

  eachout <- NULL

  genes = substring(clustGene[grep("^g", clustGene)],3)
  proteins = substring(clustGene[grep("^p", clustGene)],3)
  eachout = c(length(genes)/moduleSize , length(proteins)/moduleSize)
  if(length(pctMat) == 0){
    pctMat = eachout
  }else{
    pctMat = rbind(pctMat, eachout)
  }
}
colnames(pctMat) = c("Gene", "Protein")
rownames(pctMat)  = modList

pctMat_long <- as.data.frame(pctMat) %>%
  rownames_to_column("module") %>%
  pivot_longer(cols = c(Gene, Protein), 
               names_to = "Type", 
               values_to = "Percentage")

pctMat_long$Type <- factor(pctMat_long$Type,levels = c("Gene","Protein"))
pctMat_long$module <- factor(pctMat_long$module,levels = rownames(pctMat))

ggplot(pctMat_long, aes(x = module, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#F5F5DC","#F4A460")) +
  theme_classic() + 
  labs(y = "Region", x = "Age") +
  theme(strip.background = element_rect(color = "black", fill="black",),
        strip.text = element_text(size = 12,color = "white"),
        legend.position = "left",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
  )

cls = hclust(dist(smEigen), method = "ave")
dend <- as.dendrogram(cls)
kMax = 11
plot(dend)

myTree = fastcluster::hclust(dist(smEigen),method="ave")
groupMods = cutree(myTree, 11)

mat.st <- NULL

for(i in 1:nrow(resEigen)){
  modName = rownames(resEigen)[i]

  df = data.frame(Region = meta$Region, Age = meta$Age, Day = dayList, value = resEigen[i,], value.sm = smEigen[i, ])
  mod.t = lm(value.sm ~ Day, df)
  mod.s = lm(value.sm ~ Region, df)

  all.pval.t = (-1) * log10(overall_p(mod.t))
  all.pval.s = (-1) * log10(overall_p(mod.s))

  cor.t = cor(df$value.sm, df$Day)
  sep.pval.s = (-1) * log10(summary(mod.s)$coefficients[,4])
  names(sep.pval.s)[1] = "CERE"
  names(sep.pval.s) = gsub("Region", "", names(sep.pval.s))
  sep.pval.s = sep.pval.s[regionRef]
  sep.pval.s = t(sep.pval.s)

  eachout = data.frame(mods = modName, all.pval.t = all.pval.t, cor.t = cor.t,
                       all.pval.s = all.pval.s, sep.pval.s)
  
  if(length(mat.st) == 0){
    mat.st = eachout
  }else{
    mat.st = rbind(mat.st, eachout)
  }
}
rownames(mat.st) = mat.st$mods
mat.st = as.matrix(mat.st[, -1])

spTest <- NULL
for(i in 1:nrow(mat.st)){
  rec = mat.st[i,]
  eachout = rep(0, length(rec))

  if(rec[1] > 2){
    eachout[1] = "t"
    if(rec[2] > 0.25) eachout[2] = "plus"
    if(rec[2] < (-0.25)) eachout[2] = "minus" 
  }

  if(rec[3] > 2){
    eachout[3] = "s"
    rem = rec[4:length(rec)]
    idx.tops = findTops(rem)
    eachout[3 + idx.tops] = "spatial"
  }

  if(length(spTest) == 0){
    spTest = eachout
  }else{
    spTest = rbind(spTest, eachout)
  }
}
colnames(spTest) = colnames(mat.st)
rownames(spTest) = modList
overall_means <- rowMeans(smEigen)

Spatial2 <- matrix(0, 
                   nrow = nrow(smEigen), 
                   ncol = length(regionRef),
                   dimnames = list(rownames(smEigen), regionRef))

for(mod in rownames(smEigen)) {
  for(region in regionRef) {
    region_values <- smEigen[mod, meta$Region == region]
    region_mean <- mean(region_values)
    if(region_mean > overall_means[mod]) {
      Spatial2[mod, region] <- 1
    } else {
      Spatial2[mod, region] <- -1
    }
  }
}

for(mod in rownames(spTest)) {
  for(region in regionRef) {
    if(spTest[mod, region] != "spatial") {
      Spatial2[mod, region] <- 0 
    }
  }
}
Spatial2

pfc_ct_path <- c("/PBA/PFC/pfc_hum.spec.xls")
tmp_ct_path <- c("/PBA/TMP/tmp_hum.spec.xls")
occ_ct_path <- c("/PBA/OCC/occ_hum.spec.xls")
str_ct_path <- c("/PBA/STR/str_hum.spec.xls")
tha_ct_path <- c("/PBA/THA/tha_hum.spec.xls")
hip_ct_path <- c("/PBA/HIP/hip_hum.spec.xls")
cere_ct_path <- c("/PBA/CERE/cere_hum.spec.xls")

pfc_ct_mar <- extract_spectop(pfc_ct_path)
tmp_ct_mar <- extract_spectop(tmp_ct_path)
occ_ct_mar <- extract_spectop(occ_ct_path)
str_ct_mar <- extract_spectop(str_ct_path)
tha_ct_mar <- extract_spectop(tha_ct_path)
hip_ct_mar <- extract_spectop(hip_ct_path)
cere_ct_mar <- extract_spectop(cere_ct_path)

pfc_ct_mat <- get_ct_mat(pfc_ct_mar)
tmp_ct_mat <- get_ct_mat(tmp_ct_mar)
occ_ct_mat <- get_ct_mat(occ_ct_mar)
str_ct_mat <- get_ct_mat(str_ct_mar)
tha_ct_mat <- get_ct_mat(tha_ct_mar)
hip_ct_mat <- get_ct_mat(hip_ct_mar)
cere_ct_mat <- get_ct_mat(cere_ct_mar)

ct_ref <- Reduce(union, list(
  unique(pfc_ct_mat$ctID),
  unique(tmp_ct_mat$ctID),
  unique(occ_ct_mat$ctID),
  unique(str_ct_mat$ctID),
  unique(tha_ct_mat$ctID),
  unique(hip_ct_mat$ctID),
  unique(cere_ct_mat$ctID)
))

ct_ref_mat <- rbind(pfc_ct_mat,tmp_ct_mat,occ_ct_mat,str_ct_mat,tha_ct_mat,hip_ct_mat,cere_ct_mat)
rownames(ct_ref_mat) <- NULL
ct_ref_mat <- unique(ct_ref_mat)
add_rows <- data.frame(     
  ctID = rep("MSN",times = 2),
  geneID = c("DRD1","DRD2")
)
ct_ref_mat <- rbind(ct_ref_mat,add_rows)

recMarkers.lister = readRDS("/PBA/Lister_Cell_2022.perinatal.celltypeSpec.gene.rds")
##CERE data
df.pc = read.csv("/PBA/01MinusPC.tsv", header = T, sep = "\t")
df.gn = read.csv("/PBA/04MinusGN.tsv", header = T, sep = "\t")
recMarkers.millen = data.frame(cluster = c(rep("PC", nrow(df.pc)), rep("GN", nrow(df.gn))),
                               gene = c(df.pc$Gene, df.gn$Gene))
##STR data
recMarkers.msn = readRDS("/PBA/Stauffer_CurrentBiology2021.celltypeSpec.gene.rds")
recMarkers.msn = subset(recMarkers.msn, cluster %in% c("DRD1", "DRD2"))

##DG
load("/PBA/DG_exn_mar.RData") #gc_mar
dg_df <- data.frame(
  ctID = "DG.ExN",
  geneID = gc_mar$gene          
)

##
recMarkers = rbind(recMarkers.lister[,c(6,7)], recMarkers.millen, recMarkers.msn[, c(6,7)])
celltypeSpec = data.frame(pathID = recMarkers$cluster, 
                          geneID = recMarkers$gene)

celltypeSpec$pathID <- as.character(celltypeSpec$pathID)
celltypeSpec$pathID[celltypeSpec$pathID %in% c("Micro","Astro","OPC","Oligo")] <- "Glia"
celltypeSpec$pathID[celltypeSpec$pathID %in% c("DRD2","DRD1")] <- "MSN"
celltypeSpec$pathID[celltypeSpec$pathID == "PC"] <- "Purkinje.cell"
celltypeSpec$pathID[celltypeSpec$pathID == "GN"] <- "CB.granule.cell"

celltypeSpec <- subset(celltypeSpec,celltypeSpec$pathID %in% unique(ct_ref_mat$ctID))
colnames(celltypeSpec) <- c("ctID","geneID")

celltypeSpec <- rbind(celltypeSpec, dg_df)
ct_ref_mat <- rbind(ct_ref_mat,celltypeSpec)
table(ct_ref_mat$ctID)

celltypeEnrich <- array(0, dim = c(nrow(mat.st), length(ct_ref)))
rownames(celltypeEnrich) = rownames(mat.st)
colnames(celltypeEnrich) = ct_ref
library(clusterProfiler)
for(i in 1:length(modData)){
  rec = unlist(strsplit(as.character(modData[i]),split="\t",fixed=T))
  modName = rec[1]
  clustGene = unique(substring(rec[2:length(rec)],3))
  miTest =  enricher(clustGene, TERM2GENE = ct_ref_mat, maxGSSize = 8000, qvalueCutoff = 0.05)
  miTest.df = as.data.frame(miTest)
  if(nrow(miTest.df) > 0){
    cat(modName, length(clustGene), miTest.df$ID, "\n", sep = "\t")
    idx.sort = match(miTest.df$ID, ct_ref)
    celltypeEnrich[i, idx.sort] = 1
  }
}

transMat <- NULL
for(i in 1:length(modData)){
  rec = unlist(strsplit(as.character(modData[i]),split="\t",fixed=T))
  modName = rec[1]
  clustGene = rec[2:length(rec)]
  moduleSize = length(clustGene)
  
  ###
  eachout <- NULL
  ###
  genes = substring(clustGene[grep("^g", clustGene)],3)
  proteins = substring(clustGene[grep("^p", clustGene)],3)
  op = intersect(genes, proteins)

  idx = which(groupMods == groupMods[i])
  clustGene <- NULL
  for(k in idx){
    rec = unlist(strsplit(as.character(modData[k]),split="\t",fixed=T))
    clustGene = c(clustGene, rec[2:length(rec)])
  }
  moduleSize = length(clustGene)
  
  ###
  proteins = substring(clustGene[grep("^p", clustGene)],3)
  op2 = intersect(genes, proteins)
  
  ##
  geneNum = length(genes)+1e-6
  eachout = c(length(op)/geneNum , length(op2)/geneNum)
  if(length(transMat) == 0){
    transMat = eachout
  }else{
    transMat = rbind(transMat, eachout)
  }
}

rownames(transMat) = modList
colnames(transMat) = c("Percent.trans.mod", "Percent.trans.grp")

tfAnnot.guo = read.table("/PBA/TF-Target-information.txt", header = T, sep = "\t")
load("/PBA/tftargets.rda")

ln = readLines("/PBA/c3.tft.v2024.1.Hs.symbols.gmt")

gsea <- list()
for(i in 1:length(ln)){
  rec = unlist(strsplit(as.character(ln[i]), "\t"))
  targets = rec[-(1:2)]
  rem = unlist(strsplit(as.character(rec[1]), "_"))
  tfname = rem[1]
  if(length(gsea[[tfname]]) == 0){
    gsea[[tfname]] = rec
  }else{
    gsea[[tfname]] = unique(c(gsea[[tfname]], rec))
  }
}

tf.tot = unique(c(unique(tfAnnot.guo$TF), names(ITFP), names(Marbach2016), names(TRRUST), names(gsea)))
tf2target <- list()
for(i in 1:length(tf.tot)){
  eachtf = tf.tot[i]
  eachtarget = unique(c(ITFP[[eachtf]], Marbach2016[[eachtf]], TRRUST[[eachtf]], gsea[[eachtf]]))

  xi = which(tfAnnot.guo$TF == eachtf)
  eachtarget = unique(c(eachtarget, tfAnnot.guo$target[xi]))

  tf2target[[eachtf]] = eachtarget
}

regMat <- NULL
for(i in 1:length(modData)){
  rec = unlist(strsplit(as.character(modData[i]),split="\t",fixed=T))
  modName = rec[1]
  clustGene = rec[2:length(rec)]
  moduleSize = length(clustGene)

  eachout <- NULL
  genes = substring(clustGene[grep("^g", clustGene)],3)
  proteins = substring(clustGene[grep("^p", clustGene)],3)
  tfs = intersect(proteins, names(tf2target))
  
  targets = unique(unlist(tf2target[proteins]))
  op2 = intersect(genes, targets)

  idx = which(groupMods == groupMods[i])
  clustGene <- NULL
  for(k in idx){
    rec = unlist(strsplit(as.character(modData[k]),split="\t",fixed=T))
    clustGene = c(clustGene, rec[2:length(rec)])
  }
  moduleSize = length(clustGene)

  proteins = substring(clustGene[grep("^p", clustGene)],3)
  tfs = intersect(proteins, names(tf2target))

  targets = unique(unlist(tf2target[proteins]))
  op3 = intersect(genes, targets)

  geneNum = length(genes)+1e-6
  eachout = c(length(op2)/geneNum, length(op3)/geneNum)
  if(length(regMat) == 0){
    regMat = eachout
  }else{
    regMat = rbind(regMat, eachout)
  }
  
}
rownames(regMat) = modList
colnames(regMat) = c("Percent.reg.mod", "Percent.reg.grp")

proteinAnnot = read.table("/PBA/9606.protein.info.v11.0.txt", header = T, sep = "\t")
proteinLink = read.table("/PBA/9606.protein.links.v11.0.txt", header = T, sep = " ")
##
idx01 = match(proteinLink[,1], proteinAnnot[,1])
idx02 = match(proteinLink[,2], proteinAnnot[,1])
stringdb = cbind(proteinAnnot[idx01, 2], proteinAnnot[idx02, 2])

##
intMat <- NULL
for(i in 1:length(modData)){
  rec = unlist(strsplit(as.character(modData[i]),split="\t",fixed=T))
  modName = rec[1]
  clustGene = rec[2:length(rec)]
  moduleSize = length(clustGene)

  eachout <- NULL
  genes = substring(clustGene[grep("^g", clustGene)],3)
  clustGene = substring(clustGene,3)
  idx01 = which(stringdb[,1] %in% genes)
  idx02 = which(stringdb[,2] %in% clustGene)
  op5 = length(intersect(idx01, idx02))
  
  idx = which(groupMods == groupMods[i])
  clustGene <- NULL
  for(k in idx){
    rec = unlist(strsplit(as.character(modData[k]),split="\t",fixed=T))
    clustGene = c(clustGene, rec[2:length(rec)])
  }
  clustGene = substring(clustGene,3)
  idx02 = which(stringdb[,2] %in% clustGene)
  op6 = length(intersect(idx01, idx02))
  
  
  ##
  geneNum = length(genes)+1e-6
  eachout = c(op5/geneNum, op6/geneNum)
  if(length(intMat) == 0){
    intMat = eachout
  }else{
    intMat = rbind(intMat, eachout)
  }
  
}
rownames(intMat) = modList
colnames(intMat) = c("Norm.int.mod", "Norm.int.grp")

#####—————Transcription-Translation—————#####
load("/PBA/pg.WGCNA.eigengenes.RData")
WGCNA.info <- read.xlsx("/PBA/Table S13.xlsx")
colnames(WGCNA.info) <- WGCNA.info[1,]
WGCNA.info <- WGCNA.info[-1,]
load("/PBA/PBAtlas.geneExpr.norm.RData")
gexpr.g = gexpr
meta.g = meta

load("/PBA/PBAtlas.proteinExpr.norm.RData")
gexpr.p = combat_edata1
meta.p = meta
#
sampList = intersect(colnames(gexpr.g), colnames(gexpr.p))
x01 = match(sampList, colnames(gexpr.g))
x02 = match(sampList, colnames(gexpr.p))
gexpr2.g = gexpr.g[,x01]
gexpr2.p = gexpr.p[,x02]
meta2.g = meta.g[x01, ]
meta2.p = meta.p[x02, ]
#
genes.pc = intersect(rownames(gexpr2.g), rownames(gexpr2.p))
#
gexpr2.g = gexpr2.g[genes.pc, ]
gexpr2.p = gexpr2.p[genes.pc, ]

###RBP
RBP_ls <- read.csv("/PBA/Homo_sapiens-RBPs.csv")
cRBP_ls <- subset(RBP_ls,RBP.type == "Canonical_RBPs")
cRBP_vec <- sort(unique(cRBP_ls$Gene.symbol))
coRBP <- sort(intersect(cRBP_vec,rownames(gexpr2.p)))

###TF
tfAnnot.guo = read.table("/PBA/TF-Target-information.txt", header = T, sep = "\t")
#
load("/PBA/tftargets.rda")
#
ln = readLines("/PBA/c3.tft.v2024.1.Hs.symbols.gmt")
#
gsea <- list()
for(i in 1:length(ln)){
  rec = unlist(strsplit(as.character(ln[i]), "\t"))
  targets = rec[-(1:2)]
  rem = unlist(strsplit(as.character(rec[1]), "_"))
  tfname = rem[1]
  if(length(gsea[[tfname]]) == 0){
    gsea[[tfname]] = rec
  }else{
    gsea[[tfname]] = unique(c(gsea[[tfname]], rec))
  }
}
#
tf.tot = unique(c(unique(tfAnnot.guo$TF), names(ITFP), names(Marbach2016), names(TRRUST), names(gsea)))
tf2target <- list()
for(i in 1:length(tf.tot)){
  eachtf = tf.tot[i]
  eachtarget = unique(c(ITFP[[eachtf]], Marbach2016[[eachtf]], TRRUST[[eachtf]], gsea[[eachtf]]))
  xi = which(tfAnnot.guo$TF == eachtf)
  eachtarget = unique(c(eachtarget, tfAnnot.guo$target[xi]))
  eachtarget = eachtarget[eachtarget != eachtf]
  tf2target[[eachtf]] = eachtarget
}

for (i in 1:5) {
  cat("\nTF:", names(tf2target)[i], "\n")
  print(head(tf2target[[i]], 5))
}

clean_tf2target <- lapply(names(tf2target), function(tf_name) {
  genes <- tf2target[[tf_name]]
  genes <- genes[!grepl("https://", genes)]
  genes <- genes[!grepl("_TARGET_GENES", genes)]
  genes <- genes[genes %in% rownames(gexpr2.g)]
  genes <- genes[genes != tf_name]
})
names(clean_tf2target) <- names(tf2target)
valid_tfs <- names(clean_tf2target)[names(clean_tf2target) %in% rownames(gexpr2.p)]
tf2target <- clean_tf2target[valid_tfs]

modules <- unique(WGCNA.info$Module)
module_results <- list()

for (module in modules) {
  module_data <- WGCNA.info[WGCNA.info$Module == module, ]

  module_proteins <- module_data$`Gene/Protein Name`[module_data$Type == "Protein"]
  module_genes <- module_data$`Gene/Protein Name`[module_data$Type == "Gene"]
  gene_protein_pairs <- intersect(module_proteins, module_genes)
  if (length(gene_protein_pairs) == 0) next

  module_tfs <- intersect(names(tf2target), module_proteins)
  result_df <- data.frame(TF = character(),
                          Gene = character(),
                          Protein = character(),
                          stringsAsFactors = FALSE)
  for (tf in module_tfs) {
    module_targets <- intersect(tf2target[[tf]], module_genes)
    valid_targets <- intersect(module_targets, gene_protein_pairs)
    
    if (length(valid_targets) > 0) {
      result_df <- rbind(result_df,
                         data.frame(TF = tf,
                                    Gene = valid_targets,
                                    Protein = valid_targets,
                                    stringsAsFactors = FALSE))
    }
  }

  module_rbps <- intersect(coRBP, module_proteins)
  
  module_results[[module]] <- list(
    TF_gene_protein = result_df,
    RB_proteins = module_rbps
  )
  
}

for (module in names(module_results)) {
  cat("\nModule:", module, "\n\n")
  
  if (nrow(module_results[[module]]$TF_gene_protein) > 0 &length(module_results[[module]]$RB_proteins) > 0) {
    cat("TF-RNA-Protein: ",nrow(module_results[[module]]$TF_gene_protein))
    cat("RBP: ",length(module_results[[module]]$RB_proteins),"\n\n")
  }
}

###
correlation_results <- list()
for (module_name in names(module_results)) {
  module_data <- module_results[[module_name]]
  
  if (length(module_data$RB_proteins) == 0 || nrow(module_data$TF_gene_protein) == 0) {
    next
  }
  
  target_proteins <- unique(module_data$TF_gene_protein$Protein)
  rb_proteins <- module_data$RB_proteins
  
  
  target_proteins <- target_proteins[target_proteins %in% rownames(gexpr2.p)]
  rb_proteins <- rb_proteins[rb_proteins %in% rownames(gexpr2.p)]
  
  if (length(target_proteins) == 0 || length(rb_proteins) == 0) {
    next
  }
  
  cor_matrix <- matrix(NA, 
                       nrow = length(target_proteins),
                       ncol = length(rb_proteins),
                       dimnames = list(target_proteins, rb_proteins))
  
  
  for (target in target_proteins) {
    for (rbp in rb_proteins) {
      cor_value <- cor(gexpr2.p[target, ], gexpr2.p[rbp, ], use = "complete.obs") #row-protein, col-RBP
      cor_matrix[target, rbp] <- cor_value
    }
  }
  correlation_results[[module_name]] <- cor_matrix
}

for(m in names(correlation_results)){
  mm <- correlation_results[[m]]
  print(class(mm)) 
  mm <- as.matrix(mm)
  mm[mm == 1] <- 0
  cat(m,"\n\n")
  print(max(mm))
}
