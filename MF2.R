library(limma)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(RAPToR)
library(ComplexHeatmap)
library(base)
library(circlize)
#####————Species integration————#####
set.seed(100)

###pig
gexpr <- NULL; meta <- NULL
mydata = read.table("/PBA/geneRPKM_PBA.txt", head = T, sep = "\t",check.names = F)
gexpr = t(mydata)
meta <- read.table("/PBA/sampleInfo_PBA.txt", row.names = 1, header=TRUE, sep = '\t')
meta <- meta[colnames(gexpr), ]

idx <- NULL
idx = which(meta$Region == "PFC")

gexpr.pig = gexpr[,idx]
meta.pig = meta[idx, ]
meta2.pig = data.frame(Species = rep("pig", nrow(meta.pig)),
                       SampleId = rownames(meta.pig), 
                       Region = meta.pig$Region, Age = meta.pig$Age, Specimen = meta.pig$Piglet_number)

###NHP
gexpr <- NULL; rpkm <- NULL
load("/PBA/nhp_development_count_rpkm.RData")
pheno = pheno
rpkm = rpkm
rec = unlist(strsplit(rownames(rpkm), split = "|", fixed = T))
dim(rec) = c(2, length(rec)/2)
gsym.nhp = rec[2,]

#remove duplicates
gexpr = rpkm[!duplicated(gsym.nhp), ]
rownames(gexpr) = gsym.nhp[!duplicated(gsym.nhp)]

idx <- NULL
idx = which(pheno$Region %in% c("DFC"))

gexpr.nhp = gexpr[,idx]
meta.nhp = pheno[idx, ]
meta2.nhp = data.frame(Species = tolower(meta.nhp$Species),
                       SampleId = colnames(gexpr.nhp), 
                       Region = meta.nhp$Region, Age = meta.nhp$Age, Specimen = meta.nhp$Brain)

meta2.nhp$Species = gsub("macaque", "rhesus", meta2.nhp$Species)

###Kaessmann Nature 2019
load(file = "/PBA/Kaessmann_Nature2019.allSpecies.rpkm.RData")
idx <- NULL
idx = which(meta.human$Region == "Brain")
gexpr2.human = gexpr2.human[,idx]
meta2.human = meta.human[idx, ]

idx <- NULL
idx = which(meta.rhesus$Region == "Brain")
gexpr2.rhesus = gexpr2.rhesus[,idx]
meta2.rhesus = meta.rhesus[idx, ]

idx <- NULL
idx = which(meta.mouse$Region == "Brain")
gexpr2.mouse = gexpr2.mouse[,idx]
meta2.mouse = meta.mouse[idx, ]

idx <- NULL
idx = which(meta.rat$Region == "Brain")
gexpr2.rat = gexpr2.rat[,idx]
meta2.rat = meta.rat[idx, ]

idx <- NULL
idx = which(meta.rabbit$Region == "Brain")
gexpr2.rabbit = gexpr2.rabbit[,idx]
meta2.rabbit = meta.rabbit[idx, ]

idx <- NULL
idx = which(meta.opossum$Region == "Brain")
gexpr2.opossum = gexpr2.opossum[,idx]
meta2.opossum = meta.opossum[idx, ]

idx <- NULL
idx = which(meta.chicken$Region == "Brain")
gexpr2.chicken = gexpr2.chicken[,idx]
meta2.chicken = meta.chicken[idx, ]

###merge data
gexpr <- NULL; meta <- NULL; rpkm <- NULL
genes = Reduce(intersect, list(rownames(gexpr.pig), rownames(gexpr2.human), rownames(gexpr2.rhesus),
                               rownames(gexpr2.mouse), rownames(gexpr2.rat), rownames(gexpr2.rabbit),
                               rownames(gexpr2.opossum), rownames(gexpr2.chicken), rownames(gexpr.nhp)))

rpkm = cbind(gexpr.pig[genes, ], gexpr2.human[genes, ], gexpr2.rhesus[genes, ],
             gexpr2.mouse[genes, ], gexpr2.rat[genes, ], gexpr2.rabbit[genes, ],
             gexpr2.opossum[genes, ],  gexpr2.chicken[genes, ], gexpr.nhp[genes, ])

meta = rbind(meta2.pig, meta2.human, meta2.rhesus, meta2.mouse, meta2.rat, meta2.rabbit,
             meta2.opossum, meta2.chicken, meta2.nhp)

gexpr = log2(as.matrix(rpkm) + 1)
hiexp = which(apply(gexpr, 1, function(x) {length(which(x > 1))}) > 4)
gexpr = gexpr[hiexp, ]

gexpr = normalizeBetweenArrays(gexpr, method = "quantile")

pca = prcomp(t(gexpr))
pca$x 
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data
summary(pca)
summ <- summary(pca)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

pca.data <- cbind(pca.data, meta) 

p1 <- ggplot(data = pca.data, aes(x = pca.data$X, y = pca.data$Y,fill = pca.data$Species)) +
  geom_point(size = 4, colour= "black", shape = 21, stroke = 1) +
  scale_fill_manual(values = c('red', 'orange','yellow','green','#66FFFF','blue','#9900CC','#FF00FF'),limits=c("human","rhesus","pig","mouse","rat","rabbit","opossum", "chicken")) + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'), 
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_blank()) +
  labs(x = xlab, y = ylab, color="")

ageTrans = read.table("/PBA/PBAtlas.devevo.age.final.txt", header = F, sep = "\t")
idList1 = paste(meta$Species, meta$Age, sep = ".")
idList1 = gsub(" ", "", idList1)
idList2 = paste(ageTrans[,1], ageTrans[,2], sep = ".")
idx <-NULL
idx = match(idList1, idList2)
meta$Day = ageTrans[idx,3]

res = data.frame(pc1 = pca.data$X, Species = meta$Species, Age = meta$Age, Day = meta$Day)
speciesRef = c("human","rhesus","pig", "mouse", "rat","rabbit", "opossum", "chicken")
res$Species = factor(res$Species, rev(speciesRef))

res$DayNorm = log2(res$Day)
for(i in 1:length(speciesRef)){
  eachSpecies = speciesRef[i]
  xi = which(res$Species == eachSpecies)
  eachone = res$DayNorm[xi]/max(res$DayNorm[xi])
  res$DayNorm[xi] = eachone
}    

p2 <- ggplot(data = res, aes(x = pc1, y = Species, size = DayNorm)) +
  geom_point(aes(colour = DayNorm), shape = 15) +
  scale_colour_gradient2(low = "grey88", high = "red")+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.5),
        legend.key = element_rect(fill = 'transparent'), 
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_blank()) 

par(omi = c(0.1, 0.1, 0.1, 0.1))
plts = p1 + p2 + plot_layout(ncol = 1, heights = c(1,0.4))
print(plts)

#####————Age comparison————#####
set.seed(100)

###ctx
gexpr <- NULL; meta <- NULL
mydata = read.table("/PBA/geneRPKM_PBA.txt", head = T, sep = "\t",check.names = F)
gexpr = t(mydata)
meta <- read.table('/PBA/sampleInfo_PBA.txt', row.names = 1, header=TRUE, sep = '\t')
meta <- meta[colnames(gexpr), ]
##
idx <- NULL
idx = which(meta$Region == "PFC")
gexpr.pig = gexpr[,idx]
meta.pig = meta[idx, ]

###hsa development
load("/PBA/BrainSpan.SestanScience2018.RData")
meta = meta
rpkm = rpkm
rec = unlist(strsplit(rownames(rpkm), split = "|", fixed = T))
dim(rec) = c(2, length(rec)/2)
gsym.hsa = rec[2,]

#remove duplicates
gexpr = rpkm[!duplicated(gsym.hsa), ]
rownames(gexpr) = gsym.hsa[!duplicated(gsym.hsa)]

idx <- NULL
idx = which(meta$Regioncode %in% c("DFC", "FC"))
gexpr.hsa = gexpr[,idx]
meta.hsa = meta[idx, ]

genes = intersect(rownames(gexpr.hsa), rownames(gexpr.pig))
gexpr2.hsa = gexpr.hsa[genes, ]
gexpr2.pig = gexpr.pig[genes, ]
rownames(meta.hsa) = colnames(gexpr2.hsa)
rownames(meta.pig) = colnames(gexpr2.pig)

meta.hsa$logDays = log10(meta.hsa$Days)

mat_hsa = ge_im(X = gexpr2.hsa, p = meta.hsa,
                formula = "X ~ s(logDays, bs = 'ts')")

ref_hsa <-make_ref(
  m=mat_hsa,
  n.inter=100,
  t.unit ="hpastegg-laying",
  metadata=list("organism"="hsa",
                "profiling"="whole-organism,bulk",
                "technology"="RNAseq")
)

###predict age
pig_age_predict = ae(gexpr2.pig, ref_hsa)

df01 = data.frame(Species = rep("Human", nrow(meta.hsa)),
                  group = rep("Human_lifespan", nrow(meta.hsa)),
                  Day = meta.hsa$logDays)
df02 = data.frame(Species = rep("Pig", nrow(meta.pig)),
                  group= meta.pig$Age,
                  Day = as.numeric(pig_age_predict$age.estimates[,1]))

df = rbind(df01, df02)

df$group = factor(df$group,  c("Human_lifespan", "E76", "E85","E94", "E104", "E109", "P0","P3","P30"))

xTick=c(50,100,200,500,2000,10000);
p0 <- ggplot(df, aes(x = Day, y = group, colour = Species)) +
  labs(title = "CTX", x="Days (Log10)", y="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(),breaks=log10(xTick),labels=xTick) +
  geom_point(aes(x = Day, y = group), shape = 15, size = 2) +
  scale_colour_manual(name = "Species",values=c("blue", "red"), labels=c("Human", "Pig")) +
  geom_vline(xintercept = log10(266),colour="grey50",alpha=0.5)+
  geom_vline(xintercept = log10(c(69,111,132,167,447,1299,4648,7570)),linetype = "dashed",colour="grey50",alpha=0.5)

#####————Age regulated genes————#####
source("ref_celltypeMarkers.R")
set.seed(100)
genes.out = c("DCX", "SOX4", "SOX11", "CARNS1", "CNDP1", "TESPA1", "SNCG", "SNAP25", "NUSAP1", "TMSB15A", "HMGB3"
              , "PMP2", "HSPA2", "LINC00844", "SIRPA", "DLGAP1-AS4", "LINC00507","BCL11B", "NEUROG2","ENC1"
              ,"XPR1","FOXP2", "PTGER3", "BCL11A", "ZP2", "SLC22A31", "RAB37", "PTK2B","ALDH1A1","PRR35","KRT31")


gexpr2.hsa <- NULL
gexpr2.pig <- NULL
meta.pig <- NULL
meta.hsa <- NULL
pig_age_predict <- NULL
load("/PBA/PBAtlas.devevo.predictAge.ctx.RData")

###remove low quality brains
idx.out = grep("HSB155|HSB194", rownames(meta.hsa))
gexpr2.hsa = gexpr2.hsa[,(-1)*idx.out]
meta.hsa = meta.hsa[(-1)*idx.out, ]

###rename rownames
colnames(gexpr2.hsa) = paste(colnames(gexpr2.hsa), 
                             paste0("W", meta.hsa$Period), sep = ".")
colnames(gexpr2.pig) = paste(colnames(gexpr2.pig), meta.pig$Age, sep = ".")

###merge data
dataMat <- NULL
dataMat = cbind(gexpr2.hsa, gexpr2.pig)
pheno = data.frame(logDays = c(meta.hsa$logDays, pig_age_predict$age.estimates[,1]),
                   Age = c(meta.hsa$Age, meta.pig$Age),
                   Species = c(rep("Human", nrow(meta.hsa)), rep("Pig", nrow(meta.pig))))
idx.sort = order(pheno$logDays, decreasing = F)
dataMat = dataMat[, idx.sort]
pheno = pheno[idx.sort, ]

###log and normalize
dataMat = log2(dataMat + 1)
dataMat = normalizeBetweenArrays(dataMat, method = "quantile")

###limma
design <- model.matrix(~ logDays+Species, data = pheno)
fit <- lmFit(dataMat, design)
contrast_matrix <- makeContrasts(AgeEffect = logDays, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
tops <- topTable(fit2, number=Inf, adjust="fdr")
tops_slc <- subset(tops,abs(tops$logFC)> 2 & tops$adj.P.Val < 0.01)
tops_slc <- subset(tops_slc,!rownames(tops_slc) %in% genes.out)
ctx_tops <- tops_slc[order(tops_slc$logFC,decreasing = T),]
ctx_tops$Region <- "CTX"
ctx_tops$geneID <- rownames(ctx_tops)
colord <- c("geneID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B","Region")
ctx_tops <- ctx_tops[,colord]
##
topGenes <- NULL
topGenes = rownames(tops)[which(abs(tops$logFC)> 2 & tops$adj.P.Val < 0.01)]
topGenes = setdiff(topGenes, genes.out)

rem = dataMat[topGenes, ]
mat1 = NormalizeData(rem)
mat1 = ScaleData(mat1)
rownames(mat1) = rownames(rem)
colnames(mat1) = colnames(rem)

df <- NULL
df = data.frame(Species = pheno$Species, logDays = pheno$logDays)
ha = HeatmapAnnotation(df = df, 
                       col = list(Species = c("Human" =  "blue", "Pig" = "red"),
                                  logDays = colorRamp2(c(0, 3, 6), c("blue", "white", "red"))),
                       na_col = "grey")

x01 = which(colnames(mat1) == "HSB159.DFC.W6")
x02 = which(colnames(mat1) == "HSB121.DFC.W8")
idx.split <- c(rep("Group A", x01), rep("Group B", x02-x01-1), rep("Group C", ncol(mat1)-x02+1))

selected_genes = intersect(rownames(mat1), c(pan.genes, neun.genes, inn.genes, glia.genes))
handed_genes = c("APOE","TOP2A","CLU","AK5","CCK", "TRIM67","BHLHE22","DRAXIN","TUBB2B", "ISLR2", "ST8SIA2")
selected_genes = unique(c(selected_genes, handed_genes))

row_annot = rowAnnotation(pkGenes = anno_mark(at = match(x =selected_genes, table = rownames(mat1)), labels = selected_genes))

ht1 = Heatmap(mat1, name = "CTX age-regulated genes", cluster_columns = F, show_row_dend = F, show_column_names = T
              ,show_row_names = F,right_annotation = row_annot3
              ,top_annotation = ha, column_split = idx.split, column_gap = unit(3, "mm")
              ,row_names_gp = gpar(fontsize=8), column_names_gp = gpar(fontsize=8))
par(omi = c(0.1, 0.1, 0.1, 0.1))
draw(ht1)
