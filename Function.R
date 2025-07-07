ProcessDEanalysisResult_sRdA <- function(region,data_type,pv_cut,fc_cut) {
  if(data_type == "Gene"){
    file_paths <- paste0(dir, region, "/",region, "_", file_names, ".csv")
    
    dfs <- lapply(file_paths, function(path) {
      df <- suppressMessages(read_csv(path))
      df <- subset(df, padj < pv_cut & abs(log2FoldChange) > fc_cut)
    })
    row_counts <- sapply(dfs, nrow)
    df <- data.frame(
      File = file_names,
      RowCount = row_counts,
      Region = region
    )
    
  } else if(data_type == "Protein"){
    file_paths <- paste0(dir, region, "/",region, "_", file_names, ".csv")
    
    dfs <- lapply(file_paths, function(path) {
      df <- suppressMessages(read_csv(path))
      df <- subset(df, LRT.adj_p_value < pv_cut & abs(log2fc) > fc_cut)
    })
    row_counts <- sapply(dfs, nrow)
    df <- data.frame(
      File = file_names,
      RowCount = row_counts,
      Region = region
    )
  } else if(is.null(data_type)) {
    cat("Error: Please enter the `data_type`","\n")
  }
  return(df)
}


ProcessDEanalysisResult_sRdA_TRG <- function(region,data_type,pv_cut,fc_cut,age_group) {
  file_paths <- paste0(dir, region, "/", region, "_", file_names, ".csv")
  if(data_type == "gene"){
    dfs <- lapply(file_paths, function(path) {
      df <- read_csv(path)
      df <- subset(df, padj < pv_cut & abs(log2FoldChange) > fc_cut)
    })
    all.genes <- lapply(dfs,function(x){
      genes <- x$...1
    })
    trg.df <- data.frame(matrix(nrow = 7,ncol = 5))
    colnames(trg.df) <- c("TRGcounts","Allcounts","TRGpct","Region","Age_group")
    trg.df$Region <- region
    trg.df$Age_group <- age_group
    trg.df$Allcounts <- unlist(
      lapply(all.genes,function(x){
        length(x)
      }
      )
    )
    unique_genes_list <- data.frame(
      geneID = as.character(),
      age_group = as.character()
    )
    i <- NULL
    for (i in 1:length(age_group)){
      cur.g <- all.genes[[i]]
      others.g <- unlist(all.genes[-i])
      unique.g <- setdiff(cur.g, others.g)
      unique.df <- data.frame(
        geneID = unique.g,
        age_group = rep(age_group[i],times = length(unique.g))
      )
      unique_genes_list <- rbind(unique_genes_list,unique.df)
      trg.df$TRGcounts[trg.df$Age_group == age_group[i]] <- length(unique.g)
      trg.df$TRGpct <- trg.df$TRGcounts/trg.df$Allcounts
    }
    write.csv(unique_genes_list,file = paste0("/PBA/result_table/TRG/",region,"/trg.csv"))
    
  }else if(data_type == "protein"){
    dfs <- lapply(file_paths, function(path) {
      df <- read_csv(path)
      df <- subset(df, LRT.adj_p_value < pv_cut & abs(log2fc) > fc_cut)
    })
    all.genes <- lapply(dfs,function(x){
      genes <- x$...1
    })
    trg.df <- data.frame(matrix(nrow = 7,ncol = 5))
    colnames(trg.df) <- c("TRGcounts","Allcounts","TRGpct","Region","Age_group")
    trg.df$Region <- region
    trg.df$Age_group <- age_group
    trg.df$Allcounts <- unlist(
      lapply(all.genes,function(x){
        length(x)
      }
      )
    )
    unique_genes_list <- data.frame(
      geneID = as.character(),
      age_group = as.character()
    )
    i <- NULL
    for (i in 1:length(age_group)){
      cur.g <- all.genes[[i]]
      others.g <- unlist(all.genes[-i])
      unique.g <- setdiff(cur.g, others.g)
      unique.df <- data.frame(
        geneID = unique.g,
        age_group = rep(age_group[i],times = length(unique.g))
      )
      unique_genes_list <- rbind(unique_genes_list,unique.df)
      trg.df$TRGcounts[trg.df$Age_group == age_group[i]] <- length(unique.g)
      trg.df$TRGpct <- trg.df$TRGcounts/trg.df$Allcounts
    }
    write.csv(unique_genes_list,file = paste0("/PBA/result_table/TRP/",region,"/trp.csv"))
  }
  return(trg.df)
}

get_ct_related_hcdvg <- function(region, dvg_genes) {
  spec_path <- "/PBA/deconvolution/"
  reg_spec_path <- paste0(spec_path, region, "/", tolower(region), "_hum.spec.xls")
  reg_spec <- read.table(reg_spec_path, header = TRUE, sep = '\t', row.names = 1)
  
  idx.pk <- apply(reg_spec, 2, function(x) rownames(reg_spec)[order(x, decreasing = TRUE)[1:100]])
  rec <- table(unlist(idx.pk))
  topMarkers <- names(rec)[rec <= 2]
  
  hcdvg <- intersect(dvg_genes, topMarkers)
  return(hcdvg)
}

classify_ct_related_dvg <- function(gene, region, high_conf_dvg, spec_df) {
  if (gene %in% high_conf_dvg) {
    return("High-confidence cell-type related")
  }
  
  if (gene %in% spec_df$gene) {
    row <- spec_df[spec_df$gene == gene, -1]
    if (any(row > 0)) {
      return("Low-confidence cell-type related")
    } 
    else {
      return("Non-cell-type related")
    }
  }
  return("Unknown")
}

extract_spectop <- function(file_path) {
  spec = read.table(file_path, header = T, sep = '\t')
  rownames(spec) = as.character(spec[,1])
  idx.pk = apply(spec, 2 ,function(x) {
    y = order(x, decreasing = T)
    y = y[1:200]
  })
  
  rec <- table(rownames(spec)[as.numeric(unlist(idx.pk))])
  topMarkers = names(rec)[rec <= 3]
  sigMat.spec = spec[topMarkers, ]
  return(sigMat.spec)
}

get_ct_mat <- function(df) {
  
  max_column <- apply(df[, -1], 1, function(x) names(x)[which.max(x)])
  result_df <- data.frame(
    ctID = max_column,
    geneID = df$gene
  )
  
  return(result_df)
}

get_ndd_cor <- function(data, drg){
  resEigen <- NULL
  ndd <- NULL
  gexpr <- NULL
  ##
  gexpr = data
  ndd = drg
  ##
  for(i in 1:length(ndd)){
    eachndd = names(ndd)[i]
    genes.ndd = ndd[[i]]
    
    eachgene = intersect(genes.ndd, rownames(gexpr))
    idx = match(eachgene, rownames(gexpr))
    exprMat = gexpr[idx, ]
    
    ####remove low expressed genes and less variable gene
    mean_expr <- rowMeans(exprMat)
    keep_expr <- mean_expr > 0
    
    sd_expr <- apply(exprMat, 1, sd)
    keep_var <- sd_expr > 0
    
    exprMat <- exprMat[keep_expr & keep_var,]
    cat(dim(exprMat), "\n")
    
    ####compute eigengene
    myMat = t(exprMat)
    myEigen = moduleEigengenes(myMat, colors=rep("red",ncol(myMat)))
    geneigen = unlist(myEigen$eigengenes)
    
    #
    if(length(resEigen) == 0){
      resEigen = geneigen
    }else{
      resEigen = rbind(resEigen, geneigen)
    }
  } 
  colnames(resEigen) = colnames(gexpr)
  rownames(resEigen) = names(ndd)
  
  calCor = cor(t(resEigen))
  return(calCor)
}