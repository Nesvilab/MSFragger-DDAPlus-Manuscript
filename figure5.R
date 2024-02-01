#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.

rm(list=ls())

library(dplyr)
library(reshape2)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(ggpubr)


workdir <- "./"
if(!dir.exists(workdir)) dir.create(workdir, recursive = T)



## ---------------------------------------
#  Input file
## ---------------------------------------
manifest_file <- "Z:/yufe/results/msfragger_ddaplus_paper/PXD024427/DDA+/fragpipe-files.fp-manifest"
fp_dda_prot_file <- "Z:/yufe/results/msfragger_ddaplus_paper/PXD024427/DDA/combined_protein.tsv"
fp_ddaplus_prot_file <- "Z:/yufe/results/msfragger_ddaplus_paper/PXD024427/DDA+/combined_protein.tsv"




## ---------------------------------------
#  Prepare quant data
## ---------------------------------------
meta <- read.delim(manifest_file, header = F)
colnames(meta) <- c("sample","condition","replicate","type")
meta$sample <- sapply(meta$sample, 
                      function(x) gsub("20170612_QEP8_JaBA_SA_LT01_V1_LC12_8_[^_]*_|\\.mzML","",basename(x)))
meta$sample[meta$sample=="sample4_2nd_StageTip"] <- "sample04_2nd_StageTip"

# FragPipe dda protein quant 
fp_dda <- read.delim(fp_dda_prot_file, check.names = F)
int_cols <- colnames(fp_dda)[grepl("MaxLFQ Intensity",colnames(fp_dda))]
fp_dda <- fp_dda[,c("Protein ID","Gene",int_cols)]
colnames(fp_dda) <- gsub(" MaxLFQ Intensity","",colnames(fp_dda))
fp_dda[fp_dda==0] <- NA
fp_dda[,3:ncol(fp_dda)] <- apply(fp_dda[,3:ncol(fp_dda)], c(1,2), function(x) log2(x+1))
colnames(fp_dda)[1:2] <- c("ProteinID","Gene")

# FragPipe ddaplus protein quant 
fp_plus <- read.delim(fp_ddaplus_prot_file, check.names = F)
int_cols <- colnames(fp_plus)[grepl("MaxLFQ Intensity",colnames(fp_plus))]
fp_plus <- fp_plus[,c("Protein ID","Gene",int_cols)]
colnames(fp_plus) <- gsub(" MaxLFQ Intensity","",colnames(fp_plus))
fp_plus[fp_plus==0] <- NA
fp_plus[,3:ncol(fp_plus)] <- apply(fp_plus[,3:ncol(fp_plus)], c(1,2), function(x) log2(x+1))
colnames(fp_plus)[1:2] <- c("ProteinID","Gene")




## -----------------------------------------------------
#  Differential analysis for IDHwt-VS-IDHmut
## -----------------------------------------------------
fc <- 2 # fold change threshold
pthres <- 0.05 # pvalue threshold
useAdjustPvalue <- TRUE
min_smp_num <- 3 # minimal sample number control

dlist <- list("DDA"=fp_dda,
              "DDA+"=fp_plus)


filter.missingness <- function(dat, 
                               gc1_list, 
                               gc2_list, 
                               min_smp_num=3){
  nmn <- data.table::rbindlist(lapply(seq_len(nrow(dat)), function(x) {
    df1 <- dat[x, gc1_list]
    df2 <- dat[x, gc2_list]
    
    n1 <- sum(!is.na(df1))
    n2 <- sum(!is.na(df2))
    
    df <- data.frame(notNA1 = n1,
                     notNA2 = n2,
                     stringsAsFactors = F,
                     check.names = F)
    return(df)
  }))
  
  dat_new <- data.frame(dat, nmn, stringsAsFactors = F, check.names = F)
  keep <- which(apply(nmn, 1, function(x) x[1]>min_smp_num && x[2]>min_smp_num))
  cat(paste0("remove ",nrow(nmn)-length(keep)," features with too few sample detections\n"))
  dat_new <- dat_new[keep, colnames(dat)]
  return(dat_new)
}


limma.DE.test <- function(y, 
                          design_matrix, 
                          contrast_matrix, 
                          use.trend=FALSE){
  
  feature_list <- rownames(y)
  
  # fitting data 
  fit <- lmFit(y, design = design_matrix)
  # computing estimated coefficients and standard errors 
  fit.c <- contrasts.fit(fit=fit, contrasts = contrast_matrix)
  # computing moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
  fit.c.eb <- eBayes(fit.c, trend=use.trend)
  
  # log2fc from lmFit
  log2FC <- fit.c.eb$coefficients[, 1]

  # The moderated t-statistic and pvalue
  t.mod <- fit.c.eb$t[, 1]
  t.mod.pvalue <- fit.c.eb$p.value[, 1]
  
  # recover the ordinary t-statistics and the ordinary t-statistic p-values
  t.ord <- fit.c.eb$coefficients[, 1]/fit.c.eb$sigma/fit.c.eb$stdev.unscaled[, 1]
  t.ord.pvalue <- 2*pt(abs(t.ord), df=fit.c.eb$df.residual, lower.tail=FALSE) 
  
  # adjust pvalues
  t.mod.adj.pvalue <- p.adjust(t.mod.pvalue, method = "BH")
  t.ord.adj.pvalue <- p.adjust(t.ord.pvalue, method = "BH")
  
  # create DE result table
  results.eb <- data.frame(log2FC,
                           tvalue=t.mod,
                           tvalue.ord=t.ord,
                           pvalue=t.mod.pvalue,
                           pvalue.ord=t.ord.pvalue,
                           adj.pvalue.ord=t.ord.adj.pvalue,
                           adj.pvalue=t.mod.adj.pvalue)
  
  results.eb$ID <- feature_list
  results.eb$Label <- rep(colnames(contrast_matrix)[1], nrow(results.eb))
  return(results.eb[order(results.eb$adj.pvalue, decreasing = FALSE), ])
}


assign.dep.class <- function(dat, 
                             adj.pval=FALSE, 
                             pval.col=NULL, 
                             pthres=0.05, 
                             fc=2){
  dat$Class<-"None"
  if(is.null(pval.col)){
    if(adj.pval){
      dat[mapply(function(x,y) x & y, dat$log2FC>=log2(fc), dat$adj.pvalue<=pthres), "Class"] <- "Up"
      dat[mapply(function(x,y) x & y, dat$log2FC<=-log2(fc), dat$adj.pvalue<=pthres), "Class"] <- "Down"
    }else{
      dat[mapply(function(x,y) x & y, dat$log2FC>=log2(fc), dat$pvalue<=pthres), "Class"] <- "Up"
      dat[mapply(function(x,y) x & y, dat$log2FC<=-log2(fc), dat$pvalue<=pthres), "Class"] <- "Down"
    }
    dat$ID <- as.character(dat$ID)
    
  }else{
    dat[mapply(function(x,y) x & y, dat$log2FC>=log2(fc), as.numeric(dat[, pval.col])<=pthres), "Class"] <- "Up"
    dat[mapply(function(x,y) x & y, dat$log2FC<=-log2(fc), as.numeric(dat[, pval.col])<=pthres), "Class"] <- "Down"
  }
  
  return(dat)
}

gc.plan <- "IDHwt-VS-IDHmut" # Comparison of IDHwt and IDHmut glioma proteomes
for(i in seq_along(dlist)){
  dtype <- names(dlist)[i]
  data <- dlist[[i]]
  rownames(data) <- data$ProteinID
  
  if(dtype=="DDA"){
    dir_prefix <- "dda"
  }
  if(dtype=="DDA+"){
    dir_prefix <- "dda+"
  }
  
  print(paste0("working on ",dtype))
  outdir <- file.path(workdir, dir_prefix, gc.plan)
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  cat(paste0("working on ", gc.plan," for group comparison analysis\n"))
  gc.list <- unlist(strsplit(gc.plan,"-VS-"))
  
  gc1 <- colnames(data)[grepl("IDHwt",colnames(data))]
  gc2 <- colnames(data)[grepl("IDHmut",colnames(data))]
  cat(paste0(length(gc1)+length(gc2), " cases in total for ", gc.plan, " group comparison\n"))
  cat(paste0(length(gc1), " cases in ", gc.list[1], "\n", length(gc2), " cases in ", gc.list[2], "\n"))
  
  # Make design matrix
  group <- factor(c(rep(gc.list[1], length(gc1)), rep(gc.list[2], length(gc2))))
  design <- model.matrix( ~ 0 + group)
  colnames(design) <- gsub("group","", colnames(design))
  
  # Make contrast matrix
  plan <- paste0(gc.list[1],"-",gc.list[2])
  contrast <- makeContrasts(contrasts = plan, 
                            levels = colnames(design))
  
  cat("*** filter proteins with too many missing values per group\n")
  data_sub <- data[,c(gc1,gc2)]
  data_sub <- filter.missingness(data_sub, gc1, gc2, min_smp_num = min_smp_num)
  
  cat("*** perform DE analysis\n")
  res <- limma.DE.test(y = data_sub,
                       design_matrix = design,
                       contrast_matrix = contrast,
                       use.trend = T)
  
  cat("*** remove NA result\n")
  res <- na.omit(res)
  
  cat("*** assign DE class\n")
  res <- assign.dep.class(res, 
                          adj.pval = useAdjustPvalue, 
                          pthres = pthres, 
                          fc = fc)
  
  cat("*** write DE result to file\n")
  if (!("Gene" %in% colnames(res))) {
    res <- merge(res, data[,1:2], by.x="ID", by.y="ProteinID", all.x=T)
  } 
  write.table(res, 
              file = file.path(outdir,paste0(gc.plan,".All.tsv")), 
              sep = "\t", row.names = F, quote = F)
  write.table(res[res$Class=="Up",], 
              file = file.path(outdir, paste0(gc.plan,".Up.tsv")), 
              sep="\t", row.names = F, quote = F)
  write.table(res[res$Class=="Down",], 
              file = file.path(outdir, paste0(gc.plan,".Down.tsv")), 
              sep="\t", row.names = F, quote = F)
  write.table(res[res$Class!="None",], 
              file = file.path(outdir, paste0(gc.plan,".Up-Down.tsv")), 
              sep="\t", row.names = F, quote = F)
}




## ---------------------------------------------------------------------
#  Functional analysis for IDHwt-VS-IDHmut
## ---------------------------------------------------------------------
dda_de <- read.delim(file.path(workdir, "dda", "IDHwt-VS-IDHmut/IDHwt-VS-IDHmut.All.tsv"))
ddap_de <- read.delim(file.path(workdir, "dda+", "IDHwt-VS-IDHmut/IDHwt-VS-IDHmut.All.tsv"))

dlist <- list("dda" = dda_de,
              "dda+" = ddap_de)


# Perform GO enrichment analysis
do_enrichGO_analysis <- function(d,
                                 gene_col = "Gene",
                                 left_condition = "IDHmut",
                                 right_condition = "IDHwt",
                                 log2fc_thres = 0.1,
                                 pval_thres = 0.05,
                                 ont.type = "BP",
                                 workdir,
                                 prefix = "") {
  
  d1 <- d
  d2 <- d
  all_genes <- d[,gene_col]
  d_left <- d1[(d1$log2FC < -log2fc_thres) & (d1$pvalue < pval_thres), gene_col]
  d_right <- d2[(d2$log2FC > log2fc_thres) & (d2$pvalue < pval_thres), gene_col]
  
  ego1 <- clusterProfiler::enrichGO(
    gene          = d_left,
    universe      = all_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = 'SYMBOL',
    ont           = ont.type,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  ego1_res <- ego1@result
  ego1_res <- ego1_res[order(ego1_res$Count, decreasing = T),]
  ego1_res$label <- rep(left_condition, nrow(ego1_res))
  write.table(ego1_res, file = file.path(workdir, paste0(prefix,"_",ont.type,"_",left_condition,"_enrich_pathway.tsv")),
              row.names = F, sep="\t", quote = F)
  
  
  ego2 <- clusterProfiler::enrichGO(
    gene          = d_right,
    universe      = all_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = 'SYMBOL',
    ont           = ont.type,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  ego2_res <- ego2@result
  ego2_res <- ego2_res[order(ego2_res$Count, decreasing = T),]
  ego2_res$label <- rep(right_condition, nrow(ego2_res))
  write.table(ego2_res, file = file.path(workdir, paste0(prefix,"_",ont.type,"_",right_condition,"_enrich_pathway.tsv")),
              row.names = F, sep="\t", quote = F)
}

for(ont in c("BP","CC","MF")){
  for (i in seq_along(dlist)) {
    dtype <- names(dlist)[i]
    print(paste0("working on ", dtype, " ", ont))
    curr_dir <- file.path(workdir, dtype)
    
    do_enrichGO_analysis(
      d = dlist[[i]],
      gene_col = "Gene",
      left_condition = "IDHmut",
      right_condition = "IDHwt",
      log2fc_thres = 0.1,
      pval_thres = 0.05,
      ont.type = ont,
      workdir = curr_dir,
      prefix = dtype
    )
  }
  
}


# Plot side by side comparison 
plot_enrich_comparison <- function(enrich_list,
                                   cond1 = "",
                                   cond1_label = "DDA",
                                   desc_label = "GO Cellular Component",
                                   cond2 = "",
                                   cond2_label = "DDA+",
                                   topN = 20,
                                   xlim = 220) {
  
  cond1_ont <- enrich_list[[cond1]][1:topN, ]
  cond2_ont <- enrich_list[[cond2]][1:topN, ]
  common_ont <- intersect(cond1_ont[1:topN, ]$ID, cond2_ont[1:topN, ]$ID)
  
  cond1_ont <- cond1_ont[cond1_ont$ID %in% common_ont,]
  cond2_ont <- cond2_ont[cond2_ont$ID %in% common_ont,]
  print(paste0("common: ", length(common_ont), ", cond1: ", nrow(cond1_ont), ", cond2: ", nrow(cond2_ont)))
  
  
  cond1_ont$common_count <- 0
  cond1_ont$unique_count <- sapply(cond1_ont$geneID, function(x) {
    length(unlist(strsplit(x, "\\/")))
  })
  cond2_ont$common_count <- 0
  cond2_ont$unique_count <- sapply(cond2_ont$geneID, function(x) {
    length(unlist(strsplit(x, "\\/")))
  })
  
  for(i in seq_along(common_ont)) {
    curr_ont <- common_ont[i]
    cond1_genes <- unlist(strsplit(cond1_ont[cond1_ont$ID == curr_ont,]$geneID, "\\/"))
    cond2_genes <- unlist(strsplit(cond2_ont[cond2_ont$ID == curr_ont,]$geneID, "\\/"))
    # common count
    common_count <- length(intersect(cond1_genes, cond2_genes))
    cond1_ont[cond1_ont$ID == curr_ont,]$common_count <- common_count
    cond2_ont[cond2_ont$ID == curr_ont,]$common_count <- common_count
    # dda unique count
    cond1_count <- length(setdiff(cond1_genes, cond2_genes))
    # ddap unique count
    cond2_count <- length(setdiff(cond2_genes, cond1_genes))
    
    cond1_ont[cond1_ont$ID == curr_ont,]$unique_count <- cond1_count
    cond2_ont[cond2_ont$ID == curr_ont,]$unique_count <- cond2_count
  }
  cond1_ont <- cond1_ont[, c("ID", "Description", "Count", "common_count", "unique_count")]
  cond2_ont <- cond2_ont[, c("ID", "Description", "Count", "common_count", "unique_count")]
  
  desc_order <- cond2_ont[order(cond2_ont$Count, decreasing = F), ]$Description
  cond1_ont <- melt(cond1_ont, id.vars = c("ID", "Description")) %>%
    mutate(data = cond1_label) %>%
    filter(variable != "Count")
  cond1_ont$variable <- factor(cond1_ont$variable, levels = c("unique_count", "common_count"))
  cond1_ont$Description <- factor(cond1_ont$Description, levels = desc_order)
  
  cond2_ont <- melt(cond2_ont, id.vars = c("ID", "Description")) %>%
    mutate(data = cond2_label) %>%
    filter(variable != "Count")
  cond2_ont$variable <- factor(cond2_ont$variable, levels = c("unique_count", "common_count"))
  cond2_ont$Description <- factor(cond2_ont$Description, levels = desc_order)
  
  p1 <- ggplot(cond1_ont, aes(y = Description, x = value, fill = variable)) +
    geom_bar(position = "stack", stat = "identity") +
    xlim(0, xlim) +
    xlab("Number of genes") +
    ylab("") +
    ggtitle(cond1_label) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_x_reverse(lim = c(xlim, 0), expand = c(0, 0)) +
    scale_fill_manual(
      labels = c("Unique", "Common"),
      values = c(
        "common_count" = "#4285F4",
        "unique_count" = "#EA4335"
      )
    ) +
    guides(fill = guide_legend(title = "Type")) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(color="black",face = "plain", size = 12),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length.y = unit(.25, "cm"),
          legend.text = element_text(face = "plain", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          plot.title = element_text(hjust = 1, color = "black", face = "bold", size = 16),
          plot.margin = unit(c(1,0,1,0), 'lines'),
          panel.spacing = unit(0, "lines"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          text = element_text(family="Arial")) 
  
  p2 <- ggplot(cond2_ont, aes(y=Description, x=value, fill=variable)) +
    geom_bar(position = "stack", stat = "identity") +
    xlab("Number of genes") +
    ylab("") +
    ggtitle(cond2_label) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(lim = c(0, xlim), expand = c(0, 0)) +
    scale_fill_manual(labels = c("Unique", "Common"), values = c("common_count"="#4285F4","unique_count"="#EA4335")) +
    guides(fill=guide_legend(title="Type")) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(color="black",face = "plain", size = 12),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length.y = unit(.25, "cm"),
          legend.text = element_text(face = "plain", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          plot.title = element_text(hjust = 0, color = "black", face = "bold", size = 16),
          plot.margin = unit(c(1,0,1,0), 'lines'),
          panel.spacing = unit(0, "lines"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          text = element_text(family="Arial")) 
  
  desc_df <- data.frame("Description" = desc_order,
                        "text" = desc_order,
                        "x"=1)
  p3 <- ggplot(desc_df, aes(x=x,y=Description)) +
    geom_tile(color="white", fill="white", fontface="bold") +
    geom_text(data=desc_df,aes(y=Description, label=text)) +
    ggtitle(desc_label) +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void() +
    theme(axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, color = "black", face = "plain", size = 16),
          plot.margin = unit(c(1,0,1,0), 'lines'),
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.spacing = unit(0, "lines"),
          panel.border = element_blank(),
          text = element_text(family="Arial")) 
  
  
  plist <- list("cond1"=p1,"desc"=p3,"cond2"=p2)
  names(plist) <- c(cond1_label,"desc",cond2_label)
  return(plist)
  
}

for(enrich in c("BP","CC","MF")){
  for(grp in c("IDHwt", "IDHmut")) {
    enrich_list <- list()
    for (i in seq_along(dlist)) {
      dtype <- names(dlist)[i]
      curr_dir <- file.path(workdir, dtype)
      filename <- file.path(curr_dir, paste0(dtype, "_", enrich, "_",grp,"_enrich_pathway.tsv"))
      print(paste0("working on ", dtype, " ", filename))
      d <- read.delim(filename, check.names = F)
      enrich_list[[paste0(grp, "_", enrich, "_", dtype)]] <- d
    } 
    
    if(enrich == "BP") desc_label <- "GO Biological Processing"
    if(enrich == "CC") desc_label <- "GO Cellular Component"
    if(enrich == "MF") desc_label <- "GO Molecular Function"

    plist <- plot_enrich_comparison(enrich_list = enrich_list,
                                    cond1 = names(enrich_list)[grepl("_dda$",names(enrich_list))],
                                    cond1_label = "DDA",
                                    desc_label = desc_label,
                                    cond2 = names(enrich_list)[grepl("_dda\\+$",names(enrich_list))],
                                    cond2_label = "DDA+",
                                    topN = 20,
                                    xlim = 220)
    
    figure <- ggpubr::ggarrange(plotlist = plist,
                                ncol = 3,
                                widths = c(0.6, 0.54, 0.6),
                                common.legend = TRUE,
                                align = "hv") +
      theme(plot.background = element_rect(fill="white", colour = "white"))
    
    annotate_figure(figure,
                    top = text_grob(paste0(grp, " enriched"), hjust = 2, color = "black", face = "bold", size = 18))
    ggsave(file.path(workdir, paste0(grp,"_",enrich,"_barplot_combined.png")),
           height = 8, width = 10, units = "in", dpi = 300)
    
  }
}
