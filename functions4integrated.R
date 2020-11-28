library(ggrepel)
library(Seurat)
library(ggplot2)
library(ggradar)
library(dplyr)
library(pheatmap)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(future)
source("~/Desktop/singlecellanalysis/scRNAseqFunctions.R")

LIGHTGRAY = "#D4D3D2"
set.seed(777)
#time used for outputting plots
timeFmt <- function() {
  format(Sys.time(), "%y%m%d_%H%M%S")
}

COLORSCHEME = c("#color1", "#color2", "#color3") #equal to the number of conditions or factors in sample_col
COLORSCHEME1 = c("#55A4B0", "#5577B0", "#8E55B0") #sample1 (purple colors, as are sample2)
COLORSCHEME2 = c("#3E80A0", "#5E3EA0", "#A03E80", "#3ea05e", "#a05e3e")

#### TOOLS ####

toClipboard <- function(x) {
  clip <- pipe("pbcopy", "w")
  write.table(x, file=clip)
  close(clip)
  print("saved to clipboard")
}

zipTable <- function(x, separator = ":") {
  zipped <- c()
  for (i in 1:length(x)) {
    zipped <- c(zipped, paste0(names(x)[i], separator, " ", (round(x[i],3) * 100))) 
  }
  return(zipped)
}

zip <- function(x, separator = "") {
  zipped <- c()
  for (i in 1:length(x)) {
    zipped <- c(zipped, paste0(names(x)[i], separator, " ", x[i]))
  }
  return(zipped)
}

wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

sortCharArrayByNum <- function(x) {
  positions <- length(x)
  for (i in 1:(positions-1)) {
    for (j in 1:(positions-i)) {
      if (as.numeric(as.character(x[j]))  > as.numeric(as.character(x[j+1]))) { #flip
        temp <- x[j]
        x[j] <- x[j+1]
        x[j+1] <- temp
      } 
    }  
  }
  return(x)
}

getCellProportionFigures <- function(seurat_object, a = "seurat_clusters", b = "condition") {
  
  celltypeNums <- getCellNumbers(seurat_object, celltype_col = a, sample_col = b)
  celltypeNums <- celltypeNums + 0.0001
  celltypeNums.pct <- round(celltypeNums / rowSums(celltypeNums), 3)*100
  # celltypeNums.errbars <- getStdErrOfMean()
  
  celltypeNums.melt <- melt(celltypeNums.pct)
  colnames(celltypeNums.melt) <- c("cluster", "condition", "percent_cells_in_celltype")
  celltypeNums.melt$cluster <- factor(celltypeNums.melt$cluster)
  if (length(unique(celltypeNums.melt$condition)) == 3)  {
    celltypeNums.melt$condition <- factor(celltypeNums.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
  } 
  ggplot(data=celltypeNums.melt, aes(x=cluster, y=percent_cells_in_celltype, group=condition, fill=condition)) + geom_bar(stat="identity", position = "dodge", width=0.7) +
    # geom_errorbar(aes(ymin=ymin, ymax=ymax), position="dodge") +
    labs(y="Percent of Celltype", x="Cluster") + 
    guides(fill=guide_legend(title="Condition")) +
    scale_y_continuous(expand = c(0,0)) + #removes the annoying space between x axis label and the bottom of plot
    scale_fill_manual(values=c(COLORSCHEME2)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1), panel.background = element_rect(fill="white"), 
          panel.grid.major.y = element_line(size=0.5,color=LIGHTGRAY), panel.grid.minor = element_blank(),
          panel.border = element_blank())
  ggsave("cells_per_condition_per_cluster_NOTnorm.pdf", width=12, height=4)
  
  # Normalized
  cond.coefs <- sum(celltypeNums)/colSums(celltypeNums)
  celltypeNums.norm <- t(apply(celltypeNums , 1, function(x) return(cond.coefs * x)))
  celltypeNums.norm.pct <- round(celltypeNums.norm / rowSums(celltypeNums.norm), 3)*100
  celltypeNums.melt <- melt(celltypeNums.norm.pct)
  colnames(celltypeNums.melt) <- c("cluster", "condition", "percent_cells_in_celltype")
  celltypeNums.melt$cluster <- factor(celltypeNums.melt$cluster)
  if (length(unique(celltypeNums.melt$condition)) == 3)  {
    celltypeNums.melt$condition <- factor(celltypeNums.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
  } 
  ggplot(data=celltypeNums.melt, aes(x=cluster, y=percent_cells_in_celltype, group=condition, fill=condition)) + geom_bar(stat="identity", position="dodge", width=0.7) +
    labs(y="Percent of Celltype", x="Cluster") + 
    guides(fill=guide_legend(title="Condition")) +
    scale_fill_manual(values=c(COLORSCHEME2)) +
    scale_y_continuous(expand = c(0,0)) + #removes the annoying space between x axis label and the bottom of plot
    theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1), panel.background = element_rect(fill="white"), 
          panel.grid.major.y = element_line(size=0.5,color=LIGHTGRAY), panel.grid.minor = element_blank(),
          panel.border = element_blank())
  ggsave("cells_per_condition_per_cluster_norm.pdf", width=12, height=4)
}

toClipboard <- function(x) {
  clip <- pipe("pbcopy", "w")
  write.table(x, file=clip)
  close(clip)
  print("saved to clipboard")
}

zipTable <- function(x, separator = ":") {
  zipped <- c()
  for (i in 1:length(x)) {
    zipped <- c(zipped, paste0(names(x)[i], separator, " ", (round(x[i],3) * 100))) 
  }
  return(zipped)
}

zip <- function(x, separator = "") {
  zipped <- c()
  for (i in 1:length(x)) {
    zipped <- c(zipped, paste0(names(x)[i], separator, " ", x[i]))
  }
  return(zipped)
}

wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

sortCharArrayByNum <- function(x) {
  positions <- length(x)
  for (i in 1:(positions-1)) {
    for (j in 1:(positions-i)) {
      # print(x)
      # print(paste0(x[j], " ", x[j+1]))
      # print(as.numeric(as.character(x[j])))
      # print(as.numeric(as.character(x[j+1])))
      # print(as.double(x[j])  > as.double(x[j+1]))
      if (as.numeric(as.character(x[j]))  > as.numeric(as.character(x[j+1]))) { #flip
        temp <- x[j]
        x[j] <- x[j+1]
        x[j+1] <- temp
      } 
    }  
  }
  return(x)
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

#http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

allUMAPs <- function(seurat_object) {
  UMAPPlot(seurat_object, group.by="seurat_clusters", label=TRUE)
  ggsave("umap_byCluster.png", device = "png", units = "in", width = 7, height = 6)
  UMAPPlot(seurat_object, group.by="condition")
  ggsave("umap_byCondition.png", device = "png", units = "in", width = 7, height = 6)
  UMAPPlot(seurat_object, group.by="orig.ident")
  ggsave("umap_bySample.png", device = "png", units = "in", width = 7, height = 6)
}

clusterDE <- function(seurat_object, name, coi="seurat_clusters") {
  Idents(seurat_object) <- seurat_object@meta.data[,coi]
  seurat_object.markers <- FindAllMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, random.seed = 777)
  write.table(seurat_object.markers, paste0("markers_", name, ".txt"), sep = "\t", quote = F, col.names = NA)
  top10 <- seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(seurat_object, features = as.vector(top10$gene)) 
  height80dpi <- (length(as.vector(top10$gene))*12)*(1/80)+1
  ggsave(paste0("heatmap_", name, ".png"), height= height80dpi, width = 8)
}

#DEFINED FUNCTIONS ----
#preprocess
# seurat_object <- fibroid55_12843lane1
pp_plots <- function(seurat_object, dir=".") {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^RPS|RPL")
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^HBB|^HBA2")
  featureinfo <- seurat_object@meta.data[,c("nCount_RNA", "nFeature_RNA", "orig.ident", "percent.mt", "percent.ribo", "percent.hbb")]
  VlnPlot(seurat_object, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hbb"), pt.size = 0)
  ggsave(paste0(dir, "/", seurat_object@project.name, "_violin.pdf"), width=10,height=10)
  FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1)
  ggsave(paste0(dir, "/", seurat_object@project.name, "_featuresCounts.pdf"), width=6,height=4)
}
seurat_object <- all5_med12pos
pp <- function(seurat_object, SCT = FALSE) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^RPS|RPL")
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^HBB|^HBA2")
  featureinfo <- seurat_object@meta.data[,c("nCount_RNA", "nFeature_RNA", "orig.ident", "percent.mt", "percent.ribo", "percent.hbb")]
  #subset on filtered cells and continue
  seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7)
  if (SCT) {
    seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt") #use this to replace scaledata if it's not good enough
  } else {
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1e4)
    seurat_object <- FindVariableFeatures(object = seurat_object, selection.method = 'vst', nfeatures = 2000)
    seurat_object <- ScaleData(object = seurat_object, vars.to.regress = 'percent.mt')
  }
  seurat_object <- RunPCA(object = seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- RunUMAP(object = seurat_object, features = VariableFeatures(seurat_object))
  seurat_object <- FindNeighbors(object = seurat_object, features = VariableFeatures(seurat_object))
  seurat_object <- FindClusters(object = seurat_object, features = VariableFeatures(seurat_object), resolution = 0.4)
  message("saving...")
  saveRDS(seurat_object, paste0(seurat_object@project.name, "_raw_pp.rds"))
  return(seurat_object)
}


umap_withPies <- function(seurat_object, pie_cond = "orig.ident", pie_loc = "seurat_clusters", normalize = TRUE) {
  clusteringRes <- colnames(seurat_object@meta.data)[grep("snn_res", colnames(seurat_object@meta.data))]
  clusteringRes <- unlist(strsplit(clusteringRes, "_"))[length(unlist(strsplit(clusteringRes, "_")))]
  totalCells <- nrow(seurat_object@meta.data)
  samplePercentages <- table(seurat_object@meta.data[,pie_cond]) / totalCells
  
  #make a umap colored by sample group, but stack pie charts on top that show percentage of each sample in each cluster
  # if (class(unique(seurat_object@meta.data[,pie_loc])) %in% c("factor", "numeric")) {
  #   clusters = as.numeric(as.character(unique(seurat_object@meta.data[,pie_loc])))
  # } else {
  #   clusters = as.character(unique(seurat_object@meta.data[,pie_loc]))
  # }
  clusters <- sort(as.character(unique(seurat_object@meta.data[,pie_loc])))
  
  if (length(intersect(clusters, unique(seurat_object@meta.data[,pie_loc]))) != length(clusters)) {
    print("pie_loc not the same")
    print(clusters)
    print(unique(seurat_object@meta.data[,pie_loc]))
    return(NULL)
  }
  
  umap <- seurat_object@reductions$umap@cell.embeddings #SOlist[[x]]@dr$tsne@cell.embeddings
  umap.withmeta <- merge(umap, seurat_object@meta.data, by='row.names')
  rownames(umap.withmeta) <- umap.withmeta[,1]
  umap.withmeta <- umap.withmeta[,-1]
  umap.withmeta$pie.x <- NA
  umap.withmeta$pie.y <- NA
  umap.withmeta$pie.xmed <- NA
  umap.withmeta$pie.ymed <- NA
  for (i in clusters) { #take average of x,y for each cluster for location of pie chart
    cur.x <- umap.withmeta[which(umap.withmeta[,pie_loc] == i),'UMAP_1']
    cur.y <- umap.withmeta[which(umap.withmeta[,pie_loc] == i),'UMAP_2']
    cur.x.mean <- mean(cur.x)
    cur.y.mean <- mean(cur.y)
    cur.x.median <- median(cur.x)
    cur.y.median <- median(cur.y)
    umap.withmeta$pie.x[which(umap.withmeta[,pie_loc] == i)] <- cur.x.mean #mean
    umap.withmeta$pie.y[which(umap.withmeta[,pie_loc] == i)] <- cur.y.mean #mean
    umap.withmeta$pie.xmed[which(umap.withmeta[,pie_loc] == i)] <- cur.x.median #median
    umap.withmeta$pie.ymed[which(umap.withmeta[,pie_loc] == i)] <- cur.y.median #median
  }
  pies <- data.frame("cluster"=clusters, "pie.x"=rep(0,length(clusters)), "pie.y"=rep(0,length(clusters)),"pie.xmed"=rep(0,length(clusters)), "pie.ymed"=rep(0,length(clusters)))
  for (c in unique(seurat_object@meta.data[,pie_cond]) ) {
    pies[,c] <- rep(0, length(clusters))
  }
  
  umap.withmeta[,pie_cond] <- as.character(umap.withmeta[,pie_cond])
  for (j in clusters) { #pie data: number of cells from each clh in each cluster
    for (i in unique(umap.withmeta[,pie_cond])) {
      pies[which(pies$cluster==j),i] <- length(umap.withmeta[,pie_cond][which(umap.withmeta[,pie_cond]==i & umap.withmeta[,pie_loc]==j)])
    }
  }
  
  for (i in 1:length(clusters)) { #pie posiiton: median or mean
    pies$pie.x[i] <- unique(umap.withmeta$pie.x[which(umap.withmeta[,pie_loc]==clusters[i])])
    pies$pie.y[i] <- unique(umap.withmeta$pie.y[which(umap.withmeta[,pie_loc]==clusters[i])])
    pies$pie.xmed[i] <- unique(umap.withmeta$pie.xmed[which(umap.withmeta[,pie_loc]==clusters[i])])
    pies$pie.ymed[i] <- unique(umap.withmeta$pie.ymed[which(umap.withmeta[,pie_loc]==clusters[i])])
  }
  pies$totals <- rowSums(pies[,unique(seurat_object@meta.data[,pie_cond])]) 
  pies.norm <- pies[,unique(seurat_object@meta.data[,pie_cond])]
  cond.coefs <- sum(pies.norm)/colSums(pies.norm)
  pies.norm.f <- t(apply(pies.norm , 1, function(x) return(cond.coefs * x)))
  colnames(pies.norm.f) <- paste0(colnames(pies.norm.f), ".norm")
  pies <- cbind(pies, pies.norm.f)
  pies.norm.percentages <- pies.norm.f / rowSums(pies.norm.f)
  colnames(pies.norm.percentages) <- paste0(colnames(pies.norm.f), ".pct")
  pies <- cbind(pies, pies.norm.percentages)
  if (normalize) {
    #finally plot the data with the pie charts
    ggplot(data=umap.withmeta, aes_string(x="UMAP_1", y="UMAP_2", group=pie_loc, color=pie_loc)) + geom_point(size=0.3) +
      geom_scatterpie(data=pies, aes(x=pie.xmed, y=pie.ymed), cols=paste0(unique(seurat_object@meta.data[,pie_cond]), ".norm")) + theme_classic() + 
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10))) + 
      labs(title = paste0(seurat_object@project.name, " UMAP with ", clusteringRes), 
           subtitle = wrapper(paste0("total cells: ", nrow(seurat_object@meta.data), " and percents: ", paste(zipTable(samplePercentages), collapse="% "), "%" ), width=140) ) 
    #, r=log2(totals)/log2(sum(totals)))
  } else {
    #finally plot the data with the pie charts
    ggplot(data=umap.withmeta, aes_string(x="UMAP_1", y="UMAP_2", group=pie_loc, color=pie_loc)) + geom_point(size=0.3) +
      geom_scatterpie(data=pies, aes(x=pie.xmed, y=pie.ymed), cols=unique(seurat_object@meta.data[,pie_cond])) + theme_classic() + #, r=log2(totals)/log2(sum(totals)))
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10))) + 
      labs(title = paste0(seurat_object@project.name, " UMAP with ", clusteringRes), 
           subtitle = wrapper(paste0("total cells: ", nrow(seurat_object@meta.data), " and percents: ", paste(zipTable(samplePercentages), collapse="% "), "%"), width=140) ) 
  }
  ggsave("umap_with_pie.pdf")
  write.table(pies, "umap_pies_info.txt", sep = "\t", quote = F, col.names = NA)
  return(umap.withmeta)
} 

DEanalysis <- function(seurat_object, comps) {
  
  all_cells <- rownames(seurat_object@meta.data)
  celltypes_left <- seurat_object@meta.data$seurat_clusters #for now -- we'll just do DE between seurat_clusters
  # cell_groupings.backup <- cell_groupings
  for (c in celltypes_left) {
    print(paste0("working on celltype: ", c))
    newdirname <- paste0(c,'_subset')
    if (!(newdirname %in% list.files("."))) {
      dir.create(newdirname)
    }
    #subset data and find variable genes, then do clustering
    newsetname <- seurat_object@project.name
    
    print("making dirs...")
    gridPlots <- FALSE
    dir.create(paste0(newdirname, "/boxplot_andfeatures"))
    dir.create(paste0(newdirname, "/boxplot_andfeatures/boxplots"))
    dir.create(paste0(newdirname, "/boxplot_andfeatures/UMAPs"))
    print("doing DE and making high-rank feature plots")
    res.start <- 0.4 #first iteration of trying to get more clusters
    clusterThreshold <- 2 #this is the minimum number of clusters that will be produced by the below loop
    while (length(unique(new.set@meta.data$seurat_clusters)) < clusterThreshold) { #there isn't enough separation, let's increase the resolution by 0.2 until we have at least 2 clusters
      print(paste0("reclustering with higher res: ", res.start))
      new.set <- FindNeighbors(new.set, force.recalc = TRUE)
      new.set <- FindClusters(new.set, resolution = res.start) #this isn't working exactly as I want, but it'll do for now 
      res.start <- res.start + 0.2
    }
    new.clusters <- unique(new.set@meta.data$seurat_clusters)
    
    for (clust in new.clusters) {
      if (c == "Macrophage") {
        if (length(new.set@meta.data$celltype[which(new.set@meta.data$celltype == clust)]) < 10) {
          next
        } else print(paste0("num cells in cluster: ", length(new.set@meta.data$celltype[which(new.set@meta.data$celltype == clust)])))
      } else {
        if (length(new.set@meta.data$seurat_clusters[which(new.set@meta.data$seurat_clusters == clust)]) < 10) { #if a cluster has fewer than 50 cells, just skip...
          print(paste0("cluster: ", clust, " has fewer than 50 cells (",length(new.set@meta.data$seurat_clusters[which(new.set@meta.data$seurat_clusters == clust)]) ,"), skipping..."))
          next
        } else print(paste0("num cells in cluster: ", length(new.set@meta.data$seurat_clusters[which(new.set@meta.data$seurat_clusters == clust)])))
      }
      print(paste0("making plots for cluster: ", clust))
      dir.create(paste0(newdirname, "/boxplot_andfeatures/boxplots/cluster_", clust))
      dir.create(paste0(newdirname, "/boxplot_andfeatures/UMAPs/cluster_", clust))
      dir.create(paste0(newdirname, "/boxplot_andfeatures/VolcanoPlots"))
      
      d <- rownames(new.set@meta.data[which(new.set@meta.data$seurat_clusters == clust),]) #d <- all current cluster cells
      
      new.set <- SetIdent(new.set, cells = d, value=paste0("cluster_", clust))
      new.set <- SetIdent(new.set, cells = cell_groupings[[c]][!(cell_groupings[[c]] %in% d)], value="other_clusters")
      cur.results <- FindMarkers(new.set, ident.1 = paste0("cluster_", clust), ident.2 = "other_clusters")
      # cur.results$ranking <- ((cur.results$pct.1 + cur.results$pct.2)/2)*cur.results$avg_logFC
      cur.results$ranking <- abs(cur.results$pct.1 - cur.results$pct.2)*abs(cur.results$avg_logFC)
      write.table(cur.results, paste0(newdirname,"/boxplot_andfeatures/cluster_", clust, "_markers.txt"), quote = F, sep = "\t", col.names = NA)
      
      #volcano plots
      cur.results <- cur.results[order(cur.results$p_val_adj),]
      cur.results$labels <- NA
      cur.results$labels[1:20] <- rownames(cur.results)
      cur.results$logPadjf <- -log10(cur.results$p_val_adj + 1e-300)
      cur.results$cumulative_percent_of_cluster <- (cur.results$pct.2 + cur.results$pct.1)/2
      r <- max(abs(cur.results$avg_logFC)) + 1
      volcano <- ggplot2::ggplot(data=cur.results, aes(x = avg_logFC, y = logPadjf, size=cumulative_percent_of_cluster)) + #not yet sure how this should be colored
        geom_point(aes(alpha=0.5)) +
        geom_text_repel(aes(label = labels), point.padding = 0.25, size = 3) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.5, alpha = 0.4) +
        geom_vline(xintercept=c(-1,1), linetype="dashed", size=0.5, alpha = 0.4) +
        xlab(paste("log2 Fold change ( ", clust, " vs. rest of ", c, " subclusters )", sep="")) + scale_x_continuous(limits = c(-r, r)) +
        ylab("-log10 adjusted p-value") +
        theme(legend.key.size = unit(2, 'lines')) +
        theme_bw() + theme(text = element_text(size=10),
                           legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'),
                           legend.position = "right")
      png(paste0(newdirname, "/boxplot_andfeatures/VolcanoPlots/cluster", clust, "_vs_rest", ".png"), width = 720, height=720, units='px', type = "cairo")
      print(volcano)
      dev.off()
      
      #we just want volcanos for now --- 
      next
      
      cur.results <- head(cur.results[order(cur.results$ranking, decreasing = TRUE),], n=25) #9 is easy to plot in a grid (any n^2 number is easy)
      features.to.use <- do.call(cbind.data.frame, list("genes"=rownames(cur.results), "pct1"=cur.results$pct.1, "pct2"=cur.results$pct.2) )
      plots.boxplot <- list()
      plots.features <- list()
      for (g in seq(1,nrow(features.to.use))) {
        print(features.to.use$genes[g])
        
        dir.create(paste0(newdirname,"/boxplot_andfeatures/boxplots/cluster_", clust, "/", features.to.use$genes[g]))
        dir.create(paste0(newdirname,"/boxplot_andfeatures/UMAPs/cluster_", clust, "/", features.to.use$genes[g]))
        
        #create data table with expression values for current gene (maybe scaled & raw?
        cur.expressiondata <- new.set@assays$RNA@data[which(rownames(new.set@assays$RNA@data) == features.to.use$genes[g]),]
        cur.plotdata <- cbind(new.set@meta.data, "exprs"=cur.expressiondata)
        cur.plotdata.umap <- merge(cur.plotdata, new.set@reductions$umap@cell.embeddings, by='row.names')
        
        #regular one
        ggplot(data=cur.plotdata, aes_string(x="seurat_clusters", y="exprs", fill="condition")) + #boxplot
          geom_boxplot() +
          ggtitle(paste0(features.to.use$genes[g], " pct1: ", features.to.use$pct1[g], " pct2: ", features.to.use$pct2[g])) +
          theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
        ggsave(paste0(newdirname,"/boxplot_andfeatures/boxplots/cluster_", clust, "/", features.to.use$genes[g], "/allconds_boxplot.png"), device="png")
        
        #just make the plot, don't nee to do a grid for this...
        ggplot(data=cur.plotdata.umap, aes_string(x="UMAP_1", y="UMAP_2")) +
          geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = "lightgray",high = "purple") +
          ggtitle(paste0(features.to.use$genes[g], " pct1: ", features.to.use$pct1[g], " pct2: ", features.to.use$pct2[g])) +
          theme_bw()
        ggsave(paste0(newdirname,"/boxplot_andfeatures/UMAPs/cluster_", clust, "/", features.to.use$genes[g], "/allconds_UMAP.png"), device="png")
        
        #now do condition-wise (within the cluster, but for only boxplots, tsnes, and umaps
        for (cond in newset.conditions) { #these are not 'allConds'
          cur.cells <- rownames(new.set@meta.data[which(new.set@meta.data$condition == cond),])
          new.set.cond <- SubsetData(new.set, cells = cur.cells)
          cur.expressiondata.cond <- new.set.cond@assays$RNA@data[which(rownames(new.set.cond@assays$RNA@data) == features.to.use$genes[g]),]
          cur.plotdata.cond <- cbind(new.set.cond@meta.data, "exprs"=cur.expressiondata.cond)
          cur.plotdata.umap.cond <- merge(cur.plotdata.cond, new.set.cond@reductions$umap@cell.embeddings, by='row.names')
          
          #boxplot
          ggplot(data=cur.plotdata.cond, aes_string(x="seurat_clusters", y="exprs", fill="condition")) + #boxplot
            geom_boxplot() +
            stat_summary(fun.data = give.n, geom = "text", angle = 90, position = position_dodge2(width=0.9, preserve = "total")) +
            ggtitle(paste0(features.to.use$genes[g])) +
            theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
          ggsave(paste0(newdirname,"/boxplot_andfeatures/boxplots/cluster_", clust, "/", features.to.use$genes[g], "/", cond, "_boxplot.png"), device="png")
          
          #UMAP
          ggplot(data=cur.plotdata.umap.cond, aes_string(x="UMAP_1", y="UMAP_2")) +
            geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = "lightgray",high = "purple") +
            ggtitle(paste0(features.to.use$genes[g])) +
            theme_bw()
          ggsave(paste0(newdirname,"/boxplot_andfeatures/UMAPs/cluster_", clust, "/", features.to.use$genes[g], "/", cond, "_UMAP.png"), device="png")
          
        }
      }
    }
  }
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

#http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

