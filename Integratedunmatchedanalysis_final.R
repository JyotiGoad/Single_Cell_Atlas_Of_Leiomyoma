# code for the seurat analysis of leiomyoma paper

library(data.table)
library(Seurat)
library(DirichletReg)
library(scatterpie)
library(limma)
library(ggrepel)
library(ggplot2)
library(ggradar)
library(dplyr)
library(pheatmap)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(SingleR)
library(matrixStats)
library(monocle3)
require(gam)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
library(colortools)

LIGHTGRAY = "#D4D3D2"
GRAY <- "#e6e6e6"
PURPLE <- "#8533ff"
GREEN <- "#00cc00"
RED_ORANGE <- "#ff531a"
BRIGHT_PURPLE <- "#9445FF"
set.seed(777)

#time used for outputting plots
timeFmt <- function() {
  format(Sys.time(), "%y%m%d_%H%M%S")
}
options(future.globals.maxSize = 8000 * 1024^2)
darkcolors <- c(RColorBrewer::brewer.pal(8,"Set1"), RColorBrewer::brewer.pal(8,"Dark2"), "#FFFFFF")
hemogenes <- c("HBB", "HBD", "HBM", "HBA1", "HBA2", "HBQ1", "HBG2", "HBZ")

# working on linux computer @ ucsf with vpn 2020_08_29
savedir <- ""
homedir <- ""

#LOAD DATA

#med12 negative load ----
fibroid55_11564.data <-  Read10X(data.dir = "Fibroid-55_11564_GEXL_analysis/raw_feature_bc_matrix")
fibroid55_11564 <- CreateSeuratObject(counts = fibroid55_11564.data, project = "fibroid55_11564", min.cells = 3, min.features = 200)
saveRDS(fibroid55_11564, "raw_fibroid55_11564_unfiltered.rds")

fibroid55_12716.data <-  Read10X(data.dir = "Fibroid-55_12716_GEXL_analysis/raw_feature_bc_matrix")
fibroid55_12716 <- CreateSeuratObject(counts = fibroid55_12716.data, project = "fibroid55_12716", min.cells = 3, min.features = 200)
saveRDS(fibroid55_12716, "raw_fibroid55_12716_unfiltered.rds")

fibroid55_12745.data <-  Read10X(data.dir = "Fibroid-55_12745_GEXL_analysis/raw_feature_bc_matrix")
fibroid55_12745 <- CreateSeuratObject(counts = fibroid55_12745.data, project = "fibroid55_12745", min.cells = 3, min.features = 200)
saveRDS(fibroid55_12745, "raw_fibroid55_12745_unfiltered.rds")

#med12 positive load ----
fibroid55_12746_lane1.data <-  Read10X(data.dir = "Fibroid-55_12746_1_GEXL_analysis_1/raw_feature_bc_matrix")
fibroid55_12746_lane1 <- CreateSeuratObject(counts = fibroid55_12746_lane1.data, project = "fibroid55_12746_lane1", min.cells = 3, min.features = 200)
saveRDS(fibroid55_12746_lane1, "raw_fibroid55_12746_lane1_unfiltered.rds")

fibroid55_12746_lane2.data <-  Read10X(data.dir = "Fibroid-55_12746_2_GEXL_analysis/raw_feature_bc_matrix")
fibroid55_12746_lane2 <- CreateSeuratObject(counts = fibroid55_12746_lane2.data, project = "fibroid55_12746_lane2", min.cells = 3, min.features = 200)
saveRDS(fibroid55_12746_lane2, "raw_fibroid55_12746_lane2_unfiltered.rds")

print("loading fibroid55_12382.data")
fibroid55_12382.data <- Read10X(data.dir = paste0("Fibroid-55_12382_GEXL_analysis/raw_feature_bc_matrix/"))
fibroid55_12382 <- CreateSeuratObject(counts = fibroid55_12382.data, project = "fibroid55_12382", min.cells = 3, min.features = 200)

print("loading fibroid55_12640.data")
fibroid55_12640.data <- Read10X(data.dir = paste0("Fibroid-55_12640_GEXL_analysis/raw_feature_bc_matrix/"))
fibroid55_12640 <- CreateSeuratObject(counts = fibroid55_12640.data, project = "fibroid55_12640", min.cells = 3, min.features = 200)

print("loading fibroid55_12843lane1.data")
fibroid55_12843lane1.data <- Read10X(data.dir = paste0("Fibroid-55_12843_1_GEXL_analysis/raw_feature_bc_matrix/"))
fibroid55_12843lane1 <- CreateSeuratObject(counts = fibroid55_12843lane1.data, project = "fibroid55_12843lane1", min.cells = 3, min.features = 200)

print("loading fibroid55_12843lane2.data")
fibroid55_12843lane2.data <- Read10X(data.dir = paste0("Fibroid-55_12843_2_GEXL_analysis/raw_feature_bc_matrix/"))
fibroid55_12843lane2 <- CreateSeuratObject(counts = fibroid55_12843lane2.data, project = "fibroid55_12843lane2", min.cells = 3, min.features = 200)

#merge the two lanes
print("merging the two 12843 lanes...")
fibroid55_12843 <- merge(fibroid55_12843lane1, fibroid55_12843lane2, add.cell.ids = c("L001", "L002"))
saveRDS(fibroid55_12843, "saves/raw_fibroid55_12843_unfiltered.rds")

print("loading fibroid55_12906lane1.data")
fibroid55_12906lane1.data <- Read10X(data.dir = paste0("Fibroid55-12906-1_GEXL_analysis/raw_feature_bc_matrix/"))
fibroid55_12906lane1 <- CreateSeuratObject(counts = fibroid55_12906lane1.data, project = "fibroid55_12906lane1", min.cells = 3, min.features = 200)

print("loading fibroid55_12906lane2.data")
fibroid55_12906lane2.data <- Read10X(data.dir = paste0("Fibroid-55_12843_2_GEXL_analysis/raw_feature_bc_matrix/"))
fibroid55_12906lane2 <- CreateSeuratObject(counts = fibroid55_12906lane2.data, project = "fibroid55_12906lane2", min.cells = 3, min.features = 200)
saveRDS(fibroid55_12906lane2, "saves/raw_fibroid55_12906lane2_unfiltered.rds")

#myo
print("loading myometrium061919.data")
myometrium061919.data <- Read10X(data.dir = paste0("061919Myometrium_GEXL_fastqs_analysis/raw_feature_bc_matrix/"))
myometrium061919 <- CreateSeuratObject(counts = myometrium061919.data, project = "myometrium061919", min.cells = 3, min.features = 200)
myometrium061919 <- readRDS("saves/raw_myometrium061919_unfiltered.rds")

print("loading myometrium11911.data")
myometrium11911.data <- Read10X(data.dir = paste0("Myometrium_PBSP_GEXL_analysis/raw_feature_bc_matrix/"))
myometrium11911 <- CreateSeuratObject(counts = myometrium11911.data, project = "myomerium11911", min.cells = 3, min.features = 200)
myometrium11911 <- readRDS("saves/raw_myometrium11911_unfiltered.rds")

print("loading myometrium55_12640.data")
myometrium55_12640.data <- Read10X(data.dir = paste0("Myometrim-55_12640_GEXL_analysis/raw_feature_bc_matrix/"))
myometrium55_12640 <- CreateSeuratObject(counts = myometrium55_12640.data, project = "myometrium55_12640", min.cells = 3, min.features = 200)
myometrium55_12640 <- readRDS("saves/raw_myometrium55_12640_unfiltered.rds")

print("loading myometrium55_11564.data")
myometrium55_11564.data <- Read10X(data.dir = paste0("Myometrium-55_11564_GEXL_analysis/raw_feature_bc_matrix/"))
myometrium55_11564 <- CreateSeuratObject(counts = myometrium55_11564.data, project = "myometrium55_11564", min.cells = 3, min.features = 200) 
myometrium55_11564 <- readRDS("saves/raw_myometrium55_11564_unfiltered.rds")

print("loading myometrium55_12745.data")
myometrium55_12745.data <- Read10X(data.dir = paste0("Myometrium_55_12745_GEXL_analysis/raw_feature_bc_matrix/"))
myometrium55_12745 <- CreateSeuratObject(counts = myometrium55_12745.data, project = "myometrium55_12745", min.cells = 3, min.features = 200)
myometrium55_12745 <- readRDS("saves/myometrium55_12745_unfiltered.rds")


#combine lanes
fibroid55_12843 <- merge(x=fibroid55_12843lane1, y=fibroid55_12843lane2, add.cell.ids = c("lane1", "lane2"))
fibroid55_12843@project.name <- "fibroid55_12843"
fibroid55_12906 <- merge(x=fibroid55_12906lane2, y=fibroid55_12906lane1, add.cell.ids = c("lane2", "lane1")) #names are backwards so plots are nicer
fibroid55_12906@project.name <- "fibroid55_12906"
all5_med12pos_names <- c("fibroid55_12382", "fibroid55_12906", "fibroid55_12843", "fibroid55_12640", "fibroid061919")
all5_med12pos_samplenames <- c("fibroid55_12382", "fibroid55_12906lane1","fibroid55_12906lane2", "fibroid55_12843lane1", "fibroid55_12843lane1", "fibroid55_12640", "fibroid061919")

#combine lanes for med12 negative sample
fibroid55_12746 <- merge(x=fibroid55_12746_lane1, y=fibroid55_12746_lane2, add.cell.ids = c("lane1", "lane2"))
fibroid55_12746@project.name <- "fibroid55_12746"
all_med12neg_names <- c("fibroid55_11564", "fibroid55_12716", "fibroid55_12745", "fibroid55_12746", "fibroid55_12746_lane1", "fibroid55_12746_lane2")
all_med12neg_samplenames <- c("fibroid55_11564", "fibroid55_12716", "fibroid55_12745", "fibroid55_12746_lane1", "fibroid55_12746_lane2")
# do a single sample QC and processing here. ----
for (x in c(fibroid55_12746, fibroid55_11564, fibroid55_12716, fibroid55_12745)) {
  pp_plots(x, dir = "med12neg")
}

all5_med12pos <- merge(x=fibroid55_12382, y=list(fibroid55_12906, fibroid55_12843, fibroid55_12640, fibroid061919), add.cell.ids = all5_med12pos_names)
pp_plots(all5_med12pos)


##### Merging ####
# myo list #####
myometrium55_11911 <- readRDS("saves/raw_myometrium11911_unfiltered.rds")
myometrium55_12341 <- readRDS("saves/raw_myometrium061919_unfiltered.rds")
myometrium55_11564 <- readRDS("saves/raw_myometrium55_11564_unfiltered.rds")
myometrium55_12640 <- readRDS("saves/raw_myometrium55_12640_unfiltered.rds")
myometrium55_12745 <- readRDS("saves/raw_myometrium55_12745_unfiltered.rds")
myo_all_final.samples <- list(myometrium55_11911, myometrium55_12341, myometrium55_11564, myometrium55_12640, myometrium55_12745)
names(myo_all_final.samples) <- c("myometrium55_11911", "myometrium55_12341", "myometrium55_11564", "myometrium55_12640", "myometrium55_12745")
myo_all_final <- merge(myo_all_final.samples[[1]], y=myo_all_final.samples[2:5], add.cell.ids = names(myo_all_final.samples)) 
# saveRDS(myo_all_final, "../saves/FINAL_myo_all_raw_unfiltered.rds")

# med12pos list ####
med12pos55_12341 <- readRDS("saves/raw_fibroid061919_unfiltered.rds")
med12pos55_12882 <- readRDS("saves/raw_fibroid55_12382_unfiltered.rds")
med12pos55_12460 <- readRDS("saves/raw_fibroid55_12640_unfiltered.rds")
med12pos55_12843 <- readRDS("saves/raw_fibroid55_12843_unfiltered.rds") #both lanes combined -> both look great
med12pos55_12843@meta.data$orig.ident <- "med12pos55_12843"
med12pos55_12843 <- SetIdent(med12pos55_12843, value = "orig.ident")
med12pos55_12906 <- readRDS("saves/raw_fibroid55_12906lane2_unfiltered.rds") #taking only the 2nd lane
med12pos_all_final.samples <- list(med12pos55_12341, med12pos55_12882, med12pos55_12460, med12pos55_12843, med12pos55_12906)
names(med12pos_all_final.samples) <- c("med12pos55_12341", "med12pos55_12882", "med12pos55_12460", "med12pos55_12843", "med12pos55_12906")
med12pos_all_final <- merge(med12pos_all_final.samples[[1]], y=med12pos_all_final.samples[2:5], add.cell.ids = names(med12pos_all_final.samples)) 
# saveRDS(med12pos_all_final, "../saves/FINAL_med12pos_all_raw_unfiltered.rds")

#Fibroid 12745 not used
med12pos_all_final <- readRDS("saves/FINAL_med12pos_all_raw_unfiltered_fibroid12745removed.rds")

# med12neg list ####
med12neg55_11564 <- readRDS("saves/raw_fibroid55_11564_unfiltered.rds")
med12neg55_12716 <- readRDS("saves/raw_fibroid55_12716_unfiltered.rds")
med12neg55_12746 <- readRDS("saves/raw_fibroid55_12746_unfiltered.rds") # both lanes look ok 
med12neg_all_final.samples <- list(med12neg55_11564, med12neg55_12716,  med12neg55_12746) #med12neg55_12745,
names(med12neg_all_final.samples) <- c("med12neg55_11564", "med12neg55_12716", "med12neg55_12746") # "med12neg55_12745" sample has a TON of cells of one type, not sure what
med12neg_all_final <- merge(med12neg_all_final.samples[[1]], y=med12neg_all_final.samples[2:3], add.cell.ids = names(med12neg_all_final.samples)) 
med12neg_all_final@meta.data$orig.ident[grep("fibroid55_12746", med12neg_all_final@meta.data$orig.ident)] <- "fibroid55_12746"  
Idents(med12neg_all_final) <- med12neg_all_final@meta.data$orig.ident
saveRDS(med12neg_all_final, "saves/FINAL_med12neg_all_raw_unfiltered_new.rds")

##### Integration #### 
hemogenes <- c("HBB", "HBD", "HBM", "HBA1", "HBA2", "HBQ1", "HBG2", "HBZ")
all_int_final.list <- list(myo_all_final, med12pos_all_final, med12neg_all_final)
myo_all_final@project.name <- "myo_all_final"
med12pos_all_final@project.name <- "med12pos_all_final"
med12neg_all_final@project.name <- "med12neg_all_final"
allStats <- c()
all_int_final.list <- lapply(X = all_int_final.list, FUN = function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^RP[SL]") #pattern = "^RPS|RPL"
  curhemogenes <- intersect(hemogenes, rownames(seurat_object@assays$RNA@counts)) #because all of the hemogenes aren't always in the feature set
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, features = curhemogenes) 
  
  #subset on filtered cells and continue
  statTable <- c()
  statTable <- c(statTable, nrow(seurat_object@meta.data))
  seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  statTable <- c(statTable, nrow(seurat_object@meta.data))
  seurat_object <- subset(x = seurat_object, subset = percent.mt < 7)
  statTable <- c(statTable, nrow(seurat_object@meta.data))
  seurat_object <- subset(x = seurat_object, subset = percent.hbb < 1) 
  statTable <- c(statTable, nrow(seurat_object@meta.data))
  allStats <<- cbind(allStats, statTable)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  message(paste0("finished with: ", seurat_object@project.name))
  return(seurat_object)
})
rownames(allStats) <- c("no_filtering", "filtered_by_nfeature", "filtered_by_mito", "filtered_by_hbb")
colnames(allStats) <- c("myo_all_final", "med12pos_all_final", "med12neg_all_final")
allStats <- as.data.frame(allStats)
allStats$totals <- as.vector(rowSums(as.matrix(allStats)))
write.table(allStats, "new_filtration_stats_20200221.txt", quote = F, sep = "\t", col.names = NA)

all_int_final.anchors <- FindIntegrationAnchors(object.list = all_int_final.list, dims = 1:20)
all_int_final <- IntegrateData(anchorset = all_int_final.anchors, dims = 1:20)
all_int_final <- FindVariableFeatures(all_int_final, nfeatures=2000)
all_int_final <- ScaleData(all_int_final, vars.to.regress = "percent.mt")
all_int_final <- RunPCA(object = all_int_final, features = VariableFeatures(object = all_int_final))
all_int_final <- RunUMAP(object = all_int_final, dims=1:20)
all_int_final <- FindNeighbors(object = all_int_final, features = VariableFeatures(all_int_final))
all_int_final <- FindClusters(object = all_int_final, features = VariableFeatures(all_int_final), resolution = 0.4)
saveRDS(all_int_final, "all_int_final_pp_clustered_removal_of_hbb1percentOrMoreCells.rds")
#add conditions
all_int_final@meta.data$condition <- "med12pos"
all_int_final@meta.data$condition[grep("myometrium", all_int_final@meta.data$orig.ident)] <- "myometrium"
all_int_final@meta.data$condition[which(all_int_final@meta.data$orig.ident %in% all_med12neg_names)] <- "med12neg" #use sample list to assign condition

all_int_final.noHBB <- ScaleData(all_int_final, vars.to.regress = "percent.hbb")
all_int_final.noHBB <- FindVariableFeatures(all_int_final.noHBB)
all_int_final.noHBB <- RunPCA(object = all_int_final.noHBB, features = VariableFeatures(object = all_int_final.noHBB))
all_int_final.noHBB <- RunUMAP(object = all_int_final.noHBB, dims=1:20)
all_int_final.noHBB <- FindNeighbors(object = all_int_final.noHBB, features = VariableFeatures(all_int_final.noHBB))
all_int_final.noHBB <- FindClusters(object = all_int_final.noHBB, features = VariableFeatures(all_int_final.noHBB), resolution = 0.4)

all_int_final.noHBB@meta.data$condition <- "med12pos"
all_int_final.noHBB@meta.data$condition[grep("myometrium", all_int_final.noHBB@meta.data$orig.ident)] <- "myometrium"
all_int_final.noHBB@meta.data$condition[which(all_int_final.noHBB@meta.data$orig.ident %in% all_med12neg_names)] <- "med12neg" #use sample list to assign condition

#remove all cells with actual hbb in them
hemogenes <- c("HBB", "HBD", "HBM", "HBA1", "HBA2", "HBQ1", "HBG2", "HBZ")
all_int_final.noHBBcells <- all_int_final.noHBB
quantile(all_int_final.noHBB@meta.data$percent.hbb, 0:100/100)[c(50:101)] #to check
all_int_final.noHBBcells <- subset(all_int_final.noHBB, subset=percent.hbb < 1)
quantile(all_int_final.noHBBcells@meta.data$percent.hbb, 0:100/100)[c(50:101)] #to check
nrow(all_int_final.noHBBcells@meta.data) #33,361
all_int_final.noHBBcells <- ScaleData(all_int_final.noHBBcells) #wasn't run yet
all_int_final.noHBBcells <- FindVariableFeatures(all_int_final.noHBBcells)
all_int_final.noHBBcells <- RunPCA(object = all_int_final.noHBBcells, features = VariableFeatures(object = all_int_final.noHBBcells))
all_int_final.noHBBcells <- RunUMAP(object = all_int_final.noHBBcells, dims=1:20)
all_int_final.noHBBcells <- FindNeighbors(object = all_int_final.noHBBcells, features = VariableFeatures(all_int_final.noHBBcells))
all_int_final.noHBBcells <- FindClusters(object = all_int_final.noHBBcells, features = VariableFeatures(all_int_final.noHBBcells), resolution = 0.4)
# saveRDS(all_int_final.noHBBcells, "saves/all_int_final.noHBBcells.rds")
UMAPPlot(all_int_final.noHBBcells)
ggsave("all_int_final_nohbb_cells.pdf")
VlnPlot(all_int_final.noHBBcells, features=hemogenes, assay = "RNA")
all_int_final.noHBBcells.markers <- FindAllMarkers(all_int_final.noHBBcells, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBBcells.markers, "all_int_final.noHBBcells.markers.txt", sep = "\t", quote = F, col.names = NA)
all_int_final.noHBBcells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- top10$gene
DoHeatmap(all_int_final.noHBBcells, features = top10genes) #do this

height80dpi <- (length(top10genes)*12)*(1/80)+1 ##formula for size is (num_genes*12ptFont)*(1inch/80dpi) = height   ---- 
ggsave("nohbbcells/all_int_final_noHBBcells_heatmap_res0.4.pdf", height=height80dpi, width=8)


#### Immune hbb removal #### 
hemogenes <- c("HBB", "HBD", "HBM", "HBA1", "HBA2", "HBQ1", "HBG2", "HBZ")
all_int_final.noHBBcells <- all_int_final.noHBB
quantile(all_int_final.noHBB@meta.data$percent.hbb, 0:100/100)[c(50:101)] #to check
all_int_final.noHBBcells <- subset(all_int_final.noHBB, subset=percent.hbb < 1)
quantile(all_int_final.noHBBcells@meta.data$percent.hbb, 0:100/100)[c(50:101)] #to check
nrow(all_int_final.noHBBcells@meta.data) #33,361
all_int_final.noHBBcells <- ScaleData(all_int_final.noHBBcells) #wasn't run yet
all_int_final.noHBBcells <- FindVariableFeatures(all_int_final.noHBBcells)
all_int_final.noHBBcells <- RunPCA(object = all_int_final.noHBBcells, features = VariableFeatures(object = all_int_final.noHBBcells))
all_int_final.noHBBcells <- RunUMAP(object = all_int_final.noHBBcells, dims=1:20)
all_int_final.noHBBcells <- FindNeighbors(object = all_int_final.noHBBcells, features = VariableFeatures(all_int_final.noHBBcells))
all_int_final.noHBBcells <- FindClusters(object = all_int_final.noHBBcells, features = VariableFeatures(all_int_final.noHBBcells), resolution = 0.4)
# saveRDS(all_int_final.noHBBcells, "saves/all_int_final.noHBBcells.rds")
UMAPPlot(all_int_final.noHBBcells)
ggsave("all_int_final_nohbb_cells.pdf")
VlnPlot(all_int_final.noHBBcells, features=hemogenes, assay = "RNA")
all_int_final.noHBBcells.markers <- FindAllMarkers(all_int_final.noHBBcells, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBBcells.markers, "all_int_final.noHBBcells.markers.txt", sep = "\t", quote = F, col.names = NA)
all_int_final.noHBBcells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- top10$gene
DoHeatmap(all_int_final.noHBBcells, features = top10genes) #do this


#### Final Marker List #### 
# stromal
SMC_markers <- c("Myh11", "Tagln", "ACTA2", "CNN1", "DES", "CALD1")
fibroblast_markers <- c("VIM", "ALDH1", "CD90", "FN1", "DCN", "OGN", "MGP", "COL1A1", "COL1A2", "COL3A1" )
Myofibroblasts_markers <- c("RNF213", "POSTN", "TNC", "TREM2", "IL1B", "SPP1", "LGAL3", "CCR2", "TNFSF12")
Lymphatic_endothelial_markers <- c("STMN1","CD9","FABP5","ADM", "C10ORF10")
Endothelial_markers <- c("PECAM1", "CD31", "CDH11", "FABP4", "KDR", "VEGFR2", "FLK",  "CDH5", "VWF" )
hbb_markers <- c("HBB",   "HBA1"  ,"HBA2" , "HBB"   ,"HBA2" , "HBB" ,  "HBA1",  "HBB"  , "HBA2", "HBEGF")
stromal_markers <- list(SMC_markers, fibroblast_markers, Myofibroblasts_markers, Lymphatic_endothelial_markers, 
                        Endothelial_markers)
names(stromal_markers) <- c("SMC_markers", "fibroblast_markers", "Myofibroblasts_markers",  
                            "Lymphatic_endothelial_markers", "Endothelial_markers")

#immune 
Platelet_markers <- c("PPBP")
Plasma_markers <- c("IGLC2", "IGHG1", "IGKC", "IGHG2", "IGHG3", "IGHGP", "IGLC3", "JCHAIN", "IGHA1", "IGHG4", "IGHA2", "IGHM", "IGLV3-1", "IGLC7", "MZB1", "CD79A", "SSR4", "IL16")
NaiveCD4pT_markers <- c("IL7R", "CCR7")
CD4pT_markers <- c("CD3D", "IL7R", "CCR7")
CD8pT_markers <- c("CD3D", "CD8A")
NK_markers <- c("GNLY", "NKG7", "NCAM", "GZMA", "GZMB", "PRF1")
Macrophage_markers <- c("CSF1R")
Monocyte_markers <- c("CD14", "LYZ")
FCGR3ApMono_markers <- c("FCGR3A", "MS4A7")
B_markers <- c("MS4A1", "LTB", "CD37", "CD79B", "CD52", "HLA-DQB1", "TNFRSF13C", "TCL1A", "LINC00926", "STAG3", "IGHD", "BANK1", "IRF8", "BIRC3", "P2RX5", "RP11-693J15.5", "RP5-887A10.1", "VPREB3", "CD22", "CD74", "SELL")
Conventional_dendritic_markers <- c("FCER1A", "CST3")
Mast_markers <- c("TPSAB1")
Plasmacytoid_dendritic_markers <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")
Microfibril_markers <- c("MFAP5")
Myofibroblast_resembling_FB_markers <- c("DES","MFAP4", "OGN", "S100A4")
Vascular_endothelial_markers <- c("CKR1", "PCDH17", "FAM167B", "AQP1", "HP", "IFI27")
Naive_T_markers <- c("CCR7", "SELL") 
Activated_B_and_T_markers <- c("CREM", "CD69")
matureNeutrophilMarkers <- "ALPL, IL8RB, FCGR3B, SEMA3C, HM74, SOD2, FCGR3A, IL-8, STHM, IL8RA, FCGR2A, CSF3R, NCF2, AOAH"
matureNeutrophilMarkers <- strsplit(matureNeutrophilMarkers, ", ")[[1]]
matureNeutrophilMarkers[8] <- "IL8"
immune_markers <- list(NaiveCD4pT_markers, CD4pT_markers, CD8pT_markers, NK_markers, Macrophage_markers,
                       Monocyte_markers, FCGR3ApMono_markers, B_markers, Conventional_dendritic_markers,
                       Mast_markers, Plasmacytoid_dendritic_markers, Microfibril_markers, Myofibroblast_resembling_FB_markers,
                       Vascular_endothelial_markers, Naive_T_markers, Activated_B_and_T_markers, Plasma_markers, 
                       Platelet_markers, matureNeutrophilMarkers)
names(immune_markers) <- c("NaiveCD4pT_markers", "CD4pT_markers", "CD8pT_markers", "NK_markers", "Macrophage_markers",
                           "Monocyte_markers", "FCGR3ApMono_markers", "B_markers", "Conventional_dendritic_markers",
                           "Mast_markers", "Plasmacytoid_dendritic_markers", "Microfibril_markers", "Myofibroblast_resembling_FB_markers",
                           "Vascular_endothelial_markers", "Naive_T_markers", "Activated_B_and_T_markers", "Plasma_markers", 
                           "Platelet_markers", "matureNeutrophilMarkers")

#contaminating
Stressed_dying_markers <- c("HSPB1", "DNAJB6", "HSPH1", "GADD45B")
contaminating <- c(Stressed_dying_markers)

#### define cell clusters ####

#make feature plots - stromal
seurat_object <- all_int_final.noHBB
for (fl in 1:length(stromal_markers)) {
  message(paste0("working: ",fl))
  for (f in stromal_markers[[fl]]) {
    curdir <- paste0("hbb_regressed/features/stromal/", names(stromal_markers)[fl] )
    dir.create(curdir)
    if (!f%in%rownames(seurat_object@assays$RNA@data)) {
      message("not in integrated set in @assays$RNA@data")
      next 
    }
    p1 <- FeaturePlot(seurat_object, features=f, reduction = "umap", pt.size=0.5) 
    p2 <- UMAPPlot(seurat_object)
    plot_grid(p1, p2)
    ggsave(paste0("hbb_regressed/features/stromal/", names(stromal_markers)[fl], "/", f, ".pdf"), width=6, height=3)
  }
}

#make feature plots - immune
for (fl in 1:length(immune_markers)) {
  message(paste0("working: ",fl))
  curdir <- paste0("hbb_regressed/features/immune/", names(immune_markers)[fl] )
  dir.create(curdir)
  for (f in immune_markers[[fl]]) {
    if (!f%in%rownames(seurat_object@assays$RNA@data)) {
      message("not in integrated set in @assays$RNA@data")
      next 
    }
    p1 <- FeaturePlot(seurat_object, features=f, reduction = "umap", pt.size=0.5) 
    p2 <- UMAPPlot(seurat_object)
    plot_grid(p1, p2)
    ggsave(paste0("hbb_regressed/features/immune/", names(immune_markers)[fl], "/", f, ".pdf"), width=6, height=3)
  }
}

#make feature plots - contaminating
for (f in contaminating) {
  if (!f%in%rownames(seurat_object@assays$RNA@data)) {
    message("not in integrated set in @assays$RNA@data")
    next 
  }
  p1 <- FeaturePlot(seurat_object, features=f, reduction = "umap", pt.size=0.5) 
  p2 <- UMAPPlot(seurat_object)
  plot_grid(p1, p2)
  ggsave(paste0("integrated_all/features/contaminating/", f, ".pdf"), width=6, height=3)
}


### violin features ####
seurat_object <- all_int_final.noHBB
cell_classes <- list(stromal_markers, immune_markers, contaminating)
names(cell_classes) <- c("stromal_markers", "immune_markers", "contaminating")
for (cc in 1:length(cell_classes)) {
  print(paste0("working on class: ", cell_classes[[cc]]))
  curdir <- paste0("hbb_regressed/features/",names(cell_classes)[cc] )
  dir.create(curdir)
  for (fl in 1:length(cell_classes[[cc]])) {
    #it'll just catch the error itself if the feature isn't there
    if (names(cell_classes)[cc] == "contaminating") {
      VlnPlot(seurat_object, features=cell_classes[[cc]], pt.size = 0)
      ggsave(paste0("hbb_regressed/features/contaminating/all_contaminants.pdf"), width=11, height=7)
    } else {
      VlnPlot(seurat_object, features=cell_classes[[cc]][[fl]], pt.size = 0)
      ggsave(paste0("hbb_regressed/features/",names(cell_classes)[cc], "/", names(cell_classes[[cc]])[fl], ".pdf"), width=11, height=7)
    }
  }
}

### Assign Celltypes Major ####

#reasign full data cell types
Endos <- rownames(all_int_final.noHBB.mes.reclust@meta.data[which(all_int_final.noHBB.mes.reclust$celltype == "Endothelial"),])
SMCs <- rownames(all_int_final.noHBB.mes.reclust@meta.data[which(all_int_final.noHBB.mes.reclust$celltype == "SMC"),])
Fibros <- rownames(all_int_final.noHBB.mes.reclust@meta.data[which(all_int_final.noHBB.mes.reclust$celltype == "Fibroblast"),])
all_int_final.noHBB@meta.data$celltype[which(rownames(all_int_final.noHBB@meta.data) %in% Endos)] <- "Endothelial"
all_int_final.noHBB@meta.data$celltype[which(rownames(all_int_final.noHBB@meta.data) %in% SMCs)] <- "SMC"
all_int_final.noHBB@meta.data$celltype[which(rownames(all_int_final.noHBB@meta.data) %in% Fibros)] <- "Fibroblast"


all_int_final.noHBB@meta.data$celltype_major <- NA
all_int_final.noHBB@meta.data$celltype <- NA
#all commented because we're going to label once we get up to higher res on subsets
#Mesenchymal 
all_int_final.noHBB@meta.data$celltype_major[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(0,1,2,3,5,6,11,12,16))] <- "Mesenchymal"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(0,5,12,16))] <- "Endothelial"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(1,6))] <- "Myofibroblast"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(11))] <- "Fibroblast"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(2,3))] <- "SMC"

#Immune 
all_int_final.noHBB@meta.data$celltype_major[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(4,7,8,9,10,14,17))] <- "Immune"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(4))] <- "T"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(7,10,17))] <- "Myeloid"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(8,9))] <- "NK"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(14))] <- "B"

#Contaminating
all_int_final.noHBB@meta.data$celltype_major[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(13,15))] <- "Contaminating"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(13))] <- "Epithelial"
all_int_final.noHBB@meta.data$celltype[which(all_int_final.noHBB@meta.data$seurat_clusters %in% c(15))] <- "Erythroid"

#plot the new cell populations
UMAPPlot(all_int_final.noHBB, split.by="celltype_major", cols=c("#D84CA3",  "#DC8B47", "#84BDD7")) + theme_classic()
UMAPPlot(all_int_final.noHBB, split.by="celltype_major") + theme_classic()
UMAPPlot(all_int_final.noHBB, group.by="celltype") + theme_classic()
NKvsT <- FindMarkers(all_int_final.noHBB, ident.1 = 8, ident.2 = 9)

### celltype_major subsets ####
all_int_final.noHBB.mes <- subset(all_int_final.noHBB, subset=celltype_major=="Mesenchymal")
# saveRDS(all_int_final.noHBB.mes, "saves/all_int_final_noHBB_mesenchymal.rds")
all_int_final.nohbbcells.imm <- subset(all_int_final.nohbbells, subset=celltype_major=="Immune")
# saveRDS(all_int_final.noHBB.imm, "saves/all_int_final_noHBB_immune.rds")

### Mesenchymal subset ####
nohbbcells <- TRUE
all_int_final.noHBBcells.mes <- subset(all_int_final.noHBBcells, subset = celltype_major == "Mesenchymal")
all_int_final.noHBBcells.mes <- subset(all_int_final.noHBBcells.mes, 
                                       features = rownames(all_int_final.noHBBcells.mes@assays$RNA)[grep("^HB[AB][12]*", rownames(all_int_final.noHBBcells.mes@assays$RNA), invert = T)])
# saveRDS(possibleNKcells, "saves/possibleNKcells.rds")
possibleNKcells <- readRDS("saves/possibleNKcells.rds")
not_possibleNKcells <- rownames(all_int_final.noHBBcells.mes@meta.data[which(!(rownames(all_int_final.noHBBcells.mes@meta.data) %in% possibleNKcells)),])
all_int_final.noHBBcells.mes <- subset(all_int_final.noHBBcells.mes, cells = not_possibleNKcells) #added 20200225 to remove probable NK cells (only 73)
all_int_final.noHBBcells.mes <- FindVariableFeatures(all_int_final.noHBBcells.mes, nfeatures = 1000)
all_int_final.noHBBcells.mes <- RunPCA(all_int_final.noHBBcells.mes, features = VariableFeatures(object = all_int_final.noHBBcells.mes))
all_int_final.noHBBcells.mes <- RunUMAP(all_int_final.noHBBcells.mes, dims=1:20)
all_int_final.noHBBcells.mes <- FindNeighbors(all_int_final.noHBBcells.mes)
all_int_final.noHBBcells.mes <- FindClusters(all_int_final.noHBBcells.mes, resolution = 0.4) #this isn't working exactly as I want, but it'll do for now 
UMAPPlot(all_int_final.noHBBcells.mes, group.by="condition")
ggsave("nohbbcells/mes/all_int_final_noHBBcells_mes_umap_byCondition.pdf")
UMAPPlot(all_int_final.noHBBcells.mes, label=TRUE)
ggsave("nohbbcells/mes/all_int_final_noHBBcells_mes_umap_byCluster.pdf")
FeaturePlot(all_int_final.noHBBcells.mes, features=c("HBB", "HBA2", "HBA1"))
ggsave("nohbbcells/mes/all_int_final_noHBBcells_mes_hbb_features.pdf")
all_int_final.noHBBcells.mes.markers <- FindAllMarkers(all_int_final.noHBBcells.mes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
height80dpi <- (length(all_int_final.noHBBcells.mes.markers$gene)*12)*(1/80)+1
all_int_final.noHBBcells.mes.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- top10$gene
DoHeatmap(all_int_final.noHBBcells.mes, features = top10genes) #do this
ggsave("nohbbcells/mes/all_int_final_noHBBcells_mes_heatmap_res0.4_widertop10.pdf", height=20, width=14)
umap_withPies(all_int_final.noHBBcells.mes, pie_cond = "condition")
write.table(all_int_final.noHBBcells.mes.markers, "nohbbcells/mes/all_int_final_HBBcellsremoved_mesenchymal.txt", sep = "\t", quote = F, col.names = NA)

# saveRDS(all_int_final.noHBBcells.mes, "saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.rds") #this is the version where there's a lot of endothelial cells
all_int_final.noHBBcells.mes <- readRDS("saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.rds") #this is the version where there's a lot of endothelial cells
# all_int_final.noHBB.mes.reclust <- readRDS("saves/all_int_final_noHBB_mes_reclustered_v3.rds") #this is the version where there's a lot of endothelial cells

dir.create("nohbbcells/mes/features")
for (f in top10genes) {
  message(f)
  FeaturePlot(all_int_final.noHBBcells.mes, features = f)
  ggsave(paste0("nohbbcells/mes/features/",f,".pdf"))
}


### UPDATE 20200204 # CELL LABELING ISSUE FIX
all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters <- as.character(all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters)
all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters[which(all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters %in% c("1.0", "1.1"))] <- 1
all_int_final.noHBB.mes.reclust@meta.data$celltype[which(all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters %in% c(1, 2, 10))] <- "SMC"
all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters[which(all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters %in% c("7.0", "7.1"))] <- 7
all_int_final.noHBB.mes.reclust@meta.data$celltype[which(all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters %in% c(7,8))] <- "Fibroblast"
all_int_final.noHBB.mes.reclust@meta.data$celltype[which(all_int_final.noHBB.mes.reclust@meta.data$seurat_clusters %in% c("0.0"))] <- "Endothelial"
# saveRDS(all_int_final.noHBB.mes.reclust, "saves/all_int_final_noHBB_mes_reclustered_v4.rds") #this is the version where the smcs and fibroblast cells are actually correct...
fibros <- rownames(all_int_final.noHBB.mes.reclust@meta.data[which(all_int_final.noHBB.mes.reclust@meta.data$celltype == "Fibroblast"),])
smc <- rownames(all_int_final.noHBB.mes.reclust@meta.data[which(all_int_final.noHBB.mes.reclust@meta.data$celltype == "SMC"),])
endo <- rownames(all_int_final.noHBB.mes.reclust@meta.data[which(all_int_final.noHBB.mes.reclust@meta.data$celltype == "Endothelial"),])
all_int_final.noHBB@meta.data$celltype[which(rownames(all_int_final.noHBB@meta.data) %in% fibros)] <- "Fibroblast"
all_int_final.noHBB@meta.data$celltype[which(rownames(all_int_final.noHBB@meta.data) %in% smc)] <- "SMC"
all_int_final.noHBB@meta.data$celltype[which(rownames(all_int_final.noHBB@meta.data) %in% endo)] <- "Endothelial"
# saveRDS(all_int_final.noHBB, "saves/all_int_final_noHBB_v3.rds") #this is the version where the smcs and fibroblast cells are actually correct...


### Add clusters to mesenchymal subset ###
all_int_final.noHBB.mes <- readRDS("saves/all_int_final_noHBB_mes_reclustered.rds") #mesenchymal
seurat_object <- all_int_final.noHBB.mes
new_clusters <- c(0,1,2,7)
seurat_object@meta.data$seurat_clusters <- as.character(seurat_object@meta.data$seurat_clusters)
for (nc in new_clusters) {
  seurat_object.sub <- subset(seurat_object, subset = seurat_clusters == nc) # an existing cluster that will be reclustered  
  seurat_object.sub <- FindVariableFeatures(seurat_object.sub, nfeatures = 500)
  seurat_object.sub <- RunPCA(seurat_object.sub, features = VariableFeatures(object = seurat_object.sub))
  seurat_object.sub <- RunUMAP(seurat_object.sub, dims=1:20)
  seurat_object.sub <- FindNeighbors(seurat_object.sub)
  seurat_object.sub <- FindClusters(seurat_object.sub, resolution = 0.1)
  for (sc in unique(seurat_object.sub@meta.data$seurat_clusters)) {
    print(paste0("doing reclustering on cluster: ", sc, " in ", nc))
    cur_cells <- rownames(seurat_object.sub@meta.data[which(seurat_object.sub@meta.data$seurat_clusters == sc),])
    print(length(cur_cells))
    print(head(seurat_object@meta.data[which(rownames(seurat_object@meta.data) %in% cur_cells),]))
    seurat_object@meta.data[which(rownames(seurat_object@meta.data) %in% cur_cells), "seurat_clusters"] <- paste(nc, sc, sep = ".")
    print(head(seurat_object@meta.data[which(rownames(seurat_object@meta.data) %in% cur_cells),]))
  }
}
seurat_object@meta.data$seurat_clusters[which(seurat_object@meta.data$seurat_clusters == "2.0")] <- 2
seurat_object@meta.data$seurat_clusters <- factor(seurat_object@meta.data$seurat_clusters)
seurat_object <- SetIdent(seurat_object, value="seurat_clusters")

all_markers <- c()
for (nc in c("0.0", "0.1", "0.2", "1.0", "1.1", "7.0", "7.1")) {
  if (nc == 0.0) {
    all_markers <- FindMarkers(seurat_object, ident.1 = "0.0")
    all_markers$cluster <- "0.0"
  } else {
    curmarkers <- FindMarkers(seurat_object, ident.1 = nc)
    curmarkers$cluster <- nc
    all_markers <- rbind(all_markers, curmarkers)
  }
}
write.table(all_markers, "all_int_final_noHBB_mes_reclust_v2_markers.txt", sep = "\t", quote = F, col.names = NA)

# seurat_object@meta.data$seurat_clusters <- factor(seurat_object@meta.data$seurat_clusters, levels=c(0.0, 0.1, 0.2, 1.0, 1.1, 2, 3, 4, 5, 6, 7.0, 7.1, 8, 9, 10))
# all_int_final.noHBB.mes.reclust.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_int_final.noHBB.mes.reclust <- seurat_object
# saveRDS(all_int_final.noHBB.mes.reclust, "saves/all_int_final_noHBB_mes_reclustered_v2.rds") #mesenchymal
DoHeatmap(all_int_final.noHBB.mes.reclust, features = all_int_final.noHBB.mes.reclust.markers$gene) #do this
ggsave("all_int_final_noHBB_mes_reclust_v2_heatmap.pdf", height=30, width=12)

all_int_final.noHBB.mes.reclust.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- top10$gene
DoHeatmap(all_int_final.noHBB.mes.reclust, features = top10genes) #do this
ggsave("all_int_final_noHBB_mes_reclust_v2_top10_heatmap.pdf", height=35, width=12)

# 
umap <- seurat_object@reductions$umap@cell.embeddings 
umap.withmeta <- merge(umap, seurat_object@meta.data, by='row.names')
rownames(umap.withmeta) <- umap.withmeta[,1]
umap.withmeta <- umap.withmeta[,-1]
ggplot(data=umap.withmeta, aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size = 0.6) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank(),
        text = element_text(face="bold",size=18)) 

all_int_final.noHBB.mes.reclust@meta.data$celltype[which(all_int_final.noHBB.mes.reclust@meta.data$celltype == "Myofibroblast")] <- "SMC"


#### Endothelial cells ####
all_int_final.noHBBcells.mes.endo <- subset(all_int_final.noHBBcells.mes, subset = celltype == "Endothelial")
all_int_final.noHBBcells.mes.endo <- FindVariableFeatures(all_int_final.noHBBcells.mes.endo, nfeatures = 1000)
all_int_final.noHBBcells.mes.endo <- RunPCA(all_int_final.noHBBcells.mes.endo, features = VariableFeatures(object = all_int_final.noHBBcells.mes.endo))
all_int_final.noHBBcells.mes.endo <- RunUMAP(all_int_final.noHBBcells.mes.endo, dims=1:20)
all_int_final.noHBBcells.mes.endo@meta.data$seurat_clusters_super <- all_int_final.noHBBcells.mes.endo@meta.data$seurat_clusters
all_int_final.noHBBcells.mes.endo <- FindNeighbors(all_int_final.noHBBcells.mes.endo)
all_int_final.noHBBcells.mes.endo <- FindClusters(all_int_final.noHBBcells.mes.endo, resolution = 0.4) #this isn't working exactly as I want, but it'll do for now 
x <- UMAPPlot(all_int_final.noHBBcells.mes.endo, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) +
  guides(color=guide_legend(ncol=2, override.aes = list(size=10)))
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 4
z <- ggplot_gtable(y)
pdf("nohbbcells/mes/endo/all_int_final_noHBBcells_mes_endo_umap_byNewClustering.pdf", width=6,height = 5)
plot(z) 
dev.off()

x <- UMAPPlot(all_int_final.noHBBcells.mes.endo, pt.size = 1, group.by="condition") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white")) +
  ggtitle("by old cluster")
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 7
z <- ggplot_gtable(y)
pdf("nohbbcells/mes/endo/all_int_final_noHBBcells_mes_endo_umap_byCondition.pdf")
plot(z) 
dev.off()
all_int_final.noHBBcells.mes.endo.markers <- FindAllMarkers(all_int_final.noHBBcells.mes.endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBBcells.mes.endo.markers, "nohbbcells/mes/endo/all_int_final_noHBBcells_mes_endo_markers.txt", 
            sep = "\t", quote = F, col.names = NA)
all_int_final.noHBBcells.mes.endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.noHBBcells.mes.endo, features = top10genes) #do this
height80dpi <- (length(top10genes)*12)*(1/80)+1
ggsave("nohbbcells/mes/endo/all_int_final.noHBBcells.mes.endo_heatmap_res0.4.pdf", height=height80dpi, width=8)
saveRDS(all_int_final.noHBBcells.mes.endo, "saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.endo.rds") #this is the version where there's a lot of endothelial cells
all_int_final.noHBBcells.mes.endo <- readRDS("saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.endo.rds") #this is the version where there's a lot of endothelial cells

### Smooth muscle cells ####
all_int_final.noHBBcells.mes.smc <- subset(all_int_final.noHBBcells.mes.copy, subset = celltype == "SMC") #make SURE  that these SMC cells are actually the 1,2,10 ones.... (they are)
all_int_final.noHBBcells.mes.smc <- FindVariableFeatures(all_int_final.noHBBcells.mes.smc, nfeatures = 1000)
all_int_final.noHBBcells.mes.smc <- RunPCA(all_int_final.noHBBcells.mes.smc, features = VariableFeatures(object = all_int_final.noHBBcells.mes.smc))
all_int_final.noHBBcells.mes.smc <- RunUMAP(all_int_final.noHBBcells.mes.smc, dims=1:20)
# all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters_super <- all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters
all_int_final.noHBBcells.mes.smc <- FindNeighbors(all_int_final.noHBBcells.mes.smc)
all_int_final.noHBBcells.mes.smc <- FindClusters(all_int_final.noHBBcells.mes.smc, resolution = 0.4) #this isn't working exactly as I want, but it'll do for now 
saveRDS(all_int_final.noHBBcells.mes.smc, "saves/all_int_final.noHBBcells_mes_smc_afterMovingAfewCellsfromFibroblast.rds")
x <- UMAPPlot(all_int_final.noHBBcells.mes.smc, pt.size = 1, group.by="seurat_clusters") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white")) +
  ggtitle("by old cluster")
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 7
z <- ggplot_gtable(y)
pdf("nohbbcells/mes/smc/all_int_final_noHBBcells_mes_smc_umap_byOldClustering.pdf")
plot(z) 
dev.off()
x <- UMAPPlot(all_int_final.noHBBcells.mes.smc, pt.size = 1, group.by="orig.ident") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white")) +
  ggtitle("by old cluster")
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 7
z <- ggplot_gtable(y)
pdf("nohbbcells/mes/smc/all_int_final_noHBBcells_mes_smc_umap_bySample.pdf")
plot(z) 
dev.off()
all_int_final.noHBBcells.mes.smc.markers <- FindAllMarkers(all_int_final.noHBBcells.mes.smc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBBcells.mes.smc.markers, "nohbbcells/mes/smc_after_adding_those_previous_fibroblast_cells/all_int_final_noHBBcells_mes_smc_markers.txt", 
            sep = "\t", quote = F, col.names = NA)
all_int_final.noHBBcells.mes.smc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.noHBBcells.mes.smc, features = top10genes) #do this
height80dpi <- (length(top10genes)*12)*(1/80)+1
ggsave("nohbbcells/mes/smc_after_adding_those_previous_fibroblast_cells/all_int_final.noHBBcells.mes.smc_heatmap_res0.4.pdf", height=height80dpi, width=8)

# saveRDS(all_int_final.noHBBcells.mes.smc, "saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.smc.rds") #this is the version where there's a lot of endothelial cells
all_int_final.noHBBcells.mes.smc <- readRDS("saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.smc.rds") #this is the version where there's a lot of endothelial cells

#markers between big clusters
all_int_final.noHBB.mes.smc@meta.data$bigclusters <- 1
all_int_final.noHBB.mes.smc@meta.data$bigclusters[which(all_int_final.noHBB.mes.smc@meta.data$seurat_clusters %in% c(0,3,4,6))] <- 2
Idents(all_int_final.noHBB.mes.smc) <- all_int_final.noHBB.mes.smc@meta.data$bigclusters
all_int_final.noHBB.mes.smc.markers <- FindMarkers(all_int_final.noHBB.mes.smc, ident.1 = 1, ident.2 = 2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBB.mes.smc.markers, "FINAL_merge/hbb_regressed/mesenchymal/smc/all_int_final_noHBB_mes_smc_markers_btwn_bigclusters.txt",
            sep = "\t", quote = F, col.names = NA)
all_int_final.noHBB.mes.smc.markers <- all_int_final.noHBB.mes.smc.markers[order(all_int_final.noHBB.mes.smc.markers$avg_logFC, decreasing = T),]
top10genes <- head(as.vector(rownames(all_int_final.noHBB.mes.smc.markers)), n=20)
DoHeatmap(all_int_final.noHBB.mes.smc, features = top10genes) #do this
ggsave("FINAL_merge/hbb_regressed/mesenchymal/smc/all_int_final_noHBB_mes_smc_heatmap_bigcluster1vs2.pdf", height=30, width=8)



## Fibroblast cells ####
all_int_final.noHBBcells.mes.fib <- subset(all_int_final.noHBBcells.mes, celltype == "Fibroblast")
all_int_final.noHBBcells.mes.fib <- FindVariableFeatures(all_int_final.noHBBcells.mes.fib, nfeatures = 1000)
all_int_final.noHBBcells.mes.fib <- RunPCA(all_int_final.noHBBcells.mes.fib, features = VariableFeatures(object = all_int_final.noHBBcells.mes.fib))
all_int_final.noHBBcells.mes.fib <- RunUMAP(all_int_final.noHBBcells.mes.fib, dims=1:20)
all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters_super <- all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters
all_int_final.noHBBcells.mes.fib <- FindNeighbors(all_int_final.noHBBcells.mes.fib)
all_int_final.noHBBcells.mes.fib <- FindClusters(all_int_final.noHBBcells.mes.fib, resolution = 0.4) #this isn't working exactly as I want, but it'll do for now 
# saveRDS(all_int_final.noHBBcells.mes.fib, "saves/all_int_final_noHBBcells_mes_fib.rds")
# all_int_final.noHBBcells.mes.fib <- readRDS("saves/all_int_final_noHBBcells_mes_fib.rds")
x <- UMAPPlot(all_int_final.noHBBcells.mes.fib, pt.size = 1, group.by="seurat_clusters") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white")) 
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 7
z <- ggplot_gtable(y)
pdf("nohbbcells/mes/fibroblast/all_int_final_noHBBcells_mes_fib_umap_byNewClustering.pdf")
plot(z)
dev.off()
x <- UMAPPlot(all_int_final.noHBBcells.mes.fib, pt.size = 1, group.by="condition") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white"))
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 7
z <- ggplot_gtable(y)
pdf("nohbbcells/mes/fibroblast/all_int_final_noHBBcells_mes_fib_umap_byCondition.pdf")
plot(z)
dev.off()
all_int_final.noHBBcells.mes.fib.markers <- FindAllMarkers(all_int_final.noHBBcells.mes.fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBBcells.mes.fib.markers, "nohbbcells/mes/fibroblast/all_int_final_noHBBcells_mes_fib_markers.txt",
            sep = "\t", quote = F, col.names = NA)
all_int_final.noHBBcells.mes.fib.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
height80dpi <- (length(top10genes)*12)*(1/80)+1
DoHeatmap(all_int_final.noHBBcells.mes.fib, features = top10genes) #do this
ggsave("nohbbcells/mes/fibroblast_final/all_int_final_noHBBcells_mes_fib_heatmap_res0.4.pdf", height=height80dpi, width=8)
#markers between big clusters

# saveRDS(all_int_final.noHBBcells.mes.fib, "saves/all_int_final.noHBBcells_noHBBgenes.noNK.mes.fib.rds") #this is the version where there's a lot of endothelial cells
all_int_final.noHBBcells.mes.fib <- readRDS("saves/all_int_final.noHBBcells_noHBBgenes.noNK.mes.fib.rds") #this is the version where there's a lot of endothelial cells

### Immune subset ####
all_int_final.nohbbells.imm <- FindVariableFeatures(all_int_final.nohbbells.imm, nfeatures = 1000)
all_int_final.nohbbells.imm <- RunPCA(all_int_final.nohbbells.imm, features = VariableFeatures(object = all_int_final.nohbbells.imm))
all_int_final.nohbbells.imm <- RunUMAP(all_int_final.nohbbells.imm, dims=1:20)
all_int_final.nohbbells.imm <- FindNeighbors(all_int_final.nohbbells.imm)
all_int_final.nohbbells.imm <- FindClusters(all_int_final.nohbbells.imm, resolution = 0.4) #this isn't working exactly as I want, but it'll do for now 
UMAPPlot(all_int_final.nohbbells.imm, group.by="condition")
ggsave("all_int_final_noHBB_imm_umap_byCondition.pdf")
UMAPPlot(all_int_final.nohbbells.imm, label=TRUE)
ggsave("all_int_final_noHBB_imm_umap_byCluster.pdf")
all_int_final.noHBB.imm.allmarkers <- FindAllMarkers(all_int_final.noHBB.imm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBB.imm.allmarkers, "all_int_final_noHBB_imm_markers.txt", sep = "\t", quote = F, col.names = NA)
all_int_final.noHBB.imm.allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- top10$gene
DoHeatmap(all_int_final.noHBB.imm, features = all_int_final.noHBB.imm.allmarkers$gene) #do this

height80dpi <- (length(top10genes)*12)*(1/80)+1 ##formula for size is (num_genes*12ptFont)*(1inch/80dpi) = height   ---- 
ggsave("FINAL_merge/hbb_regressed/immune/all_int_final_noHBB_imm_heatmap_res0.4_smaller.pdf", height=height80dpi, width=8)
umap_withPies(all_int_final.noHBB.imm, pie_cond = "condition")


### Add clusters to immune subset ###
all_int_final.noHBB.imm <- readRDS("saves/all_int_final_noHBB_imm_reclustered.rds") #mesenchymal
all_int_final.noHBB.imm@meta.data$orig.ident[grep("fibroid55_12746_lane1", all_int_final.noHBB.imm@meta.data$orig.ident)] <- "fibroid55_12746"
all_int_final.noHBB.imm@meta.data$orig.ident[grep("fibroid55_12746_lane2", all_int_final.noHBB.imm@meta.data$orig.ident)] <- "fibroid55_12746"
all_int_final.noHBB.imm@meta.data$orig.ident[grep("fibroid55_12906lane2", all_int_final.noHBB.imm@meta.data$orig.ident)] <- "fibroid55_12906"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(0))] <- "Cytotoxic T"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(1))] <- "T"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(2))] <- "NK"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(3))] <- "CD8+ T"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(4,6))] <- "Macrophages"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(5))] <- "CD14 Monocotyes"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(7))] <- "FCGR3A+ Monocytes"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(8))] <- "B"
all_int_final.noHBB.imm@meta.data$celltype[which(all_int_final.noHBB.imm@meta.data$seurat_clusters %in% c(9))] <- "Dendritic"
# saveRDS(all_int_final.noHBB.imm, "saves/all_int_final_noHBB_imm_reclustered_v2.rds")
# all_int_final.noHBB.imm <- readRDS("saves/all_int_final_noHBB_imm_reclustered_v2.rds")

seurat_object <- all_int_final.noHBB.mes
new_clusters <- c(0,1,2,7)
seurat_object@meta.data$seurat_clusters <- as.character(seurat_object@meta.data$seurat_clusters)
for (nc in new_clusters) {
  seurat_object.sub <- subset(seurat_object, subset = seurat_clusters == nc) # an existing cluster that will be reclustered  
  seurat_object.sub <- FindVariableFeatures(seurat_object.sub, nfeatures = 500)
  seurat_object.sub <- RunPCA(seurat_object.sub, features = VariableFeatures(object = seurat_object.sub))
  seurat_object.sub <- RunUMAP(seurat_object.sub, dims=1:20)
  seurat_object.sub <- FindNeighbors(seurat_object.sub)
  seurat_object.sub <- FindClusters(seurat_object.sub, resolution = 0.1)
  for (sc in unique(seurat_object.sub@meta.data$seurat_clusters)) {
    print(paste0("doing reclustering on cluster: ", sc, " in ", nc))
    cur_cells <- rownames(seurat_object.sub@meta.data[which(seurat_object.sub@meta.data$seurat_clusters == sc),])
    print(length(cur_cells))
    print(head(seurat_object@meta.data[which(rownames(seurat_object@meta.data) %in% cur_cells),]))
    seurat_object@meta.data[which(rownames(seurat_object@meta.data) %in% cur_cells), "seurat_clusters"] <- paste(nc, sc, sep = ".")
    print(head(seurat_object@meta.data[which(rownames(seurat_object@meta.data) %in% cur_cells),]))
  }
}
seurat_object@meta.data$seurat_clusters <- factor(seurat_object@meta.data$seurat_clusters, levels=c(0.0, 0.1, 0.2, 1.0, 1.1, 2.0, 2.1, 3, 4, 5, 6, 7.0, 7.1, 8, 9, 10))
seurat_object <- SetIdent(seurat_object, value="seurat_clusters")
all_int_final.noHBB.mes.reclust.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_int_final.noHBB.mes.reclust.markers, "all_int_final_noHBB_mes_reclust_v2_markers.txt", sep = "\t", quote = F, col.names = NA)
all_int_final.noHBB.mes.reclust <- seurat_object
# saveRDS(all_int_final.noHBB.mes.reclust, "saves/all_int_final_noHBB_mes_reclustered_v2.rds") #mesenchymal
DoHeatmap(all_int_final.noHBB.mes.reclust, features = all_int_final.noHBB.mes.reclust.markers$gene) #do this
ggsave("all_int_final_noHBB_mes_reclust_v2_heatmap.pdf", height=40, width=8)

# saveRDS(all_int_final.noHBB.mes, "all_int_final_noHBB_mes_reclustered.rds")
# saveRDS(all_int_final.noHBB.imm, "all_int_final_noHBB_imm_reclustered.rds")

#### EASY DE ####
# makes a volcano plot for each comparison in each cluster for a given seurat object
# seurat_object <- all_int_final
seurat_object <- all_int_final.nohbbells.imm
clusters <- sort(unique(seurat_object@meta.data$seurat_clusters))
for (c in clusters) {
  clusterfolder <- paste0("cluster_", c, "_DE")
  dir.create(clusterfolder)
  
  cluster.cells <- rownames(seurat_object@meta.data[which(seurat_object@meta.data$seurat_clusters == c),])
  cluster.set <- subset(seurat_object, cells = cluster.cells) #only cells in this cluster
  for (comp in c("med12pos_vs_myometrium", "med12pos_vs_med12neg", "med12neg_vs_myometrium")) { #
    dir.create(paste0("cluster_", c, "_DE/", comp))
    
    splitted <- unlist(strsplit(comp, "_vs_"))
    cond1 <- splitted[1]
    cond2 <- splitted[2]
    cond1.cells <- rownames(cluster.set@meta.data[which(cluster.set@meta.data$condition == cond1),])
    cond2.cells <- rownames(cluster.set@meta.data[which(cluster.set@meta.data$condition == cond2),])
    if (length(cond1.cells) <= 3| length(cond2.cells) <=3 ) {
      print(paste0("comp ", comp, " in cluster ", c, " has too few cells. Skipping"))
      next
    }
    
    cluster.set <- SetIdent(cluster.set, cells = cond1.cells, value=cond1)
    cluster.set <- SetIdent(cluster.set, cells = cond2.cells, value=cond2)
    cur.de <- FindMarkers(cluster.set, ident.1 = cond1, ident.2 = cond2, min.pct = 0.25, logfc.threshold = 0.25) #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
    write.table(cur.de, paste0(clusterfolder, "/", comp, "/DEmarkers_", c, "_", comp, ".txt"), quote = F, sep = "\t", col.names = NA)
    
    # volcano plot ----
    cur.de$logPadjf <- -log10(cur.de$p_val_adj + 1e-300)
    volcdat <- cur.de[,grep("avg_logFC|logPadjf", colnames(cur.de))]
    cur.de$labels <- NA
    topLeft <- c(-10,300)
    topRight <- c(10,300)
    for (i in 1:nrow(cur.de)) {
      lDist <- dist(rbind(topLeft, volcdat[i,]), method = "euclidean")
      rDist <- dist(rbind(topRight, volcdat[i,]), method = "euclidean")
      cur.de$lDist[i] <- lDist
      cur.de$rDist[i] <- rDist
    }
    closestRight = cur.de[order(cur.de$lDist, decreasing = F),]
    closestLeft = cur.de[order(cur.de$lDist, decreasing = F),]
    closestRight.genes <- rownames(head(closestRight, n=8))
    closestLeft.genes <- rownames(head(closestLeft, n=8))
    
    if (nrow(cur.de) < 15) {
      cur.de$labels <- rownames(cur.de)
    } else {
      cur.de$labels[which(rownames(cur.de) %in% c(closestRight.genes, closestLeft.genes))] <- rownames(cur.de[which(rownames(cur.de) %in% c(closestRight.genes, closestLeft.genes)),])
    }
    
    cur.de$cumulative_percent_of_cluster <- (cur.de$pct.2 + cur.de$pct.1)/2
    r <- max(abs(cur.de$avg_logFC)) + 1
    volcano <- ggplot2::ggplot(data=cur.de, aes(x = avg_logFC, y = logPadjf, size=cumulative_percent_of_cluster)) + #not yet sure how this should be colored
      geom_point(aes(alpha=0.5, color=cumulative_percent_of_cluster)) + 
      scale_color_viridis() + 
      geom_text_repel(aes(label = labels), point.padding = 0.25, size = 3) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.5, alpha = 0.4) +
      geom_vline(xintercept=c(-1,1), linetype="dashed", size=0.5, alpha = 0.4) +
      xlab(paste0("log2 Fold change (", cond1[1], "/", cond2[2], ")")) + scale_x_continuous(limits = c(-r, r)) +
      ylab("-log10 adjusted p-value") +
      theme(legend.key.size = unit(2, 'lines')) +
      ggtitle(paste0(comp))+
      theme_bw() + theme(text = element_text(size=10),
                         legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'),
                         legend.position = "right")
    pdf(paste0(clusterfolder, "/", comp, "/volcano_", c, "_", comp, ".pdf")) #width = 720, height=720, units='px', type = "cairo" (png)
    print(volcano)
    dev.off()
  }  
}


#### MYO subsets ####
#before subsetting smc, need to make adjustment to clusters (oops)
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters == 2)] <- 1
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters <- factor(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters) 
Idents(all_int_final.noHBBcells.mes.smc) <- all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters

setwd("nohbbcells/mes/myometrium_ONLY/")
all_int_final.nohbbcells.mes.fib.myo <- subset(all_int_final.noHBBcells.mes.fib, subset=condition == "myometrium")
all_int_final.nohbbcells.mes.smc.myo <- subset(all_int_final.noHBBcells.mes.smc, subset=condition == "myometrium")
all_int_final.nohbbcells.mes.endo.myo <- subset(all_int_final.noHBBcells.mes.endo, subset=condition == "myometrium")
all_int_final.nohbbcells.mes.fib.myo@project.name <- "all_int_final.nohbbcells.mes.fib.myo"
all_int_final.nohbbcells.mes.smc.myo@project.name <- "all_int_final.nohbbcells.mes.smc.myo"
all_int_final.nohbbcells.mes.endo.myo@project.name <- "all_int_final.nohbbcells.mes.endo.myo"
saveRDS(all_int_final.nohbbcells.mes.fib.myo, "../../../saves/all_int_final.nohbbcells.mes.fib.myo.rds")
saveRDS(all_int_final.nohbbcells.mes.smc.myo, "../../../saves/all_int_final.nohbbcells.mes.smc.myo.rds")
saveRDS(all_int_final.nohbbcells.mes.endo.myo, "../../../saves/all_int_final.nohbbcells.mes.endo.myo.rds")

myo.objects <- list(all_int_final.nohbbcells.mes.fib.myo, all_int_final.nohbbcells.mes.smc.myo, all_int_final.nohbbcells.mes.endo.myo)
for (o in myo.objects) {
  actual.name <- strsplit(o@project.name, "\\.")[[1]][4]
  cur.dir <- paste0(getwd(), "/", actual.name, "/")
  dir.create(cur.dir)
  
  #clustering 
  UMAPPlot(o)
  ggsave(paste0(cur.dir,"/umap_byCluster.pdf"))
  UMAPPlot(o, group.by="condition")
  ggsave(paste0(cur.dir,"/umap_byCondition.pdf"))
  UMAPPlot(o, group.by="orig.ident")
  ggsave(paste0(cur.dir,"/umap_bySample.pdf"))
  
  #markers and heatmap
  cur.markers <- FindAllMarkers(o,  logfc.threshold = 0.25, only.pos = TRUE) #assay = RNA
  write.table(cur.markers, paste0(cur.dir,"/", o@project.name, "_markers.txt"), quote = F, sep = "\t", col.names = NA)
  height80dpi <- (length(cur.markers$gene)*12)*(1/80)+1
  cur.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
  top10genes <- top10$gene
  DoHeatmap(o, features = top10genes) #do this
  ggsave(paste0(cur.dir, "/heatmap.pdf"), height=20, width=14)
  
}


COLORSCHEME2 = c("#3E80A0", "#5E3EA0", "#A03E80", "#3ea05e", "#a05e3e", "#A1903F", "#a23f3f" ,"#3fa27a", "#a23f68")
DARKBLUE = "#000099"
all.myo.jg <- readRDS("saves/all_int_final_myometriumonly_JG.rds") #wrong one
all_int_final.nohbbcells <- readRDS("saves/all_int_final.nohbbells_orEpithelial.rds")
all.myo <- subset(all_int_final.nohbbcells, subset = condition == "myometrium")
ident.names <- c("Endothelial", "Endothelial", "Smooth Muscle", "Smooth Muscle", "T", "Endothelial", "Lymphatic Endothelial", "Myeloid", "NK", "NK", "Myeloid", "Fibroblast",
                 "Endothelial", "Epithelial", "B", "Erythroid", "Endothelial", "Myeloid")

# load object, change celltype to above ident.names, remove epi and erythro, make object, then subset on myo only and make figures
all_int_final_noHBB_v3 <- readRDS("saves/all_int_final_noHBB_v3.rds") 
all_int_final_noHBB_v3@meta.data$celltype <- NA
for (i in 1:length(unique(Idents(all_int_final_noHBB_v3)))) {
  cluster = i-1
  all_int_final_noHBB_v3@meta.data$celltype[which(all_int_final_noHBB_v3@meta.data$seurat_clusters == cluster)] <- ident.names[i]
}
all_int_final_noHBB_v3 <- subset(all_int_final_noHBB_v3, subset = celltype != "Epithelial" & celltype != "Erythroid")


#### MED12+ subsets####
all_int_final.noHBBcells.mes.fib <- readRDS("saves/all_int_final.nohbbcells.mes.fib.onlyclusters6and11.rds")
all_int_final.noHBBcells.mes.smc <- readRDS("saves/all_int_final.noHBBcells_mes_smc.rds")
all_int_final.noHBBcells.mes.endo <- readRDS("saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.endo.rds")
all_int_final.noHBBcells.mes <- readRDS("saves/all_int_final_nohbbcells_noNK_mes_adjusted_by_moving_some_fibrobasts_intoSMC.rds")

#before subsetting smc, need to make adjustment to clusters 
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters == 2)] <- 1
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters <- factor(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters) 
Idents(all_int_final.noHBBcells.mes.smc) <- all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters

setwd("nohbbcells/mes/med12pos_ONLY/")
all_int_final.nohbbcells.mes.fib.med12p <- subset(all_int_final.noHBBcells.mes.fib, subset=condition == "med12pos")
all_int_final.nohbbcells.mes.smc.med12p <- subset(all_int_final.noHBBcells.mes.smc, subset=condition == "med12pos")
all_int_final.nohbbcells.mes.endo.med12p <- subset(all_int_final.noHBBcells.mes.endo, subset=condition == "med12pos")
all_int_final.noHBBcells.mes.med12p <- subset(all_int_final.noHBBcells.mes, subset=condition == "med12pos")

all_int_final.nohbbcells.mes.fib.med12p@project.name <- "all_int_final.nohbbcells.mes.fib.med12p"
all_int_final.nohbbcells.mes.smc.med12p@project.name <- "all_int_final.nohbbcells.mes.smc.med12p"
all_int_final.nohbbcells.mes.endo.med12p@project.name <- "all_int_final.nohbbcells.mes.endo.med12p"
all_int_final.noHBBcells.mes.med12p@project.name <- "all_int_final.noHBBcells.mes.med12p"

saveRDS(all_int_final.nohbbcells.mes.fib.med12p, "../../../saves/all_int_final.nohbbcells.mes.fib.med12p.rds")
saveRDS(all_int_final.nohbbcells.mes.smc.med12p, "../../../saves/all_int_final.nohbbcells.mes.smc.med12p.rds")
saveRDS(all_int_final.nohbbcells.mes.endo.med12p, "../../../saves/all_int_final.nohbbcells.mes.endo.med12p.rds")
saveRDS(all_int_final.nohbbcells.mes.med12p, "../../../saves/all_int_final.nohbbcells.mes.med12p.rds")

med12p.objects <- list(all_int_final.nohbbcells.mes.fib.med12p, all_int_final.nohbbcells.mes.smc.med12p, all_int_final.nohbbcells.mes.endo.med12p)
for (o in med12p.objects) {
  actual.name <- strsplit(o@project.name, "\\.")[[1]][4]
  cur.dir <- paste0(getwd(), "/", actual.name, "/")
  dir.create(cur.dir)
  
  #clustering 
  UMAPPlot(o)
  ggsave(paste0(cur.dir,"/umap_byCluster.pdf"))
  UMAPPlot(o, group.by="condition")
  ggsave(paste0(cur.dir,"/umap_byCondition.pdf"))
  UMAPPlot(o, group.by="orig.ident")
  ggsave(paste0(cur.dir,"/umap_bySample.pdf"))
  
  #markers and heatmap
  cur.markers <- FindAllMarkers(o,  logfc.threshold = 0.25, only.pos = TRUE) #assay = RNA
  write.table(cur.markers, paste0(cur.dir,"/", o@project.name, "_markers.txt"), quote = F, sep = "\t", col.names = NA)
  height80dpi <- (length(cur.markers$gene)*12)*(1/80)+1
  cur.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
  top10genes <- top10$gene
  DoHeatmap(o, features = top10genes) #do this
  ggsave(paste0(cur.dir, "/heatmap.pdf"), height=20, width=14)
  
}

getCellProportionFigures(all_int_final.noHBBcells.mes.med12p)

#region to slice: chrX:71119403-71119405
##### CLONALITY with freemuxlet ####
### 12906_1 & 2 , 12843, 12882; AND 20200526_UPDATE: 12640 & 061919 #### "Fibroid12843",
variant.files <- c("Fibroid55-12906-1_GEXL",  "Fibroid55-12906-2_GEXL",
                   "Fibroid-55_12882_GEXL_analysis", "061919Fibroid_GEXL_fastqs", 
                   "Fibroid-55_12640_GEXL_analysis")
barcode_tables <- list()
processFreemuxlet <- function(cel.file, plp.file, var.file, name, refsnp = NULL) {
  
  cells <- read.table(cel.file,  sep = "\t")
  colnames(cells) <- c("DROPLET_ID",	"BARCODE",	"NUM.READ",	"NUM.UMI",	"NUM.UMIwSNP",	"NUM.SNP")
  
  pileup <- read.table(plp.file, sep = "\t")
  colnames(pileup) <- c("DROPLET_ID", "SNP_ID", "ALLELES", "BASEQS")
  
  variants <- read.table(var.file, sep = "\t")
  colnames(variants) <- c("SNP_ID", "CHROM", "POS", "REF", "ALT", "AF")
  
  clonality.table <- merge(cells, pileup, by="DROPLET_ID") #now we see which droplets have the variant according to freemuxlet
  clonality.table <- merge(clonality.table, variants, by="SNP_ID") #added info about snp
  if (!is.null(refsnp)) {
    clonality.table <- clonality.table[grep(refsnp, clonality.table[,grep("POS", colnames(clonality.table))]),]
  } 
  barcodes_with_alt <- gsub("-1", "", as.character(clonality.table$BARCODE))
  return(barcodes_with_alt)
}

all_int_final.noHBBcells.mes@meta.data$med12_variant <- "UNKNOWN"
all_int_final.noHBBcells.mes@meta.data$bcs <- 
  as.vector(paste0(vapply(strsplit(rownames(all_int_final.noHBBcells.mes@meta.data), "_"), 
                          "[", length(strsplit(rownames(all_int_final.noHBBcells.mes@meta.data), "_")[[1]]), 
                          FUN.VALUE = character(1)) ))
all_int_final.noHBBcells.mes@meta.data$bcs[grep("L00[1234]", all_int_final.noHBBcells.mes@meta.data$bcs)] <- 
  as.vector(paste0(vapply(strsplit(rownames(all_int_final.noHBBcells.mes@meta.data), "_"), 
                          "[", 4, FUN.VALUE = character(1)) ))[grep("L00[1234]", all_int_final.noHBBcells.mes@meta.data$bcs)]

for (i in 1:length(variant.files)) {
  curdir <- variant.files[i]
  num.start <- unlist(gregexpr('[0-9]+', curdir))[1] 
  if (substr(curdir, num.start, num.start+1) == '55') {
    snum <- substr(curdir, num.start+3, num.start+7)
  } else {
    snum <- substr(curdir, num.start, num.start+4)
  }
  print(curdir)
  
  cel.file <- paste0( curdir ,"/all.cel") #variantAnalysis/",
  var.file <- paste0( curdir ,"/all.var")
  plp.file <- paste0( curdir ,"/all.plp")
  barcodes_with_alt <- processFreemuxlet(cel.file, plp.file, var.file, name = curdir, refsnp = "71119405")
  print(paste0(length(barcodes_with_alt), " potential new bcs"))
  
  #add to existing meta data
  all_int_final.noHBBcells.mes@meta.data[
    which(all_int_final.noHBBcells.mes@meta.data$bcs %in% barcodes_with_alt & 
            grepl(snum, all_int_final.noHBBcells.mes@meta.data$orig.ident) &
            all_int_final.noHBBcells.mes@meta.data$condition == "med12pos"),"med12_variant"] <- curdir
  added_bcs <- length(all_int_final.noHBBcells.mes@meta.data[
    which(all_int_final.noHBBcells.mes@meta.data$bcs %in% barcodes_with_alt & 
            grepl(snum, all_int_final.noHBBcells.mes@meta.data$orig.ident) &
            all_int_final.noHBBcells.mes@meta.data$condition == "med12pos"),"med12_variant"])
  print(paste0("added: ", added_bcs, " matching barcodes"))
}
# UMAPPlot(all_int_final.noHBBcells.mes, group.by="med12_variant", cols = c(COLORSCHEME2[1:3], LIGHTGRAY), pt.size = 0.1)
cur.plotdata.umap <- all_int_final.noHBBcells.mes@reductions$umap@cell.embeddings
cur.plotdata.umap <- merge(cur.plotdata.umap, all_int_final.noHBBcells.mes@meta.data, by = "row.names")
rownames(cur.plotdata.umap) <- cur.plotdata.umap[,1]
cur.plotdata.umap <- cur.plotdata.umap[,-1]
cur.plotdata.umap$med12_variant <- factor(cur.plotdata.umap$med12_variant, levels = c("Fibroid-55_12882_GEXL_analysis",
                                                                                      "Fibroid55-12906-2_GEXL",
                                                                                      #"Fibroid12843",
                                                                                      "061919Fibroid_GEXL_fastqs",
                                                                                      "Fibroid-55_12640_GEXL_analysis",
                                                                                      "UNKNOWN"))

DefaultAssay(all_int_final.noHBBcells.mes) <- "RNA"
cellswithmed12 <- Cells(subset(all_int_final.noHBBcells.mes, subset = MED12 >= 1, slot = 'counts'))
cur.plotdata.umap$hasmed12 <- 0
cur.plotdata.umap$hasmed12[rownames(cur.plotdata.umap) %in% cellswithmed12] <- 1
all_int_final.noHBBcells.mes@meta.data$hasvariant <- "no"
all_int_final.noHBBcells.mes@meta.data$hasvariant[which(all_int_final.noHBBcells.mes@meta.data$med12_variant != "UNKNOWN")] <- "yes"
withVariant <- subset(all_int_final.noHBBcells.mes, subset=hasvariant == "yes")
VlnPlot(withVariant, features = "nCount_RNA", group.by = "condition")

cur.plotdata.umap <- cur.plotdata.umap[order(cur.plotdata.umap$med12_variant),]
variant_colors = c("Fibroid-55_12882_GEXL_analysis"=COLORSCHEME2[1], "Fibroid55-12906-2_GEXL"=COLORSCHEME2[2],
                   #"Fibroid12843"=COLORSCHEME2[3],
                   "Fibroid-55_12640_GEXL_analysis"=COLORSCHEME2[4],
                   "061919Fibroid_GEXL_fastqs"=COLORSCHEME2[5], "UNKNOWN"=LIGHTGRAY)
p <- ggplot(data=subset(cur.plotdata.umap, cur.plotdata.umap$med12_variant == "UNKNOWN"), 
            aes_string(x="UMAP_1", y="UMAP_2", color = "med12_variant")) +
  geom_point(size=0.3) + 
  geom_point(data = subset(cur.plotdata.umap, cur.plotdata.umap$hasmed12 == 1), 
             aes(x=UMAP_1, y=UMAP_2), size=0.5, color = "black", shape = 1, stroke=1,) + 
  geom_point(data = subset(cur.plotdata.umap, cur.plotdata.umap$med12_variant == "Fibroid-55_12882_GEXL_analysis"), 
             aes(x=UMAP_1,y=UMAP_2, color=med12_variant), size=2) + 
  geom_point(data = subset(cur.plotdata.umap, cur.plotdata.umap$med12_variant == "Fibroid55-12906-2_GEXL"), 
             aes(x=UMAP_1,y=UMAP_2, color=med12_variant), size=2) + 
  geom_point(data = subset(cur.plotdata.umap, cur.plotdata.umap$med12_variant == "Fibroid12843"), 
             aes(x=UMAP_1,y=UMAP_2, color=med12_variant), size=2) + 
  geom_point(data = subset(cur.plotdata.umap, cur.plotdata.umap$med12_variant == "Fibroid-55_12640_GEXL_analysis"), 
             aes(x=UMAP_1,y=UMAP_2, color=med12_variant), size=2) + 
  geom_point(data = subset(cur.plotdata.umap, cur.plotdata.umap$med12_variant == "061919Fibroid_GEXL_fastqs"), 
             aes(x=UMAP_1,y=UMAP_2, color=med12_variant), size=2) + 
  # geom_point(data = cur.plotdata.umap[which(cur.plotdata.umap$hasmed12 == 1 & cur.plotdata.umap$med12_variant != "UNKNOWN"),], 
  #            aes(x=UMAP_1,y=UMAP_2), size=0.5, color = "red", shape = 1, stroke = 1) + 
  scale_color_manual(name = "med12 variant source", values = variant_colors) + 
  # scale_fill_identity(guide = "legend" name = "med12 expression", labels = c("med12 counts > 1", "med12 counts > 1 and variant detected"="#FF0000")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) 
ggsave("new_freemuxlet_outs_umap_withAllMed12_fixed_v2.pdf",  width=10, height=7)


# clonality DE ---- 

all_int_final.noHBBcells.mes@meta.data$hasmed12 <- "no"
all_int_final.noHBBcells.mes@meta.data$hasmed12[which(
  rownames(all_int_final.noHBBcells.mes@meta.data) %in% cellswithmed12)] <- "yes"
all_int_final.noHBBcells.mes@meta.data$hasmed12[which(
  all_int_final.noHBBcells.mes@meta.data$hasvariant == "yes")] <- "variant"
Idents(all_int_final.noHBBcells.mes) <- all_int_final.noHBBcells.mes@meta.data$hasmed12
variant.markers <- FindMarkers(all_int_final.noHBBcells.mes, 
                               ident.1 = "variant", ident.2 = "yes")


loadTableAsDF <- function(tableFile) {
  ourTable <- read.table(tableFile, sep = "\t")
  ourTable <- as.matrix(ourTable)
  ourTable[1,1] <- "placeholder"
  colnames(ourTable) <- ourTable[1,]
  rownames(ourTable) <- ourTable[,1]
  ourTable <- ourTable[-1,]
  ourTable <- ourTable[,-1]
  ourTable <- as.data.frame(ourTable)
  return(ourTable)
}

#### Num Cells ####
getCellNumbers <- function(seurat_object, celltype_col = "seurat_clusters", sample_col = "orig.ident") {
  sampleNames <- unique(seurat_object@meta.data[,sample_col])
  meta <- seurat_object@meta.data
  celltypes <- unique(meta[,celltype_col])
  if (!is.na(as.numeric(celltypes[1]))) {
    celltypes <- sort(as.numeric(as.character(celltypes)))
  } else print("chosen a non-default column, non-numeric")
  
  clh.incelltype <- c()
  cur.vec <- c()
  for (c in celltypes) {
    cur.vec <- c()
    for (i in 1:length(sampleNames)) {
      if (i==1) {
        cur.vec <- c(length(meta[which(meta[,sample_col] == sampleNames[i] & meta[,celltype_col] == c), celltype_col]))
      }
      cur.vec <- c(cur.vec,  length(meta[which(meta[,sample_col] == sampleNames[i] & meta[,celltype_col] == c), celltype_col]))
    }
    clh.incelltype <- rbind(clh.incelltype, cur.vec)
  }
  clh.incelltype <- clh.incelltype[,-1]
  colnames(clh.incelltype) <- sampleNames
  rownames(clh.incelltype) <- celltypes
  clh.incelltype.melt <- melt(clh.incelltype)
  colnames(clh.incelltype.melt) <- c("celltype", "sample", "number_cells")
  clh.incelltype.melt$celltype <- as.factor(clh.incelltype.melt$celltype)
  write.table(clh.incelltype, paste0("cellNumbers_", seurat_object@project.name, ".txt"), quote = F, col.names = NA, sep = "\t")
  
  ## barplot1 
  ggplot(data=clh.incelltype.melt,aes(x=celltype, y=number_cells, group=sample, color=sample)) + geom_bar(stat="identity", position = "dodge", width=0.7) +
    labs(y="Number of cells", x="celltype") + 
    guides(fill=guide_legend(title="Sample")) +
    scale_y_continuous(expand = c(0,0)) + #removes the annoying space between x axis label and the bottom of plot
    theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1), panel.background = element_rect(fill="#E5E5E5"), 
          panel.grid.major = element_line(size=0.5,color="white"), panel.grid.minor = element_line(size=0.2,color="white"),
          panel.border = element_blank() )
  ggsave("cells_per_sample_per_celltype.pdf")
  
  ### barplot2
  ggplot(data=clh.incelltype.melt,aes(x=sample, y=number_cells, group=celltype, color=celltype)) + geom_bar(stat="identity", position = "dodge", width=0.7) +
    labs(y="Number of cells", x="celltype") + 
    guides(fill=guide_legend(title="Sample")) +
    scale_y_continuous(expand = c(0,0), limits=c(0,3000)) + #removes the annoying space between x axis label and the bottom of plot
    theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1), panel.background = element_rect(fill="#E5E5E5"), 
          panel.grid.major = element_line(size=0.5,color="white"), panel.grid.minor = element_line(size=0.2,color="white"),
          panel.border = element_blank() )
  ggsave("cells_per_celltype_per_sample.pdf")
  
  return(clh.incelltype)
  print("done")
}


all_markers <- FindAllMarkers(object = all_int_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, paste0("markers_all_int.txt"), quote = F, col.names = NA, sep = "\t")

umap_data <- umap_withPies(all_int_final.noHBB, pie_cond = "condition", pie_loc = "celltype", normalize = T) #test
umap_withPies(seurat_object, "condition", normalize=T)

all_int_final.noHBBcells.mes@project.name <- "all_int_final.noHBBcells.mes"
all_int_final.noHBBcells.mes.endo@project.name <- "all_int_final.noHBBcells.mes.endo"
all_int_final.noHBBcells.mes.smc@project.name <- "all_int_final.noHBBcells.mes.smc"
all_int_final.noHBBcells.mes.fib@project.name <- "all_int_final.noHBBcells.mes.fib"
for (i in list(all_int_final.noHBBcells.mes, all_int_final.noHBBcells.mes.endo, all_int_final.noHBBcells.mes.smc, all_int_final.noHBBcells.mes.fib)) {
  x <- UMAPPlot(i, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
          panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) +
    guides(color=guide_legend(ncol=2, override.aes = list(size=10)))
  y <- ggplot_build(x)
  y$data[[2]]$fontface <- 2
  y$data[[2]]$size <- 4
  z <- ggplot_gtable(y)
  pdf(paste0("nohbbcells/", i@project.name, ".pdf"), width=7, height=5)
  plot(z) 
  dev.off()
}


#### MANUSCRIPT FIGURES ####
#beware: below also overwrites current folder's umap_with_pies.pdf plot
umap_withMeta <- umap_withPies(all_int_final.noHBBcells.mes.smc, "condition", normalize = T)
umap_withMeta$orig.ident[grep("fibroid55_12746_lane[12]", umap_withMeta$orig.ident)] <- "fibroid55_12746"

COLORSCHEME1 = c("#55A4B0", "#5577B0", "#8E55B0")
COLORSCHEME2 = c("#3E80A0", "#5E3EA0", "#A03E80", "#3ea05e", "#a05e3e", "#A1903F", "#a23f3f" ,"#3fa27a", "#a23f68")

# Fig 1 umap by cluster ----
# regular umap labeled by cluster #
x <- UMAPPlot(all_int_final.nohbbcells, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.nohbbcells, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.6, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=2, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig1_umap_clustering.pdf", width=7, height=5)
plot(z) 
dev.off()

# Fig 1 umap split by condition  ----
umap_withMeta <- umap_withPies(all_int_final.nohbbcells, "condition", normalize = F, onlyumapwithmeta = T)

p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  labs(title = "MED12+")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(),
        axis.title = element_blank(),#element_text(hjust = 0), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        # axis.line.x.bottom = element_line(color = "black", size=0.8,,arrow = arrow(angle=35, length = unit(0.1, "inches"))),
        panel.grid = element_blank()) +
  labs(title = "MED12-")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/FINAL/fig1_splitByCondition_v2.pdf", width=4, height=12)

#umap with cells labeled by type
UMAPPlot(seurat_object, group.by = "celltype", pt.size = 0.2) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) +
  guides(color = guide_legend(override.aes = list(size = 10))) 
ggsave("manuscript/FINAL/fig1_celltype_labeled.pdf", height=5, width=7)


# Fig 1 dotplot ----
all_int_final.noHBBcells.markers <- read.table("nohbbcells/all_int_final_markers.txt",  sep = "\t", row.names = 1)
colnames(all_int_final.noHBB.markers) <- as.vector(unlist(all_int_final.noHBB.markers[1,]))
all_int_final.noHBB.markers <- all_int_final.noHBB.markers[-1,]
all_int_final.noHBB.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC) -> top10
top10s <- top10$gene

# seurat_object <- all_int_final.noHBB
bcell.feat <- c("MS4A1", "CD79A")
nk.feat <- c("CCL5", "GNLY", "PRF1", "NKG7")
myeloid.feat <- c("FCER1G", "CSF1R", "CD68", "CD14")
t.feat <- c("CD3E", "CD3D")
fibro.feat <- c("THY1", "DCN", "LUM")
smc.feat <- c("COL6A2", "ACTA2", "CNN1", "MYH11")
endo.feat <- c("VIM", "PDPN", "EGFL7", "CCL21", "ACKR1", "PECAM1", "VWF")
feats <- c(endo.feat, smc.feat, fibro.feat, t.feat, nk.feat, myeloid.feat, bcell.feat)
seurat_object <- SetIdent(seurat_object, value="celltype")
Idents(seurat_object) <- factor(Idents(seurat_object), levels=rev(c("B", "Myeloid", "NK", "T", "Fibroblast", "SMC", "Endothelial")))
# seurat_object <- BuildClusterTree(seurat_object, assay = "integrated", features = feats, reorder = T)
DotPlot(seurat_object, assay = "integrated", features=feats, cols = "RdYlBu", dot.scale = 8) + #scale_colour_gradient(low = LIGHTGRAY, high = "#009933") + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=12), axis.title = element_text(face="bold"))
ggsave("manuscript/FINAL/fig1_dotplot.pdf", width=12, height=5)

# Fig 1 supp feature plots ----
nice_features <- c("MYH11", "CNN1", "ACTA2", "CD3D", "CD3E", "NKG7", "PRF1", "GNLY", "CCL5", "LUM", "DCN", "COL6A2", "CD79A", "MS4A1", "HBB", "HBA2", 
                   "VWF", "PECAM1", "ACKR1", "CD14", "CD68", "CSF1R", "FCER1G", "CCL21", "EGFL7", "PDPN", "VIM", "THY1")

plot.list <- list()
LIGHTGRAY = "#F0F0F0"
for (f in nice_features) {
  plot.list[[f]] <- FeaturePlot(all_int_final.nohbbcells, cols = c("#F0F0F0", "#990099"), features = f) +
    theme(title = element_text(face="bold", size=12),
          axis.line = element_blank(),
          axis.text = element_blank(), 
          axis.title = element_blank(),
          axis.ticks = element_blank(), 
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank()) +
    labs(title = f)
}
plot_grid(plotlist = plot.list, nrow = 7, ncol = 4)
ggsave("manuscript/FINAL/featureplot_all.pdf", width=20, height=32)

# mes #
nice_features <- c("MYH11", "CNN1", "ACTA2", "CD3D", "CD3E", "NKG7", "PRF1", "GNLY", "CCL5", "LUM", "DCN", "COL6A2", "CD79A", "MS4A1", "HBB", "HBA2", 
                   "VWF", "PECAM1", "ACKR1", "CD14", "CD68", "CSF1R", "FCER1G", "CCL21", "EGFL7", "PDPN", "VIM", "THY1")
FeaturePlot(all_int_final.noHBB.mes.reclust, features = nice_features)
ggsave("manuscript/fig1_supp/featureplot_mes.pdf", width=20, height=32)


# Fig 1 supp. sample umaps ---- 
#QC umaps
colorscheme_celltype <- c("B"="#E9FD72", "Endothelial"="#BE9B33", "Fibroblast"="#6BB134", "Myeloid"="#53BD97", 
                          "NK"="#4EB4E6", "SMC"="#A18DF8", "T"="#EA6DD2")
umap_withMeta <- umap_withPies(all_int_final.nohbbcells, "condition", normalize = F, onlyumapwithmeta = T)
plots=TRUE
if (plots) {
  p_myo1 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium11911"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 1")
  p_myo2 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium061919"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 2")
  p_myo3 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium55_11564"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 3")
  p_myo4 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium55_12640"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 4")
  p_myo5 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium55_12745"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 5")
  p_myo <- as_grob(plot_grid(p_myo1, p_myo2, p_myo3, p_myo4, p_myo5, ncol=5, nrow=1)) 
  
  p_med12pos1 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid061919"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 1")
  p_med12pos2 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12382"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 2")
  p_med12pos3 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12640"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 3")
  p_med12pos4 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "med12pos55_12843"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 4")
  p_med12pos5 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12906"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 5")
  p_med12pos <- as_grob(plot_grid(p_med12pos1, p_med12pos2, p_med12pos3, p_med12pos4, p_med12pos5, ncol=5, nrow=1))
  
  p_med12neg1 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_11564"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Negative 1")
  p_med12neg2 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12716"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Negative 2")
  p_med12neg3 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12746"),], aes(x=UMAP_1, y=UMAP_2, color=celltype)) + geom_point(size=0.3) +
    scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Negative 3")
  # p_legend <- 
  p_med12neg <- as_grob(plot_grid(p_med12neg1, p_med12neg2, p_med12neg3, "", "", ncol=5, nrow=1))
  
  plot_grid(p_myo, p_med12pos, p_med12neg, nrow=3, ncol=1)
  ggsave("manuscript/FINAL/fig1supp_c.pdf", width=18, height=12)
  
}

# Fig 1 heatmap ----
Idents(all_int_final.nohbbcells) <- all_int_final.nohbbcells@meta.data$celltype
markers.all <- FindAllMarkers(all_int_final.nohbbcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.all, "manuscript/FINAL/fig1/heatmap_all/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.all$gene)*12)*(1/80)+1
markers.all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.nohbbcells, features = top10genes, draw.lines = T, lines.width = 80) 
ggsave("manuscript/FINAL/fig1/heatmap_all/heatmap.pdf", height=20, width=15)



#ratios of celltype in each sample
celltypeNums <- getCellNumbers(all_int_final.nohbbcells, celltype_col = "celltype", write_output = F, 
                               outdir = paste0(homedir))
celltypeNums.pctOfSmple <- round(t(celltypeNums) / colSums(celltypeNums), 3)*100 
celltypeNums.pctOfSmple.melt <- melt(celltypeNums.pctOfSmple)
colnames(celltypeNums.pctOfSmple.melt) <- c( "sample", "celltype","number_of_cells")
condition_ordering <- as.vector(sapply(celltypeNums.pctOfSmple.melt$sample, function(x) {
  return(unique(all_int_final.nohbbcells@meta.data$condition[which(all_int_final.nohbbcells@meta.data$orig.ident == x)]))}))
covariates <- data.frame("condition"=condition_ordering) #when celltypeNums uses orig.ident
celltypeNums.wconds <- cbind(celltypeNums.pctOfSmple.melt, covariates)
ggplot(data = celltypeNums.wconds, aes(x=sample, y=number_of_cells, fill=celltype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~condition, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))
ggsave("celltype_percent_of_sample_barplot.pdf")
write.table(celltypeNums.wconds, "celltypenums_used_for_barplot.txt", 
            sep = "\t", quote = F, col.names = NA)


# Fig 2 smc umap ----
all_int_final.noHBBcells.mes.smc <- readRDS("saves/all_int_final.noHBBcells_mes_smc_clust1and2merged.rds")
#change the cluster numbers so no number-skipping
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters <- as.character(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters)
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters == 3)] <- 2
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters == 4)] <- 3
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters == 5)] <- 4
all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters <- as.factor(all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters)
x <- UMAPPlot(all_int_final.noHBBcells.mes.smc, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.noHBBcells.mes.smc, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.6, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig2_smc_umap_clustering.pdf", width=6.5, height=5)
plot(z)
dev.off()

# Fig 2 umap split by condition  ----
umap_withMeta <- umap_withPies(all_int_final.noHBBcells.mes.smc, "condition", normalize = F, onlyumapwithmeta = T)
p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  labs(title = "MED12+")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(),
        axis.title = element_blank(),#element_text(hjust = 0), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        # axis.line.x.bottom = element_line(color = "black", size=0.8,,arrow = arrow(angle=35, length = unit(0.1, "inches"))),
        panel.grid = element_blank()) +
  labs(title = "MED12-")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/FINAL/fig2_splitByCondition_v2.pdf", width=4, height=12)

# Fig 2 smc heatmap ----
# cur.markers <- FindAllMarkers(all_int_final.noHBBcells.mes.smc,  logfc.threshold = 0.25, only.pos = TRUE) #assay = RNA
Idents(all_int_final.noHBBcells.mes.smc) <- all_int_final.noHBBcells.mes.smc@meta.data$seurat_clusters
markers.smc <- FindAllMarkers(all_int_final.noHBBcells.mes.smc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.smc, "manuscript/FINAL/fig2/heatmap/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.smc$gene)*12)*(1/80)+1
markers.smc %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.noHBBcells.mes.smc, features = top10genes, draw.lines = T, lines.width = 80) 
ggsave("manuscript/FINAL/fig2/heatmap/heatmap.pdf", height=20, width=15)

# Fig 2 smc barplot ----
samples = c("fibroid061919", "fibroid55_12382", "fibroid55_12640", "med12pos55_12843", "fibroid55_12906", "fibroid55_11564",
            "fibroid55_12716", "fibroid55_12746", "myometrium061919", "myometrium11911", "myometrium55_11564", 
            "myometrium55_12640", "myometrium55_12745")
barPlotsWmedian(seurat_object = all_int_final.noHBBcells.mes.smc, stats = "all",
                name = "SMC", ordering = samples,
                celltype_col = "seurat_clusters", type = "mean", outdir = "manuscript/FINAL")
#check --
smc <- getCellNumbers(all_int_final.noHBBcells.mes.smc)
smc.p <- round(t(smc) / colSums(smc), 3)*100 
condition_ordering <- as.vector(sapply(rownames(smc.p), function(x) {
  return(unique(all_int_final.noHBBcells.mes.smc@meta.data$condition[which(all_int_final.noHBBcells.mes.smc@meta.data$orig.ident == x)]))}))
covariates <- data.frame("condition"=condition_ordering) #when celltypeNums uses orig.ident
smc.p.wconds <- cbind(smc.p, covariates)
smc.p.wconds <- smc.p.wconds[order(smc.p.wconds$condition),]
smc.p.wconds$means <- 999
smc.p.wconds$means[which(smc.p.wconds$condition == "med12neg")] <- 
  paste0(colMeans(smc.p.wconds[which(smc.p.wconds$condition == "med12neg"),grep("^[0-9]$", colnames(smc.p.wconds))]), collapse = ", ")
smc.p.wconds$means[which(smc.p.wconds$condition == "med12pos")] <- 
  paste0(colMeans(smc.p.wconds[which(smc.p.wconds$conditiontion == "med12pos"),grep("^[0-9]$", colnames(smc.p.wconds))]), collapse = ", ")
smc.p.wconds$means[which(smc.p.wconds$condition == "myometrium")] <- 
  paste0(colMeans(smc.p.wconds[which(smc.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(smc.p.wconds))]), collapse = ", ")
smc.p.wconds$stderr[which(smc.p.wconds$condition == "myometrium")] <- 
  paste0(colMeans(smc.p.wconds[which(smc.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(smc.p.wconds))]), collapse = ", ")
smc.p.wconds$stderr <- 999
smc.p.wconds$stderr[which(smc.p.wconds$condition == "med12neg")] <- paste0(apply(smc.p.wconds[which(smc.p.wconds$condition == "med12neg"),grep("^[0-9]$", colnames(smc.p.wconds))], 2,
                                                                                 function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")
smc.p.wconds$stderr[which(smc.p.wconds$condition == "med12pos")] <- paste0(apply(smc.p.wconds[which(smc.p.wconds$condition == "med12pos"),grep("^[0-9]$", colnames(smc.p.wconds))], 2,
                                                                                 function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")
smc.p.wconds$stderr[which(smc.p.wconds$condition == "myometrium")] <- paste0(apply(smc.p.wconds[which(smc.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(smc.p.wconds))], 2,
                                                                                   function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")

# Ratio of SMC to Fib ----
smc.persample <- table(all_int_final.nohbbcells@meta.data$orig.ident[which(all_int_final.nohbbcells@meta.data$celltype == "SMC")])
fibroblast.persample <- table(all_int_final.nohbbcells@meta.data$orig.ident[which(all_int_final.nohbbcells@meta.data$celltype == "Fibroblast")])
ratio.smc2fib <- smc.persample / fibroblast.persample
ratio.tab <- do.call(cbind, list(smc.persample, fibroblast.persample, ratio.smc2fib))
colnames(ratio.tab) <- c("smc", "fibroblast", "ratio")
write.table(ratio.tab, "smc2fib_ratio.txt", sep = "\t", quote = F, col.names = NA)


# Fig 3 fibro umap ----
all_int_final.noHBBcells.mes.fib <- readRDS("saves/all_int_final.nohbbcells.mes.fib.onlyclusters6and11.rds")
#change the cluster numbers so no number-skipping
all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters <- as.character(all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters)
all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters == 6)] <- 0
all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters == 11)] <- 1
all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters <- as.factor(all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters)
x <- UMAPPlot(all_int_final.noHBBcells.mes.fib, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.noHBBcells.mes.fib, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.75, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig3_fib_umap_clustering.pdf", width=6.5, height=5)
plot(z)
dev.off()

# Fig 3 fib umap split by condition  ----
umap_withMeta <- umap_withPies(all_int_final.noHBBcells.mes.fib, "condition", normalize = F, onlyumapwithmeta = T)
p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  labs(title = "MED12+")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(),
        axis.title = element_blank(),#element_text(hjust = 0), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        # axis.line.x.bottom = element_line(color = "black", size=0.8,,arrow = arrow(angle=35, length = unit(0.1, "inches"))),
        panel.grid = element_blank()) +
  labs(title = "MED12-")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/FINAL/fig3_splitByCondition_v2.pdf", width=4, height=12)

# Fig 3 fib heatmap ----
# cur.markers <- FindAllMarkers(all_int_final.noHBBcells.mes.smc,  logfc.threshold = 0.25, only.pos = TRUE) #assay = RNA
Idents(all_int_final.noHBBcells.mes.fib) <-all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters
markers.fib <- FindAllMarkers(all_int_final.noHBBcells.mes.fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.fib, "manuscript/FINAL/fig3/heatmap/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.fib$gene)*12)*(1/80)+1
markers.fib %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.noHBBcells.mes.fib, features = top10genes, draw.lines = T, lines.width = 20) 
ggsave("manuscript/FINAL/fig3/heatmap/heatmap.pdf", height=20, width=15)

# Fig 3 fib barplot ----
samples = c("fibroid061919", "fibroid55_12382", "fibroid55_12640", "med12pos55_12843", "fibroid55_12906", "fibroid55_11564",
            "fibroid55_12716", "fibroid55_12746", "myometrium061919", "myometrium11911", "myometrium55_11564", 
            "myometrium55_12640", "myometrium55_12745")
Idents(all_int_final.noHBBcells.mes.fib) <- all_int_final.noHBBcells.mes.fib@meta.data$seurat_clusters
barPlotsWmedian(seurat_object = all_int_final.noHBBcells.mes.fib, stats = "all",
                name = "fibro", ordering = samples,
                celltype_col = "seurat_clusters", type = "mean", outdir = "manuscript/FINAL")
#check --
fib <- getCellNumbers(all_int_final.noHBBcells.mes.fib)
fib.p <- round(t(fib) / colSums(fib), 3)*100 
condition_ordering <- as.vector(sapply(rownames(fib.p), function(x) {
  return(unique(all_int_final.noHBBcells.mes.fib@meta.data$condition[which(all_int_final.noHBBcells.mes.fib@meta.data$orig.ident == x)]))}))
covariates <- data.frame("condition"=condition_ordering) #when celltypeNums uses orig.ident
fib.p.wconds <- cbind(fib.p, covariates)
fib.p.wconds <- fib.p.wconds[order(fib.p.wconds$condition),]
fib.p.wconds$means <- 999
fib.p.wconds$means[which(fib.p.wconds$condition == "med12neg")] <- 
  paste0(colMeans(fib.p.wconds[which(fib.p.wconds$condition == "med12neg"),grep("^[0-9]$", colnames(fib.p.wconds))]), collapse = ", ")
fib.p.wconds$means[which(fib.p.wconds$condition == "med12pos")] <- 
  paste0(colMeans(fib.p.wconds[which(fib.p.wconds$condition == "med12pos"),grep("^[0-9]$", colnames(fib.p.wconds))]), collapse = ", ")
fib.p.wconds$means[which(fib.p.wconds$condition == "myometrium")] <- 
  paste0(colMeans(fib.p.wconds[which(fib.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(fib.p.wconds))]), collapse = ", ")
fib.p.wconds$stderr[which(fib.p.wconds$condition == "myometrium")] <- 
  paste0(colMeans(fib.p.wconds[which(fib.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(fib.p.wconds))]), collapse = ", ")
fib.p.wconds$stderr <- 999
fib.p.wconds$stderr[which(fib.p.wconds$condition == "med12neg")] <- paste0(apply(fib.p.wconds[which(fib.p.wconds$condition == "med12neg"),grep("^[0-9]$", colnames(fib.p.wconds))], 2,
                                                                                 function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")
fib.p.wconds$stderr[which(fib.p.wconds$condition == "med12pos")] <- paste0(apply(fib.p.wconds[which(fib.p.wconds$condition == "med12pos"),grep("^[0-9]$", colnames(fib.p.wconds))], 2,
                                                                                 function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")
fib.p.wconds$stderr[which(fib.p.wconds$condition == "myometrium")] <- paste0(apply(fib.p.wconds[which(fib.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(fib.p.wconds))], 2,
                                                                                   function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")

# Fig 4 endo umap ----
all_int_final.noHBBcells.mes.endo <- readRDS("saves/all_int_final.noHBBcells_noHBBgenes_noNK.mes.endo.rds")
x <- UMAPPlot(all_int_final.noHBBcells.mes.endo, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.noHBBcells.mes.endo, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.75, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig4_endo_umap_clustering.pdf", width=6.5, height=5)
plot(z)
dev.off()

# Fig 4 endo umap split by condition  ----
umap_withMeta <- umap_withPies(all_int_final.noHBBcells.mes.endo, "condition", normalize = F, onlyumapwithmeta = T)
p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  labs(title = "MED12+")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(),
        axis.title = element_blank(),#element_text(hjust = 0), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        # axis.line.x.bottom = element_line(color = "black", size=0.8,,arrow = arrow(angle=35, length = unit(0.1, "inches"))),
        panel.grid = element_blank()) +
  labs(title = "MED12-")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/FINAL/fig4_splitByCondition_v2.pdf", width=4, height=12)

# Fig 4 sample umap ----
x <- UMAPPlot(all_int_final.noHBBcells.mes.endo, pt.size = 0.3, split.by="orig.ident", group.by="seurat_clusters") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig4_endo_umap_bysample.pdf", width=20, height=4)
plot(x)
dev.off()

umap_withMeta <- umap_withPies(all_int_final.noHBBcells.mes.endo, "condition", normalize = F, onlyumapwithmeta = T)
plots=TRUE
if (plots) {
  p_myo1 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium11911"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 1")
  p_myo2 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium061919"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 2")
  p_myo3 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium55_11564"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 3")
  p_myo4 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium55_12640"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 4")
  p_myo5 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "myometrium55_12745"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Myometrium 5")
  p_myo <- as_grob(plot_grid(p_myo1, p_myo2, p_myo3, p_myo4, p_myo5, ncol=5, nrow=1)) 
  
  p_med12pos1 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid061919"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 1")
  p_med12pos2 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12382"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 2")
  p_med12pos3 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12640"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 3")
  p_med12pos4 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "med12pos55_12843"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 4")
  p_med12pos5 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12906"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Positive 5")
  p_med12pos <- as_grob(plot_grid(p_med12pos1, p_med12pos2, p_med12pos3, p_med12pos4, p_med12pos5, ncol=5, nrow=1))
  
  p_med12neg1 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_11564"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Negative 1")
  p_med12neg2 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12716"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Negative 2")
  p_med12neg3 <- ggplot(data=umap_withMeta[which(umap_withMeta$orig.ident == "fibroid55_12746"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters)) + 
    geom_point(size=0.3) +
    # scale_color_manual(values = colorscheme_celltype) +
    theme(axis.title = element_blank(), axis.text = element_blank(), title = element_text(face="bold", size=12), legend.position = "none",
          axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
    labs(title = "Med12 Negative 3")
  # p_legend <- 
  p_med12neg <- as_grob(plot_grid(p_med12neg1, p_med12neg2, p_med12neg3, "", "", ncol=5, nrow=1))
  
  plot_grid(p_myo, p_med12pos, p_med12neg, nrow=3, ncol=1)
  ggsave("manuscript/FINAL/fig4_endo_bysample.pdf", width=18, height=12)
  
}

# Fig 4 endo heatmap ----
# cur.markers <- FindAllMarkers(all_int_final.noHBBcells.mes.endo,  logfc.threshold = 0.25, only.pos = TRUE) #assay = RNA
Idents(all_int_final.noHBBcells.mes.endo) <- all_int_final.noHBBcells.mes.endo@meta.data$seurat_clusters
markers.endo <- FindAllMarkers(all_int_final.noHBBcells.mes.endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.endo, "manuscript/FINAL/fig4/heatmap/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.endo$gene)*12)*(1/80)+1
markers.endo %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.noHBBcells.mes.endo, features = top10genes, draw.lines = T, lines.width = 80) 
ggsave("manuscript/FINAL/fig4/heatmap/heatmap.pdf", height=20, width=15)

# Fig 4 endo barplot ----
samples = c("fibroid061919", "fibroid55_12382", "fibroid55_12640", "med12pos55_12843", "fibroid55_12906", "fibroid55_11564",
            "fibroid55_12716", "fibroid55_12746", "myometrium061919", "myometrium11911", "myometrium55_11564", 
            "myometrium55_12640", "myometrium55_12745")
barPlotsWmedian(seurat_object = all_int_final.noHBBcells.mes.endo, stats = "all",
                name = "endo", ordering = samples, celltype_col = "seurat_clusters",
                type = "mean", outdir = "manuscript/FINAL/fig4")
#check --
endo <- getCellNumbers(all_int_final.noHBBcells.mes.endo)
endo.p <- round(t(endo) / colSums(endo), 3)*100 
condition_ordering <- as.vector(sapply(rownames(endo.p), function(x) {
  return(unique(all_int_final.noHBBcells.mes.endo@meta.data$condition[which(all_int_final.noHBBcells.mes.endo@meta.data$orig.ident == x)]))}))
covariates <- data.frame("condition"=condition_ordering) #when celltypeNums uses orig.ident
endo.p.wconds <- cbind(endo.p, covariates)
endo.p.wconds <- endo.p.wconds[order(endo.p.wconds$condition),]
endo.p.wconds$means <- 999
endo.p.wconds$means[which(endo.p.wconds$condition == "med12neg")] <- 
  paste0(colMeans(endo.p.wconds[which(endo.p.wconds$condition == "med12neg"),grep("^[0-9]$", colnames(endo.p.wconds))]), collapse = ", ")
endo.p.wconds$means[which(endo.p.wconds$condition == "med12pos")] <- 
  paste0(colMeans(endo.p.wconds[which(endo.p.wconds$condition == "med12pos"),grep("^[0-9]$", colnames(endo.p.wconds))]), collapse = ", ")
endo.p.wconds$means[which(endo.p.wconds$condition == "myometrium")] <- 
  paste0(colMeans(endo.p.wconds[which(endo.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(endo.p.wconds))]), collapse = ", ")
endo.p.wconds$stderr[which(endo.p.wconds$condition == "myometrium")] <- 
  paste0(colMeans(endo.p.wconds[which(endo.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(endo.p.wconds))]), collapse = ", ")
endo.p.wconds$stderr <- 999
endo.p.wconds$stderr[which(endo.p.wconds$condition == "med12neg")] <- paste0(apply(endo.p.wconds[which(endo.p.wconds$condition == "med12neg"),grep("^[0-9]$", colnames(endo.p.wconds))], 2,
                                                                                   function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")
endo.p.wconds$stderr[which(endo.p.wconds$condition == "med12pos")] <- paste0(apply(endo.p.wconds[which(endo.p.wconds$condition == "med12pos"),grep("^[0-9]$", colnames(endo.p.wconds))], 2,
                                                                                   function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")
endo.p.wconds$stderr[which(endo.p.wconds$condition == "myometrium")] <- paste0(apply(endo.p.wconds[which(endo.p.wconds$condition == "myometrium"),grep("^[0-9]$", colnames(endo.p.wconds))], 2,
                                                                                     function(x) return(sd(as.numeric(x)) / sqrt(length(x)))), collapse = ", ")

# Fig 5 myo smc only ----
all_int_final.noHBBcells.mes.smc.myo <- subset(all_int_final.noHBBcells.mes.smc, subset=condition == "myometrium")
all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper <- all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident
all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident == "myometrium061919")] <- "Myometrium 1"
all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident == "myometrium11911")] <- "Myometrium 2"
all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident == "myometrium55_11564")] <- "Myometrium 3"
all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident == "myometrium55_12640")] <- "Myometrium 4"
all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident == "myometrium55_12745")] <- "Myometrium 5"
samplerename_reference <- all_int_final.noHBBcells.mes.smc.myo@meta.data[!duplicated(all_int_final.noHBBcells.mes.smc.myo@meta.data$orig.ident4paper),]
write.table(samplerename_reference,"samplerename_reference.txt", quote = F, sep = "\t", col.names = NA)

# Fig 5 myo smc umap by cluster ----
x <- UMAPPlot(all_int_final.noHBBcells.mes.smc.myo, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.noHBBcells.mes.smc.myo, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.75, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig5_smc_myo_umap_clustering.pdf", width=6.5, height=5)
plot(z)
dev.off()

# Fig 5 myo smc umap by sample ----
x <- UMAPPlot(all_int_final.noHBBcells.mes.smc.myo, pt.size = 0.3, group.by="orig.ident4paper") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig5_smc_myo_umap_bysample.pdf", width=7.5, height=5)
plot(x)
dev.off()

# Fig 5 myo smc heatmap ----
markers.smc.myo <- FindAllMarkers(all_int_final.noHBBcells.mes.smc.myo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.smc.myo, "manuscript/FINAL/fig5/heatmap/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.smc.myo$gene)*12)*(1/80)+1
markers.smc.myo %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
Idents(all_int_final.noHBBcells.mes.smc.myo) <- all_int_final.noHBBcells.mes.smc.myo@meta.data$seurat_clusters
DoHeatmap(all_int_final.noHBBcells.mes.smc.myo, features = top10genes, draw.lines = T, lines.width = 30) 
ggsave("manuscript/FINAL/fig5/heatmap_smc/heatmap.pdf", height=20, width=15)


# Fig 5 myo fib only ----
all_int_final.noHBBcells.mes.fib.myo <- subset(all_int_final.noHBBcells.mes.fib, subset=condition == "myometrium")
all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident4paper <- all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident
all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident == "myometrium061919")] <- "Myometrium 1"
all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident == "myometrium11911")] <- "Myometrium 2"
all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident == "myometrium55_11564")] <- "Myometrium 3"
all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident == "myometrium55_12640")] <- "Myometrium 4"
all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident4paper[
  which(all_int_final.noHBBcells.mes.fib.myo@meta.data$orig.ident == "myometrium55_12745")] <- "Myometrium 5"
#sample reference sheet already in folder
all_int_final.noHBBcells.mes.fib.myo@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.fib.myo@meta.data$seurat_clusters == 6)] <- 0
all_int_final.noHBBcells.mes.fib.myo@meta.data$seurat_clusters[which(all_int_final.noHBBcells.mes.fib.myo@meta.data$seurat_clusters == 11)] <- 1

# Fig 5 myo fib umap by cluster ----
x <- UMAPPlot(all_int_final.noHBBcells.mes.fib.myo, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.noHBBcells.mes.fib.myo, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.75, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig5_fib_myo_umap_clustering.pdf", width=6.5, height=5)
plot(z)
dev.off()

# Fig 5 myo fib umap by sample ----
# umap by cluster
x <- UMAPPlot(all_int_final.noHBBcells.mes.fib.myo, pt.size = 0.3, group.by="orig.ident4paper") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig5_smc_myo_umap_bysample.pdf", width=7.5, height=5)
plot(x)
dev.off()

# Fig 5 myo fib heatmap ----
Idents(all_int_final.noHBBcells.mes.fib.myo) <- all_int_final.noHBBcells.mes.fib.myo@meta.data$seurat_clusters
markers.fib.myo <- FindAllMarkers(all_int_final.noHBBcells.mes.fib.myo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.fib.myo, "manuscript/FINAL/fig5/heatmap_fib/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.fib.myo$gene)*12)*(1/80)+1
markers.fib.myo %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.noHBBcells.mes.fib.myo, features = top10genes, draw.lines = T, lines.width = 10) 
ggsave("manuscript/FINAL/fig5/heatmap_fib/heatmap.pdf", height=20, width=15)


# Fig 6 Immune umap by clusters----
all_int_final.imm <- readRDS("saves/all_int_final_noHBB_imm_reclustered_v2.rds")
all_int_final.imm@meta.data$condition[grep("12746", all_int_final.imm@meta.data$orig.ident)] <- "med12neg"

x <- UMAPPlot(all_int_final.imm, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
labels <- y$data[[2]]
z <- UMAPPlot(all_int_final.imm, pt.size = 0.3, group.by="seurat_clusters") +
  # geom_label(data=labels, aes(x=x, y=y, label = label), fontface = "bold", fill = LIGHTGRAY, alpha = 0.7) +
  geom_point(data=labels, aes(x=x, y=y), color=LIGHTGRAY, alpha=0.75, size=6) +
  geom_text(data= labels, aes(x=x,y=y,label=label, fontface = "bold")) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig6_immune_umap_clustering.pdf", width=6.5, height=5)
plot(z)
dev.off()

# Fig 6 immune celltypes  ----
x <- UMAPPlot(all_int_final.imm, pt.size = 0.3, group.by="celltype") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
y <- ggplot_build(x)
y$data[[1]]$colour[which(y$data[[1]]$group == 9)] <- "#C96132"
y$data[[1]]$colour[which(y$data[[1]]$group == 4)] <- "#F57ADD"
colors <- unique(y$data[[1]]$colour)
z <- UMAPPlot(all_int_final.imm, pt.size = 0.3, group.by="celltype", cols=colors) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18))
pdf("manuscript/FINAL/fig6_immune_umap_celltypes.pdf", width=7.5, height=5)
plot(z)
dev.off()

# Fig 6 immune umap split by condition  ----
umap_withMeta <- umap_withPies(all_int_final.imm, "condition", normalize = F, onlyumapwithmeta = T)

p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  labs(title = "MED12+")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(title = element_text(face="bold", size=14),
        axis.text = element_blank(),
        axis.title = element_blank(),#element_text(hjust = 0), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        # axis.line.x.bottom = element_line(color = "black", size=0.8,,arrow = arrow(angle=35, length = unit(0.1, "inches"))),
        panel.grid = element_blank()) +
  labs(title = "MED12-")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/FINAL/fig6_splitByCondition_v2.pdf", width=4, height=12)

# Fig 6 imm heatmap ----
Idents(all_int_final.imm) <- all_int_final.imm@meta.data$celltype
markers.imm <- FindAllMarkers(all_int_final.imm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.imm, "manuscript/FINAL/fig6/heatmap/heatmap_markers.txt", quote = F, sep = "\t", col.names = NA)
height80dpi <- (length(markers.imm$gene)*12)*(1/80)+1
markers.imm %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(all_int_final.imm, features = top10genes, draw.lines = T, lines.width = 10) 
ggsave("manuscript/FINAL/fig6/heatmap/heatmap.pdf", height=20, width=15)


# Fig 6 imm barplot ----
samples = c("fibroid061919", "fibroid55_12382", "fibroid55_12640", "med12pos55_12843", "fibroid55_12906", "fibroid55_11564",
            "fibroid55_12716", "fibroid55_12746", "myometrium061919", "myometrium11911", "myometrium55_11564", 
            "myometrium55_12640", "myometrium55_12745")
barPlotsWmedian(seurat_object = all_int_final.imm, stats = "all",
                name = "imm", ordering = samples, celltype_col = "seurat_clusters",
                type = "mean", outdir = "manuscript/FINAL/fig6")

# QC Violin figures ----
myometrium <- subset(all_int_final.nohbbcells, subset=condition=="myometrium")
med12pos <- subset(all_int_final.nohbbcells, subset=condition=="med12pos")
med12neg <- subset(all_int_final.nohbbcells, subset=condition=="med12neg")
p1.1 <- VlnPlot(myometrium, group.by = "celltype", 
                features = c("nFeature_RNA"), pt.size = 0)
p1.2 <- VlnPlot(myometrium, group.by = "celltype", 
                features = c("nCount_RNA"), pt.size = 0)
p1.3 <- VlnPlot(myometrium, group.by = "celltype", 
                features = c("percent.mt"), pt.size = 0)
p2.1 <- VlnPlot(med12pos, group.by = "celltype", 
                features = c("nFeature_RNA"), pt.size = 0)
p2.2 <- VlnPlot(med12pos, group.by = "celltype", 
                features = c("nCount_RNA"), pt.size = 0)
p2.3 <- VlnPlot(med12pos, group.by = "celltype", 
                features = c("percent.mt"), pt.size = 0)
p3.1 <- VlnPlot(med12neg, group.by = "celltype", 
                features = c("nFeature_RNA"), pt.size = 0)
p3.2 <- VlnPlot(med12neg, group.by = "celltype", 
                features = c("nCount_RNA"), pt.size = 0)
p3.3 <- VlnPlot(med12neg, group.by = "celltype", 
                features = c("percent.mt"), pt.size = 0)
plot_grid(p1.1, p1.2, p1.3, p2.1, p2.2, p2.3, p3.1, p3.2, p3.3, nrow=3, ncol=3)
ggsave("manuscript/FINAL/violin_QC_byCondition.pdf", width=18, height=15)

# endo 2,8 clusters ----
samples = c("fibroid061919", "fibroid55_12382", "fibroid55_12640", "med12pos55_12843", "fibroid55_12906", "fibroid55_11564",
            "fibroid55_12716", "fibroid55_12746", "myometrium061919", "myometrium11911", "myometrium55_11564", 
            "myometrium55_12640", "myometrium55_12745")
barPlotsWmedian(seurat_object = all_int_final.nohbbcells, 
                name = "all", ordering = samples, celltype_col = "celltype",
                type = "mean", outdir = "presentation/")
endo2_8 <- FindMarkers(all_int_final.noHBBcells.mes.endo, min.pct = 0.1, logfc.threshold = 0.1, ident.1 = 2, ident.2 = 8)
endo2_8.best <- rownames(head(endo2_8[order(abs(endo2_8$avg_logFC), decreasing = T),], 12))
all_int_final.noHBBcells.mes.endo2_8 <- subset(all_int_final.noHBBcells.mes.endo, subset = seurat_clusters %in% c(2,8))
VlnPlot(all_int_final.noHBBcells.mes.endo2_8, features = endo2_8.best, cols = c("#C601C3", "#01BFC4"), pt.size = 0.05)
ggsave("endo2_8_best.pdf")
feats2_8 <- rownames(all_int_final.noHBBcells.mes.endo2_8@assays$RNA@counts)
hmgfeats <- feats2_8[grep("^HMG", feats2_8)]
VlnPlot(all_int_final.noHBBcells.mes.endo2_8, assay = "RNA", features = hmgfeats, cols = c("#C601C3", "#01BFC4"), pt.size = 0.05)
ggsave("endo2_8_hmgs.pdf", width = 20, height = 20)

f <- FetchData(all_int_final.noHBBcells.mes, vars = "MED12")
f.f <- subset(f, f$rna_MED12 >= 1)
f.f.dat <- subset(all_int_final.noHBBcells.mes, cells = rownames(f.f))
VlnPlot(f.f.dat, features = c("nCount_RNA", "percent.mt", "nFeature_RNA"), group.by = "orig.ident")
ggsave("violin_med12_cells.png", width = 12, height = 6)


## Supp ---- 
# condition split umaps for mesenchymal and immune  #
#mesenchymal
umap_withMeta <- umap_withPies(all_int_final.noHBB.mes.reclust, "condition", normalize = T)
umap_withMeta$orig.ident[grep("fibroid55_12746_lane[12]", umap_withMeta$orig.ident)] <- "fibroid55_12746"

p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Med12 positive")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Med12 Negative")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/fig3supp_mes.pdf", width=4, height=12)

#immune
umap_withMeta <- umap_withPies(all_int_final.noHBB.imm, "condition", normalize = T)
umap_withMeta$orig.ident[grep("fibroid55_12746_lane[12]", umap_withMeta$orig.ident)] <- "fibroid55_12746"

p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Med12 positive")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Med12 Negative")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/fig3supp_imm.pdf", width=4, height=12)

#endothelial
umap_withMeta <- umap_withPies(all_int_final.noHBB.mes.endo, "condition", normalize = T)
umap_withMeta$orig.ident[grep("fibroid55_12746_lane[12]", umap_withMeta$orig.ident)] <- "fibroid55_12746"
umap_withMeta$orig.ident[grep("fibroid55_12906lane2", umap_withMeta$orig.ident)] <- "fibroid55_12906"

p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "myometrium"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[1]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Myometrium")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12pos"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[2]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Med12 positive")

p3 <- ggplot(data=umap_withMeta[which(umap_withMeta$condition == "med12neg"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.3, color=COLORSCHEME2[3]) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Med12 Negative")

plot_grid(p1, p2, p3, nrow=3, ncol=1)
ggsave("manuscript/fig3/fig3supp_endo.pdf", width=4, height=12)


## B ##
#without gray
p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$celltype_major == "Mesenchymal"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.18, color="#CB975C") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Mesenchymal")

p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$celltype_major == "Immune"),], aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.18, color="#5CC7CB") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Immune")

plot_grid(p1, p2, nrow=2, ncol=1)
ggsave("manuscript/fig3a.pdf", width=3, height=6)

#with gray
umap_withMeta$celltype_major <- factor(umap_withMeta$celltype_major, levels=c("Mesenchymal", "Immune", "Contaminating"))
p1 <- ggplot(data=umap_withMeta[which(umap_withMeta$celltype_major != "Contaminating"),], aes(x=UMAP_1, y=UMAP_2, color=celltype_major)) + geom_point(size=0.18) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Mesenchymal") + scale_color_manual(values=c("#CB975C", LIGHTGRAY))
umap_withMeta$celltype_major <- factor(umap_withMeta$celltype_major, levels=c("Immune", "Mesenchymal", "Contaminating"))
p2 <- ggplot(data=umap_withMeta[which(umap_withMeta$celltype_major != "Contaminating"),], aes(x=UMAP_1, y=UMAP_2, color=celltype_major)) + geom_point(size=0.18) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  labs(title = "Immune") + scale_color_manual(values=c("#5CC7CB", LIGHTGRAY))

plot_grid(p1, p2, nrow=2, ncol=1)
ggsave("manuscript/fig3a.pdf", width=5.3, height=8)

#single plot 
#with gray
umap_withMeta$celltype_major <- factor(umap_withMeta$celltype_major, levels=c("Mesenchymal", "Immune", "Contaminating"))
ggplot(data=umap_withMeta[which(umap_withMeta$celltype_major != "Contaminating"),], aes(x=UMAP_1, y=UMAP_2, color=celltype_major)) + geom_point(size=0.18) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
  scale_color_manual(values=c( "#5CC7CB", "#CB975C"))
ggsave("manuscript/fig3a_both.pdf", width=5.3, height=5.3)


#### Fig 4 Volcanoes #### 
# load 4 lists of markers
# MUST HAVE ALREADY RUN 'EASY_DE' SECTION
# mesenchymal med12p v myo
# volcdir <- "nohbbcells/mes/smc_final/DE_inclusters"
#
if (nohbbcells) { #smc, fibroblast, endo
  # smc med12n v myo
  setwd(volcdir)
  deList_smc_med12p <- list()
  for (f in list.files("R_Analysis/DE/smc")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"))) {
      deList_smc_med12p[[paste0(f)]] <- NULL
      next
    }
    deList_smc_med12p[[paste0(f)]] <- read.table(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_smc_med12p[[paste0(f)]]) <- as.character(unlist(deList_smc_med12p[[paste0(f)]][1,]))
    deList_smc_med12p[[paste0(f)]] <- deList_smc_med12p[[paste0(f)]][-1,]
  }
  # smc med12n v myo
  setwd("R_Analysis/DE/smc") #load all of the med12pvmyo & med12nvmyo
  deList_smc_med12n <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"))) {
      deList_smc_med12n[[paste0(f)]] <- NULL
      next
    }
    deList_smc_med12n[[paste0(f)]] <- read.table(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_smc_med12n[[paste0(f)]]) <- as.character(unlist(deList_smc_med12n[[paste0(f)]][1,]))
    deList_smc_med12n[[paste0(f)]] <- deList_smc_med12n[[paste0(f)]][-1,]
  }
  
  # fibroblast med12p v myo
  setwd("R_Analysis/DE/fib") #load all of the med12pvmyo & med12nvmyo
  deList_fib_med12p <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"))) {
      deList_fib_med12p[[paste0(f)]] <- NULL
      next
    }
    deList_fib_med12p[[paste0(f)]] <- read.table(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_fib_med12p[[paste0(f)]]) <- as.character(unlist(deList_fib_med12p[[paste0(f)]][1,]))
    deList_fib_med12p[[paste0(f)]] <- deList_fib_med12p[[paste0(f)]][-1,]
  }
  
  # fibroblast med12n v myo
  setwd("R_Analysis/DE/fib") #load all of the med12pvmyo & med12nvmyo
  deList_fib_med12n <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"))) {
      deList_fib_med12n[[paste0(f)]] <- NULL
      next
    }
    deList_fib_med12n[[paste0(f)]] <- read.table(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_fib_med12n[[paste0(f)]]) <- as.character(unlist(deList_fib_med12n[[paste0(f)]][1,]))
    deList_fib_med12n[[paste0(f)]] <- deList_fib_med12n[[paste0(f)]][-1,]
  }
  
  # endothelial med12p v myo
  setwd("nohbbcells/mes/endo/DE_inclusters") #load all of the med12pvmyo & med12nvmyo
  deList_endo_med12p <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"))) {
      deList_endo_med12p[[paste0(f)]] <- NULL
      next
    }
    deList_endo_med12p[[paste0(f)]] <- read.table(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_endo_med12p[[paste0(f)]]) <- as.character(unlist(deList_endo_med12p[[paste0(f)]][1,]))
    deList_endo_med12p[[paste0(f)]] <- deList_endo_med12p[[paste0(f)]][-1,]
  }
  
  # endothelial med12n v myo
  setwd("nohbbcells/mes/endo/DE_inclusters") #load all of the med12pvmyo & med12nvmyo
  deList_endo_med12n <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"))) {
      deList_endo_med12n[[paste0(f)]] <- NULL
      next
    }
    deList_endo_med12n[[paste0(f)]] <- read.table(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_endo_med12n[[paste0(f)]]) <- as.character(unlist(deList_endo_med12n[[paste0(f)]][1,]))
    deList_endo_med12n[[paste0(f)]] <- deList_endo_med12n[[paste0(f)]][-1,]
  }
  
} else {
  
  setwd("FINAL_merge/hbb_regressed/mesenchymal/DE_inclusters") #load all of the med12pvmyo & med12nvmyo
  deList_mes_med12p <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    deList_mes_med12p[[paste0(f)]] <- read.table(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_mes_med12p[[paste0(f)]]) <- as.character(unlist(deList_mes_med12p[[paste0(f)]][1,]))
    deList_mes_med12p[[paste0(f)]] <- deList_mes_med12p[[paste0(f)]][-1,]
  }
  # mesenchymal med12n v myo
  setwd("FINAL_merge/hbb_regressed/mesenchymal/DE_inclusters") #load all of the med12pvmyo & med12nvmyo
  deList_mes_med12n <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    deList_mes_med12n[[paste0(f)]] <- read.table(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_mes_med12n[[paste0(f)]]) <- as.character(unlist(deList_mes_med12n[[paste0(f)]][1,]))
    deList_mes_med12n[[paste0(f)]] <- deList_mes_med12n[[paste0(f)]][-1,]
  }
  
  # immune med12p v myo
  setwd("R_Analysis/DE/imm") #load all of the med12pvmyo & med12nvmyo
  deList_imm_med12p <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"))) {
      deList_imm_med12p[[paste0(f)]] <- NULL
      next
    }
    deList_imm_med12p[[paste0(f)]] <- read.table(paste0(f, "/med12pos_vs_myometrium/DEmarkers_", clustnum,"_med12pos_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_imm_med12p[[paste0(f)]]) <- as.character(unlist(deList_imm_med12p[[paste0(f)]][1,]))
    deList_imm_med12p[[paste0(f)]] <- deList_imm_med12p[[paste0(f)]][-1,]
  }
  
  # immune med12n v myo
  setwd("R_Analysis/DE/imm") #load all of the med12pvmyo & med12nvmyo
  deList_imm_med12n <- list()
  for (f in list.files(".")) { #load the list of DE tables
    clustnum <- unlist(strsplit(f, "_"))[2]
    if (!file.exists(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"))) {
      deList_imm_med12n[[paste0(f)]] <- NULL
      next
    }
    deList_imm_med12n[[paste0(f)]] <- read.table(paste0(f, "/med12neg_vs_myometrium/DEmarkers_", clustnum,"_med12neg_vs_myometrium.txt"),  sep = "\t", row.names = 1)
    colnames(deList_imm_med12n[[paste0(f)]]) <- as.character(unlist(deList_imm_med12n[[paste0(f)]][1,]))
    deList_imm_med12n[[paste0(f)]] <- deList_imm_med12n[[paste0(f)]][-1,]
  }
}


fib_med12p.genes <- c("RBP1", "ALDH1A2", "NR4A1", "GAS1", "VCAM1", "PGF", "COL1A2", "TKT", "COMP", "MMP11")
fib_med12n.genes <- c("MT1M", "MT1A", "MT1G", "KRT8", "FABP4", "HYAL2")
smc_med12p.genes <- c("DUSP1", "GJA1", "CTGF", "HPGD", "IGFBP2", "MAP1B", "POSTN", "SFRP1", "IGF2", "ICD24")
smc_med12n.genes <- c("GPC3", "STRA6", "ADH1B", "CRABP2", "HSPG2", "LRP1", "DAB2", "SFRP1", "SERPINE1", "APOD")
imm_med12p.genes <- c("GNLY", "RHOB" , "TRBV7-9", "TRBV9", "TRAV38-2DV", "TYROBP", "CCL21", "CD96")
imm_med12n.genes <- c("IGLV3???11", "IGLV3???21", "FCER1G1")

if (nohbbcells) {
  
  #trim the lists & add -log10FDR column
  deLists <- list("deList_smc_med12p"=deList_smc_med12p, "deList_smc_med12n"=deList_smc_med12n, "deList_fib_med12p"=deList_fib_med12p, "deList_fib_med12n"=deList_fib_med12n) #, "deList_endo_med12p"=deList_endo_med12p, "deList_endo_med12n"=deList_endo_med12n
  for (x in 1:length(deLists)) {
    for (c in 1:length(deLists[[x]])) { #goes through as many clusters as each comparison has 
      if (is.null(deLists[[x]][[c]])) {
        print(paste0("skipping delist ", names(deLists[[x]][c]), "")) 
        next
      }
      print(names(deLists[[x]])[c])
      print(paste0(nrow(deLists[[x]][[c]])))
      deLists[[x]][[c]]$avg_logFC <- as.numeric(as.character(deLists[[x]][[c]]$avg_logFC))
      deLists[[x]][[c]] <- deLists[[x]][[c]][order(abs(deLists[[x]][[c]]$avg_logFC), decreasing = T),]
      deLists[[x]][[c]]$p_val_adj <- as.numeric(as.character(deLists[[x]][[c]]$p_val_adj))
      deLists[[x]][[c]] <- deLists[[x]][[c]][which(deLists[[x]][[c]]$p_val_adj < 0.05),]
      deLists[[x]][[c]]$log10pvalAdj <- -log10(deLists[[x]][[c]]$p_val_adj)
      # deLists[[x]][[c]] <- head(deLists[[x]][[c]]) #instead just make it gray
      print(paste0(nrow(deLists[[x]][[c]])))
    }
  }
  
  #now make the volcano data frames & plot
  # setwd("nohbbcells/mes/DE_all/")
  setwd("R_Analysis/DE/stacked_volcanoes")
  for (x in 1:length(deLists)) {
    deFrame <- data.frame()
    for (c in 1:length(deLists[[x]])) {
      print(names(deLists[[x]])[c])
      curframe <- deLists[[x]][[c]]
      curframe$cluster <- "gray"
      curframe$gene <- rownames(curframe)
      if (length(curframe$cluster) < 6) {
        curframe$cluster <- gsub("_DE", "", names(deLists[[x]])[c])
      } else {
        # curframe$cluster[1:6] <- gsub("_DE", "", names(deLists[[x]])[c])  
        #added 2020 09 16 so we can have named, specific volcano plots
        if (grepl("smc", names(deLists)[x])) {
          curframe$cluster[which(rownames(curframe) %in% c(smc_med12p.genes, smc_med12n.genes))] <- 
            gsub("_DE", "", names(deLists[[x]])[c]) 
          print(curframe[which(rownames(curframe) %in% c(smc_med12p.genes, smc_med12n.genes)),])
          # print(paste0("^", paste(smc_med12p.genes,  collapse = "$|^"), "$"))
        } else if (grepl("fib", names(deLists)[x])) {
          curframe$cluster[which(rownames(curframe) %in% c(fib_med12p.genes, fib_med12n.genes))] <- 
            gsub("_DE", "", names(deLists[[x]])[c])  
          # print(paste0("^", paste(fib_med12p.genes,  collapse = "$|^"), "$"))
        } 
      }
      # curframe <- cbind(curframe, "cluster"=gsub("_DE", "", names(deLists[[x]])[c]))
      curframe <- curframe[,c(2,3,4,6,7,8)] #avg_logFC, pct.1, pct.2, log10FDR, cluster, and gene
      if (c==1) {
        deFrame <- curframe
      } else {
        deFrame <- rbind(deFrame, curframe)
      }
    }
    #plot volcanoes
    num_clusts_detected <- length(deLists[[x]]) #length(levels(deFrame$cluster)) - 1
    correct_clusts <- vapply(strsplit(as.character(names(deLists[[x]])), "_"),'[', 2, FUN.VALUE=character(1))
    # if (grepl("fib", names(deLists)[x])) {
    deFrame$cluster <- factor(deFrame$cluster, levels=c(paste("cluster", correct_clusts, sep = "_"), "gray") )
    # } else {
    #   deFrame$cluster <- factor(deFrame$cluster, levels=c(paste("cluster", c(0:(num_clusts_detected - 1)), sep = "_"), "gray") )
    # }
    if (grepl("smc", names(deLists)[x])) { #this is just for naming in the plot
      pop = "Smooth_Muscle"
      print(pop)
      if (grepl("med12p", names(deLists)[x])) {
        comp = "Med12+ vs. Myometrium"
      } else {
        comp = "Med12- vs. Myometrium"
      }
      if (num_clusts_detected <= 3) {
        myColors <- c(brewer.pal(3,"Dark2")[1:num_clusts_detected], LIGHTGRAY )
      } else {
        myColors <- c(brewer.pal(9,"Set1")[1:(num_clusts_detected-3)], brewer.pal(3,"Dark2"), LIGHTGRAY )  
      }
    } else if (grepl("fib", names(deLists)[x])) {
      pop = "Fibroblast"
      print(pop)
      if (grepl("med12p", names(deLists)[x])) {
        comp = "Med12+ vs. Myometrium"
      } else {
        comp = "Med12- vs. Myometrium"
      }
      if (num_clusts_detected <= 3) {
        myColors <- c(brewer.pal(3,"Dark2")[1:num_clusts_detected], LIGHTGRAY )
      } else {
        myColors <- c(brewer.pal(9,"Set1")[1:(num_clusts_detected-3)], brewer.pal(3,"Dark2"), LIGHTGRAY )  
      }
      
    } else if (grepl("endo", names(deLists)[x])) {
      pop = "Endothelial"
      print(pop)
      if (grepl("med12p", names(deLists)[x])) {
        comp = "Med12+ vs. Myometrium"
      } else {
        comp = "Med12- vs. Myometrium"
      }
      if (num_clusts_detected <= 3) {
        myColors <- c(brewer.pal(3,"Dark2")[1:num_clusts_detected], LIGHTGRAY )
      } else {
        myColors <- c(brewer.pal(9,"Set1")[1:(num_clusts_detected-3)], brewer.pal(3,"Dark2"), LIGHTGRAY )  
      }
    }
    
    # print(pop)
    print(levels(deFrame$cluster))
    print(myColors)
    # print(length(levels(deFrame$cluster)))
    # print(length(myColors))
    
    
    names(myColors) <- levels(deFrame$cluster)
    deFrame$label <- NA
    deFrame$label[which(deFrame$cluster != "gray")] <- deFrame[which(deFrame$cluster != "gray"),]$gene
    p <- ggplot2::ggplot(data=deFrame[which(deFrame$cluster=="gray"),], aes(x = avg_logFC, y = log10pvalAdj)) + 
      geom_point(aes(color=cluster), size=5.5) +
      geom_point(data=deFrame[which(deFrame$cluster!="gray"),], aes(x=avg_logFC, y=log10pvalAdj, color = cluster), size=5.5) +
      geom_text_repel(data=deFrame[which(deFrame$cluster!="gray"),], aes(label = label), point.padding = 0.25, size = 5.5) +
      xlab(paste("log2 Fold change ( ",comp, " )", sep="")) + 
      ylab("-log10 FDR adjusted p-value") +
      scale_colour_manual(name = "cluster", values = myColors) +
      labs(title = paste0(pop, " ", comp)) + 
      guides(color=guide_legend(override.aes = list(size=10))) +
      theme(legend.text=element_text(size=12, face = "bold"), legend.key.height = unit(2, 'lines'), 
            panel.background = element_rect(fill = "white", color="black"), 
            axis.title = element_text(face="bold"), title = element_text(face="bold")) 
    
    print(p)
    ggsave(paste0(names(deLists)[x], "_volc.pdf"), dpi = 300)
    
    #and write plot data to table 
    write.table(deFrame, paste0(names(deLists)[x], "_volcDATA.txt"), sep = "\t", quote = F, col.names = NA)
  }
} else {
  
  #trim the lists & add -log10FDR column
  deLists <- list("deList_imm_med12p"=deList_imm_med12p, "deList_imm_med12n"=deList_imm_med12n) #"deList_mes_med12p"=deList_mes_med12p, "deList_mes_med12n"=deList_mes_med12n, 
  for (x in 1:length(deLists)) {
    for (c in 1:length(deLists[[x]])) {
      print(names(deLists[[x]])[c])
      print(paste0(nrow(deLists[[x]][[c]])))
      deLists[[x]][[c]]$avg_logFC <- as.numeric(as.character(deLists[[x]][[c]]$avg_logFC))
      deLists[[x]][[c]] <- deLists[[x]][[c]][order(abs(deLists[[x]][[c]]$avg_logFC), decreasing = T),]
      deLists[[x]][[c]]$p_val_adj <- as.numeric(as.character(deLists[[x]][[c]]$p_val_adj))
      deLists[[x]][[c]] <- deLists[[x]][[c]][which(deLists[[x]][[c]]$p_val_adj < 0.05),]
      deLists[[x]][[c]]$log10pvalAdj <- -log10(deLists[[x]][[c]]$p_val_adj)
      # deLists[[x]][[c]] <- head(deLists[[x]][[c]]) #instead just make it gray
      print(paste0(nrow(deLists[[x]][[c]])))
    }
  }
  
  #now make the volcano data frames & plot
  # setwd("manuscript/fig4/")
  setwd("R_Analysis/DE/stacked_volcanoes")
  for (x in 1:length(deLists)) {
    deFrame <- data.frame()
    for (c in 1:length(deLists[[x]])) {
      print(names(deLists[[x]])[c])
      curframe <- deLists[[x]][[c]]
      curframe$cluster <- "gray"
      if (length(curframe$cluster) < 6) {
        curframe$cluster <- gsub("_DE", "", names(deLists[[x]])[c])
      } else {
        # curframe$cluster[1:6] <- gsub("_DE", "", names(deLists[[x]])[c]) 
        #added 2020 09 16 so we can have named, specific volcano plots
        curframe$cluster[grep(paste0("^", paste(imm_med12p.genes,  collapse = "$|^"), "$"), rownames(curframe))] <- 
          gsub("_DE", "", names(deLists[[x]])[c])  
      }
      # curframe <- cbind(curframe, "cluster"=gsub("_DE", "", names(deLists[[x]])[c]))
      curframe <- curframe[,c(2,3,4,6,7)] #avg_logFC, pct.1, pct.2, log10FDR, cluster
      if (c==1) {
        deFrame <- curframe
      } else {
        deFrame <- rbind(deFrame, curframe)
      }
    }
    #plot
    if (grepl("imm", names(deLists)[x])) {
      deFrame$cluster <- factor(deFrame$cluster, levels=c("cluster_0", "cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6",
                                                          "cluster_7", "cluster_8", "cluster_9", "gray"))
    } else {
      deFrame$cluster <- factor(deFrame$cluster, levels=c("cluster_0", "cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6",
                                                          "cluster_7", "cluster_8", "cluster_9", "cluster_10", "gray"))  
    }
    
    if (grepl("mes", names(deLists)[x])) {
      pop = "Mesenchymal"
      if (grepl("med12p", names(deLists)[x])) {
        comp = "Med12+ vs. Myometrium"
      } else {
        comp = "Med12- vs. Myometrium"
      }
    } else {
      pop = "Immune"
      if (grepl("med12p", names(deLists)[x])) {
        comp = "Med12+ vs. Myometrium"
      } else {
        comp = "Med12- vs. Myometrium"
      }
    }
    
    if (pop == "Mesenchymal") {
      myColors <- c(brewer.pal(8,"Set1"), brewer.pal(3,"Dark2"), LIGHTGRAY )
    } else {
      myColors <- c(brewer.pal(7,"Set1"), brewer.pal(3,"Dark2"), LIGHTGRAY )
    }
    
    print(levels(deFrame$cluster))
    print(myColors)
    print(length(levels(deFrame$cluster)))
    print(length(myColors))
    
    names(myColors) <- levels(deFrame$cluster)
    # deFrame$label <- NA
    # deFrame$label[which(deFrame$cluster != "gray")] <- rownames(deFrame[which(deFrame$cluster != "gray"),])
    p <- ggplot2::ggplot(data=deFrame[which(deFrame$cluster=="gray"),], aes(x = avg_logFC, y = log10pvalAdj)) + 
      geom_point(aes(color=cluster), size=3) +
      geom_point(data=deFrame[which(deFrame$cluster!="gray"),], aes(x=avg_logFC, y=log10pvalAdj, color = cluster), size=3) +
      # geom_text_repel(data=deFrame[which(deFrame$cluster!="gray"),], aes(label = label), point.padding = 0.25, size = 3) +
      xlab(paste("log2 Fold change ( ",comp, " )", sep="")) + 
      ylab("-log10 FDR adjusted p-value") +
      scale_colour_manual(name = "cluster", values = myColors) +
      labs(title = paste0(pop, " ", comp)) + 
      theme(legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'), # panel.grid.major.x = element_line(linetype="dashed",color=LIGHTGRAY), #, panel.grid = element_blank()
            panel.background = element_rect(fill = "white", color="black"), 
            axis.title = element_text(face="bold"), title = element_text(face="bold")) 
    
    print(p)
    ggsave(paste0(names(deLists)[x], "_volc.pdf"))
    
    #and write plot data to table 
    write.table(deFrame, paste0(names(deLists)[x], "_volcDATA.txt"), sep = "\t", quote = F, col.names = NA)
  }
}


#### sd median barchart ~ OK #####
barPlotsWmedian <- function(seurat_object, name, likeUCpaper=TRUE) {
  
  if (likeUCpaper) {
    celltypeNums <- getCellNumbers(seurat_object, celltype_col = "seurat_clusters", sample_col = "orig.ident")
    celltypeNums <- celltypeNums[,c("myometrium061919", "myometrium11911", "myometrium55_11564", "myometrium55_12640","myometrium55_12745",
                                    "fibroid061919", "fibroid55_12382", "fibroid55_12640", "fibroid55_12906", "med12pos55_12843",
                                    "fibroid55_11564", "fibroid55_12716", "fibroid55_12746")]
    celltypeNums.pctOfSmple <- round(t(celltypeNums) / colSums(celltypeNums), 3)*100
    covariates <- data.frame("condition"=c(rep("myometrium", 5), rep("med12pos", 5), rep("med12neg", 3))) #when celltypeNums uses orig.ident
    celltypeNums.wconds <- cbind(celltypeNums.pctOfSmple, covariates)
    celltypeNums.wconds.melt <- melt(celltypeNums.wconds)
    colnames(celltypeNums.wconds.melt) <- c("condition", "cluster",  "percent_of_sample")
    celltypeNums.wconds.melt$cluster <- factor(celltypeNums.wconds.melt$cluster)  
    if (length(unique(celltypeNums.wconds.melt$condition)) == 3)  {
      celltypeNums.wconds.melt$condition <- factor(celltypeNums.wconds.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
    }
    #add stats (mean and median)
    celltypeNums.wconds.melt$median <- NA
    celltypeNums.wconds.melt$mean <- NA
    celltypeNums.wconds.melt$total <- NA
    celltypeNums.wconds.melt$condition <- as.character(celltypeNums.wconds.melt$condition)
    celltypeNums.wconds.melt$cluster <- as.character(celltypeNums.wconds.melt$cluster)
    for (i in unique(celltypeNums.wconds.melt$cluster)) {
      for (j in unique(celltypeNums.wconds.melt$condition)) {
        cur.frame <- celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),]
        median <- median(cur.frame$percent_of_sample)
        mean <- mean(cur.frame$percent_of_sample)
        total <- sum(cur.frame$percent_of_sample)
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"median"] <- median
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"mean"] <- mean
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"total"] <- total
      }
    }
    
    celltypeNums.stats <- summarySEwithin(celltypeNums.wconds.melt, measurevar="percent_of_sample", withinvars=c("condition","cluster"))
    celltypeNums.stats <- celltypeNums.stats[order(celltypeNums.stats$condition,celltypeNums.stats$cluster),]
    celltypeNums.wconds.melt$se <- NA
    for (i in unique(celltypeNums.wconds.melt$cluster)) {
      for (j in unique(celltypeNums.wconds.melt$condition)) {
        cur.frame <- celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),]
        se <- celltypeNums.stats[which(celltypeNums.stats$condition == j & celltypeNums.stats$cluster == i), "se"]
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"se"] <- se
      }
    }
    
    celltypeNums.wconds.melt$condition <- factor(celltypeNums.wconds.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
    ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, "median")
    ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, "mean")
    ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, "total")
    
    write.table(celltypeNums.wconds, paste0("nohbbcells/mes/correct_bar_charts/celltypePctsWconds_", name, ".txt"), quote=F, col.names = NA, sep = "\t")
  }else {
    
    celltypeNums <- getCellNumbers(seurat_object, celltype_col = "seurat_clusters", sample_col = "orig.ident")
    cond.coefs <- sum(celltypeNums)/colSums(celltypeNums)
    celltypeNums.norm <- t(apply(celltypeNums , 1, function(x) return(cond.coefs * x)))
    celltypeNums.norm.pct <- round(celltypeNums.norm / rowSums(celltypeNums.norm), 3)*100
    celltypeNums.norm.pct <- celltypeNums.norm.pct[,c("myometrium061919", "myometrium11911", "myometrium55_11564", "myometrium55_12640","myometrium55_12745",
                                                      "fibroid061919", "fibroid55_12382", "fibroid55_12640", "fibroid55_12906", "med12pos55_12843",
                                                      "fibroid55_11564", "fibroid55_12716", "fibroid55_12746")]
    covariates <- data.frame("condition"=c(rep("myometrium", 5), rep("med12pos", 5), rep("med12neg", 3))) #when celltypeNums uses orig.ident
    celltypeNums.wconds <- cbind(t(celltypeNums.norm.pct), covariates)
    celltypeNums.wconds.melt <- melt(celltypeNums.wconds)
    colnames(celltypeNums.wconds.melt) <- c("condition", "cluster",  "percent_cells_in_celltype")
    celltypeNums.wconds.melt$cluster <- factor(celltypeNums.wconds.melt$cluster)
    if (length(unique(celltypeNums.wconds.melt$condition)) == 3)  {
      celltypeNums.wconds.melt$condition <- factor(celltypeNums.wconds.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
    }
    #add stats (mean and median)
    celltypeNums.wconds.melt$median <- NA
    celltypeNums.wconds.melt$mean <- NA
    celltypeNums.wconds.melt$total <- NA
    celltypeNums.wconds.melt$condition <- as.character(celltypeNums.wconds.melt$condition)
    celltypeNums.wconds.melt$cluster <- as.character(celltypeNums.wconds.melt$cluster)
    for (i in unique(celltypeNums.wconds.melt$cluster)) {
      for (j in unique(celltypeNums.wconds.melt$condition)) {
        cur.frame <- celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),]
        median <- median(cur.frame$percent_cells_in_celltype)
        mean <- mean(cur.frame$percent_cells_in_celltype)
        total <- sum(cur.frame$percent_cells_in_celltype)
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"median"] <- median
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"mean"] <- mean
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"total"] <- total
      }
    }
    
    celltypeNums.stats <- summarySEwithin(celltypeNums.wconds.melt, measurevar="percent_cells_in_celltype", withinvars=c("condition","cluster"))
    celltypeNums.stats <- celltypeNums.stats[order(celltypeNums.stats$condition,celltypeNums.stats$cluster),]
    celltypeNums.wconds.melt$se <- NA
    for (i in unique(celltypeNums.wconds.melt$cluster)) {
      for (j in unique(celltypeNums.wconds.melt$condition)) {
        cur.frame <- celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),]
        se <- celltypeNums.stats[which(celltypeNums.stats$condition == j & celltypeNums.stats$cluster == i), "se"]
        celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"se"] <- se
      }
    }
    
    celltypeNums.wconds.melt$condition <- factor(celltypeNums.wconds.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
    ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, "median")
    # ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, "mean")
    # ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, "total")
    
    write.table(celltypeNums.wconds, paste0("celltypePctsWconds_", name, ".txt"), quote=F, col.names = NA, sep = "\t")
  }
}
ggplotAndSaveBarplots <- function(plotdata, name, type) {
  prevdir <- getwd()
  if (grepl("endo", name)) {
    setwd("nohbbcells/mes/correct_bar_charts/endo")
  } else if (grepl("fib", name)) {
    setwd("nohbbcells/mes/correct_bar_charts/fib")
  } else if (grepl("smc", name)) {
    setwd("nohbbcells/mes/correct_bar_charts/smc")
  } else if (grepl("imm", name)) {
    setwd("nohbbcells/mes/correct_bar_charts/imm")
  } else { print("unable to change dir")}
  
  if (type == "median") {
    ggplot(plotdata, aes(x=cluster, y=median, fill = condition)) + 
      geom_bar(stat="identity", position = position_dodge(width=0.7),  width = 0.6) +
      geom_errorbar(aes(ymin=ifelse(median-se < 0, 0,median-se), ymax=median+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.7)) + 
      scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("Median ", name)) +
      theme(axis.text.x = element_text(vjust = 1), 
            axis.title.x = element_text(margin = margin(t=1,r=0.5,b=1,l=0.5)),
            panel.background = element_rect(fill="white"), 
            panel.grid = element_blank(),
            # panel.grid.major.y = element_line(size=0.5,color=LIGHTGRAY), panel.grid.minor = element_blank(),
            panel.border = element_blank()) + scale_y_continuous(expand = c(0,0))
    ggsave(paste0("barplot_median_", name, ".pdf"), width=12, height=5, dpi=300)
    
    #now make the median box plots
    ggplot(plotdata, aes(x=cluster, y=percent_of_sample, fill = condition)) + 
      geom_boxplot() + 
      scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("Median ", name)) +
      theme(axis.text.x = element_text(vjust = 1), panel.background = element_rect(fill="white"), 
            panel.grid = element_blank(),
            panel.border = element_blank()) + scale_y_continuous(expand = c(0,0))
    ggsave(paste0("boxplot_median_", name, ".pdf"), width=12, height=5, dpi=300)
    
  } else if (type == "mean") {
    #bar plot
    ggplot(plotdata, aes(x=cluster, y=mean, fill = condition)) + 
      geom_bar(stat="identity", position = position_dodge(width=0.7),  width = 0.6) +
      geom_errorbar(aes(ymin=ifelse(mean-se < 0, 0,mean-se), ymax=mean+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.7)) + 
      scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("Mean ", name)) +
      theme(axis.text.x = element_text(vjust = 1), 
            axis.title.x = element_text(margin = margin(t=1,r=0.5,b=1,l=0.5)),
            panel.background = element_rect(fill="white"), 
            panel.grid = element_blank(),
            panel.border = element_blank()) + scale_y_continuous(expand = c(0,0))
    ggsave(paste0("barplot_mean_", name, ".pdf"), width=12, height=5, dpi=300)
    
  } else if (type == "total") {
    ggplot(plotdata, aes(x=cluster, y=total, fill = condition)) + 
      geom_bar(stat="identity", position = position_dodge(width=0.7),  width = 0.6) +
      geom_errorbar(aes(ymin=ifelse(total-se < 0, 0,total-se), ymax=total+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.7)) + 
      scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("total ", name)) +
      theme(axis.text.x = element_text(vjust = 1), 
            axis.title.x = element_text(margin = margin(t=1,r=0.5,b=1,l=0.5)),
            panel.background = element_rect(fill="white"), 
            panel.grid = element_blank(),
            panel.border = element_blank()) + scale_y_continuous(expand = c(0,0))
    ggsave(paste0("barplot_total_", name, ".pdf"), width=12, height=5, dpi=300)
  }
  setwd(prevdir)
}


seurat.objects <- list(all_int_final.noHBBcells.mes.fib, all_int_final.noHBBcells.mes.smc, all_int_final.noHBBcells.mes.endo)
for (o in seurat.objects) {
  name <- o@project.name
  print(name)
  barPlotsWmedian(o, name, likeUCpaper = TRUE)
}

