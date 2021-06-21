# library(Seurat)
# library(SingleR)
# library(ggplot2)
HEART.data <- Read10X_h5("HEART.h5")
HEART <- CreateSeuratObject(counts = HEART.data, project = "IMMUNE_stim", min.cells = 5)
HEART <- NormalizeData(object = HEART, verbose = FALSE)
HEART <- FindVariableFeatures(object = HEART, selection.method = "vst", nfeatures = 2000)
immune.combined=HEART
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern = "^mt-")
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
immune.combined <- subset(immune.combined, 
                          subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 3000 & 
                            percent.mt < 5)
dim(immune.combined)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#3.3 
immune.combined <- RunPCA(immune.combined, npcs =20, verbose = FALSE)
#3.4 
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20,seed.use = 3
)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)
DimPlot(immune.combined, reduction = "umap",label =F)
#3.5 Find markers
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(immune.combined,immune.combined.markers,cellchat,file="FOR-ZX20210606-2.Rdata")

immune.combined_for_SingleR <- GetAssayData(immune.combined, slot="data")
clusters=immune.combined@meta.data$seurat_clusters
#3.6 ANNOTATION
load("F:/XU/20210119/20210119/ref_Mouse_imm.RData")
mouseImmu <- ref_IGD
load("F:/XU/20210119/20210119/ref_Mouse_all.RData")
mouseRNA <- ref_Mouse
pred.mouseImmu <- SingleR(test = immune.combined_for_SingleR, ref = mouseImmu, labels = mouseImmu$label.fine,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.mouseRNA <- SingleR(test = immune.combined_for_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType=data.frame(ClusterID=levels(immune.combined@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels )
head(cellType)
immune.combined@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'mouseRNA']
immune.combined@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'mouseImmu']
table(immune.combined@active.ident)

#markers.to.plot=c("Cd1c","Clec9a","Lilra4","Cd11c","Pcsk9","Cst3","Fcer1g","Lyz")
markers.to.plot=c("Cd83","Ccr7","H2-Eb1","H2-Aa","Cd80","Cd86","Cd28","Cd40","Cd40lg")
FeaturePlot(immune.combined, features = "Pcsk9" ,label = T, min.cutoff = "q9")

VlnPlot(immune.combined, features = markers.to.plot, slot = "counts", log =F, pt.size = 0
        #,cols = col
)
markers.to.plot=c("Cd3g","Cd4","Cd8a","Gzmk","Cd79a","Clec9a","Cst3", "C1qc", "Lyz2", "S100a8","Ncr1")
DotPlot(immune.combined, features = markers.to.plot,dot.scale = 10) + RotatedAxis()

#immune.combined <- subset(immune.combined, idents =c(2:11))#get subset


# 
# immune.combined <- RenameIdents(immune.combined,`0` = "C0_macrophage")#
# immune.combined <- RenameIdents(immune.combined,`1` = "C1_macrophage")#
# immune.combined <- RenameIdents(immune.combined,`2` = "C2_NK")#
# immune.combined <- RenameIdents(immune.combined,`3` = "C3_T")#
# immune.combined <- RenameIdents(immune.combined,`4` = "C4_T")#
# immune.combined <- RenameIdents(immune.combined,`5` = "C5_Granulocytes")#
# immune.combined <- RenameIdents(immune.combined,`6` = "C6_T")#
# immune.combined <- RenameIdents(immune.combined,`7` = "C7_T")#
# immune.combined <- RenameIdents(immune.combined,`8` = "C8_T")#
# immune.combined <- RenameIdents(immune.combined,`9` = "C9_Monocytes")#
# immune.combined <- RenameIdents(immune.combined,`10` = "C10_DC")#
# immune.combined <- RenameIdents(immune.combined,`11` = "C11_B")#
# immune.combined <- RenameIdents(immune.combined,`12` = "C12_T")#

top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(immune.combined, features = top10$gene) + NoLegend ()
#Comparison analysis of multiple datasets using CellChat
library(CellChat)
library(patchwork)

library(Seurat)
library(SeuratData)
library(cowplot)
library(sdtoolkit)
library(patchwork)
library(magrittr)
library(dplyr)
library(SingleR)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(tidyverse)
library(ggalluvial)

DimPlot(immune.combined, reduction = "umap",label = F)#带有数字标记的分???,即F1的p2
Idents(immune.combined) <- "stim"
immune.combined <- subset(immune.combined, idents =c("HEART"))#get subset
Idents(immune.combined) <- "singleR"
DimPlot(immune.combined, reduction = "umap",label = F)#带有数字标记的分???,即F1的p2
#immune.combined$cellclustername=immune.combined@active.ident
cellchat <- createCellChat(object = immune.combined, group.by = "ident")
groupSize <- as.numeric(table(cellchat@idents)) # 后面有用

##设置参考数据库
# 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.mouse  
# 使用"Secreted Signaling"用于细胞通讯分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# table(CellChatDB$interaction$annotation)
# Cell-Cell Contact       ECM-Receptor Secreted Signaling 
# 378                432               1211 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 


# 将数据库传递给cellchat对象
cellchat@DB <- CellChatDB.use 
##配体-受体分析
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
##推测细胞通讯网络
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#
save(cellchat,groupSize,file="cellchat-groupSize-20210606Cell-Cell-Contact.Rdata")
#细胞通讯网络系统分析及可视化
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, signaling = "MHC-I", color.heatmap = "Reds")

levels(cellchat@idents)            #查看细胞顺序
vertex.receiver = c(1, 6)          #指定靶细胞的索引
vertex.receiver = c(1,2,3,4,5)          #指定靶细胞的索引
vertex.receiver = c(1,2,3,4)

cellchat@netP$pathways             #查看富集到的信号通路
[1] "ALCAM"       "APP"         "BST2"        "CD200"       "CD22"       
[6] "CD226"       "CD39"        "CD40"        "CD45"        "CD48"       
[11] "CD52"        "CD6"         "CD80"        "CD86"        "CD96"       
[16] "CDH"         "CDH1"        "CLEC"        "ICAM"        "ITGAL-ITGB2"
[21] "JAM"         "LAIR1"       "LCK"         "MHC-I"       "MHC-II"     
[26] "PD-L1"       "PDL2"        "PVR"         "SELPLG"      "SEMA4"      
[31] "SEMA7"       "THY1"        "TIGIT"

pathways.show <- "MHC-I"             #指定需要展示的通路
netVisual_aggregate(cellchat, signaling = "CD40",
                    # vertex.receiver = c(5,6)
                    , layout = "circle"
                    , vertex.size = groupSize)
ccc=cellchat@netP$pathways
# for (i in 1:18) 
# {pp=netVisual_aggregate(cellchat, signaling = ccc[i], layout = "circle", 
#       vertex.size = groupSize);print(pp);
# ggplot2::ggsave(filename = paste0(ccc[i],'_PATHWAY.pdf'))}
for (i in 1:18) 
{netVisual_aggregate(cellchat, signaling = ccc[i], layout = "circle", 
                     vertex.size = groupSize)}



# Hierarchy plot
png(filename = "sig_pathway_hierarchy.png", width = 1000, height = 650)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
dev.off()
# Circle plot
png(filename = "sig_pathway_cricle.png", width = 650, height = 600)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
dev.off()
# 计算配体-受体对信号网络的贡献???
png(filename = "sig_pathway_L-R.png", width = 800, height = 600)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
# 分析细胞在信号网络中角色

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")#######
netAnalysis_signalingRole_network(cellchat, slot.name = "netP")#######
# png(filename = "sig_pathway_role.png", width = 800, height = 600)
# netVisual_signalingRole(cellchat, signaling = pathways.show)
# dev.off()


##细胞通讯模式和信号网???
nPatterns = 5   #
cellchat <- identifyCommunicationPatterns(cellchat, slot.name = "netP", pattern = "outgoing", k = nPatterns)
# river plot
p = netAnalysis_river(cellchat, pattern = "outgoing")
ggsave("com_pattern_outgoing_river.png", p, width = 12, height = 6)
# dot plot
p = netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave("com_pattern_outgoing_dot.png", p, width = 9, height = 6)

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
p = netAnalysis_river(cellchat, pattern = "incoming")
ggsave("com_pattern_incoming_river.png", p, width = 12, height = 6)
# dot plot
p = netAnalysis_dot(cellchat, pattern = "incoming")
ggsave("com_pattern_incoming_dot.png", p, width = 9, height = 6)


##信号网络聚类
# 按功能相似性聚???

cellchat <- computeNetSimilarity(cellchat, type = "functional")
# Manifold learning of the signaling networks for a single dataset 
# Error in runUMAP(Similarity, min.dist = 0.3, n.neighbors = k) : 
#   Cannot find UMAP, please install through pip 
# (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).
#immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20,seed.use = 2)
cellchat <- netEmbedding(cellchat, type = "functional")#
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
p = netVisual_embedding(cellchat, type = "functional")
ggsave("custer_pathway_function.png", p, width = 9, height = 6)
p = netVisual_embeddingZoomIn(cellchat, type = "functional")
ggsave("custer_pathway_function2.png", p, width = 8, height = 6)

# 按结构相似性聚???
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
p = netVisual_embedding(cellchat, type = "structural")
ggsave("custer_pathway_structure.png", p, width = 9, height = 6)
p = netVisual_embeddingZoomIn(cellchat, type = "structural")
ggsave("custer_pathway_structure2.png", p, width = 8, height = 6)

save(cellchat, file = "cellchat.rds")
