# library("DESeq2")
#读取基因表达值矩阵
exprSet <- read.csv("exprSet.csv", row.names=1) 
sampletable <- read.csv("sampletable.csv", row.names=1)

#PCA
#推荐使用 log 转化后的基因表达值，降低不同基因表达水平数量级相差过大的问题
gene=log(exprSet+1)
gene <- t(gene)

gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)

colnames(exprSet)
colnames(sampletable)
sampletable$group
View(exprSet)
View(sampletable)
rownames(sampletable)==colnames(exprSet)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = sampletable,
                              design = ~group)
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)
dds <- DESeq(dds)
save(dds,exprSet,sampletable,file="A1.Rdata")
#a 火山图
nrDEG<- results(dds, contrast=c("group","A","B"),lfcThreshold=1,alpha = 0.05)
nrDEG=as.data.frame(nrDEG)
colnames(nrDEG)
# colnames(nrDEG)=c("baseMean","logFC","lfcSE","stat","pvalue","adj.P.Val")
# colnames(nrDEG)
nrDEG = nrDEG[,c(2,6)]
colnames(nrDEG)=c("logFC","adj.P.Val")
colnames(nrDEG)
nrDEG$group <- as.factor(ifelse(nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC) > 1,ifelse(nrDEG$logFC > 1,'up-regulated','down-regulated'),'not-significant'))
markers=c("Cd40","Cd80","Cd86")
nrDEG$sign <- ifelse(rownames(nrDEG)%in%markers,rownames(nrDEG),NA)
nrDEG$adj.P.Val <- ifelse(nrDEG$adj.P.Val < 1e-150,1e-150,nrDEG$adj.P.Val)
ggplot( nrDEG, aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point(alpha=0.6, size = 1.5) +
  coord_cartesian(xlim =c(-12, 12), ylim = c(0, 150))+#
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept=1 ,linetype=4) +
  geom_vline(xintercept=c(-0.25,0.25) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "blue", "grey90"), limits = c("up-regulated", "down-regulated", "not-significant")) +
  geom_label_repel(aes(label=sign), fontface="bold", color="black",  max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                   box.padding=unit(0.5, "lines"), point.padding=unit(0.25, "lines"), 
                   segment.colour = "grey50")


aa <- nrDEG[which((nrDEG$group== 'down-regulated')), ]
gene <- rownames(aa)
gene.df <- bitr(gene, fromType ="SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.5)#######
#head(kk)[,1:6]

dotplot(kk, showCategory=15) 
kktable=as.data.frame(kk)
write.csv(kktable, file = "CD8-DOWN-KEGG1.csv")

#####################################GSEA##############
# Molecular Signatures Database

msigdbr_show_species() #支持的物种
#"Homo sapiens"             "Mus musculus"            
Dm_msigdbr <- msigdbr(species="Mus musculus")
head(Dm_msigdbr, 2) %>% as.data.frame
DmGO <- msigdbr(species="Mus musculus",category="C2") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
head(DmGO)

geneList=nrDEG$logFC
names(geneList)=rownames(nrDEG)
geneList = sort(geneList, decreasing = TRUE)


#GSEA分析
KEGG <- GSEA(geneList,TERM2GENE=DmGO[,c(1,3)])
#KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
dotplot(KEGG) #出点图 
dotplot(KEGG,color="pvalue")  #按p值出点图
dotplot(KEGG,split=".sign")+facet_grid(~.sign) #出点图，并且分面激活和抑制
gseaplot2(KEGG,1,color="red",pvalue_table = T)

KEGG<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 0.05,organism = "mmu")
ccc=as.data.frame(KEGG)
write.csv(ccc, file = "GSEA.csv")

gseaplot2(KEGG,'LINDSTEDT_DENDRITIC_CELL_MATURATION_B',color="red",pvalue_table = T,)
gseaplot2(KEGG,'LINDSTEDT_DENDRITIC_CELL_MATURATION_A',color="red",pvalue_table = T,)
BIOCARTA_NFKB_PATHWAY
gseaplot2(KEGG,'BIOCARTA_NFKB_PATHWAY',color="red",pvalue_table = T,)