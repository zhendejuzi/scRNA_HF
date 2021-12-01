#######################################################
#########integrating the HF samples####################
#######################################################
#R-4.1.0
library("Seurat")
library("plyr")
library("ggplot2")
library("patchwork")
DF.data <- Read10X(data.dir = "./")
DF <- CreateSeuratObject(counts = DF.data, project = "HF", min.cells = 10, min.features = 200)
DF$group=ldply(strsplit(colnames(DF),split="-"))[,2]
DF[["percent.mt"]] <- PercentageFeatureSet(object = DF, pattern = "^MT-")
DF_QC <- subset(x = DF, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 10 )
ggplot(DF_QC@meta.data,aes(x=nFeature_RNA,y=percent.mt))+geom_point(size=1)+facet_wrap(~group)+theme_bw()
DF.list <- SplitObject(DF_QC, split.by = "group")
for(i in names(DF.list)){DF.list[[i]] <- SCTransform(DF.list[[i]])}
DF.features <- SelectIntegrationFeatures(object.list = DF.list, nfeatures = 1800)
DF.list <- PrepSCTIntegration(object.list = DF.list, anchor.features = DF.features)
DF.anchors <- FindIntegrationAnchors(object.list = DF.list, normalization.method = "SCT",anchor.features = DF.features)
DF.integrated <- IntegrateData(anchorset = DF.anchors, normalization.method = "SCT",verbose = FALSE, dims = 1:30)  # but do not run the ScaleData function after integration. for SCT
DF.integrated <- RunPCA(DF.integrated, verbose = FALSE, npcs = 50)
ElbowPlot(object = DF.integrated,ndim=50)
npc=30
DF.integrated <- RunUMAP(DF.integrated, dims = 1:npc,return.model=TRUE)
DF.integrated <- FindNeighbors(object = DF.integrated, dims = 1:npc,k.param=20)
DF.integrated <- FindClusters(object = DF.integrated, resolution = 0.3,n.start=100)
DimPlot(object = DF.integrated, reduction = "umap",label = TRUE,label.size=7,pt.size = 0.8)

#######################################################
#########extract subset cells for subcluster (PR)######
#######################################################
DF=subset(DF.integrated,ident=5)
DefaultAssay(DF)="RNA"
DF.list<-SplitObject(DF, split.by = "group")
for(i in names(DF.list)){DF.list[[i]] <- SCTransform(DF.list[[i]],method = "glmGamPoi")}
DF.features <- SelectIntegrationFeatures(object.list = DF.list, nfeatures = 600)
DF.list <- PrepSCTIntegration(object.list = DF.list, anchor.features = DF.features)
DF.anchors <- FindIntegrationAnchors(object.list = DF.list, normalization.method = "SCT",anchor.features = DF.features)
DF.integrated <- IntegrateData(anchorset = DF.anchors, normalization.method = "SCT",verbose = FALSE, dims = 1:30) 
DF.integrated <- RunPCA(DF.integrated, verbose = FALSE, npcs = 50)
ElbowPlot(object = DF.integrated,ndim=50)
npc=15
DF.integrated <- RunUMAP(DF.integrated, dims = 1:npc,return.model=TRUE)
DF.integrated <- FindNeighbors(object = DF.integrated, dims = 1:npc,k.param=20)
DF.integrated <- FindClusters(object = DF.integrated, resolution = 0.05,n.start=100)
DimPlot(object = DF.integrated, reduction = "umap",label = TRUE,label.size=7,pt.size = 0.8)

#######################################################
#########extract subset cells for subcluster (Immu)######
#######################################################
DF=subset(DF.integrated,ident=6)
DefaultAssay(DF)="RNA"
DF.list<-SplitObject(DF, split.by = "group")
for(i in names(DF.list)){DF.list[[i]] <- SCTransform(DF.list[[i]],method = "glmGamPoi")}
DF.features <- SelectIntegrationFeatures(object.list = DF.list, nfeatures = 100)
DF.list <- PrepSCTIntegration(object.list = DF.list, anchor.features = DF.features)
DF.anchors <- FindIntegrationAnchors(object.list = DF.list, normalization.method = "SCT",anchor.features = DF.features,k.score=7,dims = 1:7)
DF.integrated <- IntegrateData(anchorset = DF.anchors, normalization.method = "SCT",verbose = FALSE, dims = 1:7,k.weight = 7)  
DF.integrated <- RunPCA(DF.integrated, verbose = FALSE, npcs = 50)
ElbowPlot(object = DF.integrated,ndim=50)
npc=15
DF.integrated <- RunUMAP(DF.integrated, dims = 1:npc,return.model=TRUE)
DF.integrated <- FindNeighbors(object = DF.integrated, dims = 1:npc,k.param=20)
DF.integrated <- FindClusters(object = DF.integrated, resolution = 0.05,n.start=100)
DimPlot(object = DF.integrated, reduction = "umap",label = TRUE,label.size=7,pt.size = 0.8)
K=list();for(i in 1:3){K[[i]]<-FindMarkers(DF.integrated, ident.1 = i-1,only.pos=TRUE ,logfc.threshold=0.25,min.pct=0.25,min.diff.pct=0.25)}
enrichGO.list=list();for(i in 1:3){enrichGO.list[[i]]=enrichGO(rownames(K[[i]]),ont="BP",OrgDb=org.Hs.eg.db,keyType = "SYMBOL")}

library("SingleR")
library("celldex")
ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(test=expr, ref=ref, labels=ref$label.main)
pred$labels=paste("Immu",sam$ident,sep="")
plotScoreHeatmap(pred)
table(pred$pruned.labels,sam$ident)

###################################################
#####integrating the F31 data to reference#########
###################################################
DF.data <- Read10X(data.dir = "./")
DF <- CreateSeuratObject(counts = DF.data, project = "HF", min.cells = 10, min.features = 200)
DF$group=ldply(strsplit(colnames(DF),split="-"))[,2]
DF[["percent.mt"]] <- PercentageFeatureSet(object = DF, pattern = "^MT-")
DF_QC <- subset(x = DF, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 10 )
ggplot(DF_QC@meta.data,aes(x=nFeature_RNA,y=percent.mt))+geom_point(size=1)+facet_wrap(~group)+theme_bw()
F31.list <- SplitObject(DF_QC, split.by = "group")
for(i in names(F31.list)){F31.list[[i]] <- SCTransform(F31.list[[i]])}
for(i in 1:length(F31.list)){
  anchors <- FindTransferAnchors(reference = HF4.integrated, query = F31.list[[i]],dims = 1:30, normalization.method = "SCT",reference.reduction = "pca",recompute.residuals=FALSE,n.trees = 200,nn.method="rann")
  F31.list[[i]]<- MapQuery(anchorset = anchors, reference = HF4.integrated, query = F31.list[[i]],refdata =HF4.integrated.r3$seurat_cluster, reference.reduction = "pca", reduction.model = "umap")
  F31.list[[i]]@active.ident=as.factor(as.numeric(F31.list[[i]]$predicted.id))
  names(F31.list[[i]]@active.ident)=colnames(F31.list[[i]])
}
p1=DimPlot(object = HF4.integrated, reduction = "umap",label = TRUE,label.size=5,pt.size = 0.7)
p2=DimPlot(object = F31.list[[1]],reduction = "ref.umap",label = TRUE,label.size=5,pt.size = 0.7)
p3=DimPlot(object = F31.list[[2]],reduction = "ref.umap",label = TRUE,label.size=5,pt.size = 0.7)
p1 + p2 + p3 + plot_layout(guides = "collect")

#########################################################################
#############find geneMarkers for cellType###############################
#########################################################################
#pwd /local/data/wusijie/SingleCellProject/HF_BamCR31/ZL_HF_F31pair/outs/filtered_feature_bc_matrix
library("scran")
library("clusterProfiler")
library("org.Hs.eg.db")
markers=findMarkers(DF@assays$SCT@data,DF@active.ident,block=DF$group,direction="up",pval.type="some",test.type="wilcox")
enrichGO.list=list()
enrichGO.QC=list()
g=list()
markers.QC=list()
for(i in 1:length(markers)){
  print(i)
  markers[[i]]$mean=apply(markers[[i]][,-c(1:3)],1,mean)
  markers[[i]]$pct=avgPCT[rownames(markers[[i]]),i]}
  markers.QC[[i]]=markers[[i]][markers[[i]][,"mean"]>0.7 & markers[[i]][,"FDR"]<0.01,]
  enrichGO.list[[i]]=enrichGO(rownames(markers.QC[[i]])[1:50],ont="BP",OrgDb=org.Hs.eg.db,keyType = "SYMBOL")
  enrichGO.QC[[i]]=simplify(enrichGO.list[[i]], cutoff = 0.7, by = "p.adjust")
  print(head(enrichGO.QC[[i]]))
}
#########################################################################
#############find geneMarkers for SC#####################################
#########################################################################
library("scran")
library("clusterProfiler")
library("org.Hs.eg.db")
load("analysis_SC.genemarker.Rdata")
DF=subset(HF6,ident=c("0","3","5","S0","S1","S2","S3"))
markers=list()
for(i in 4:7){
  subD=subset(DF,ident=c("0","3","5",levels(DF@active.ident)[i]))
  markers[[i-3]]=findMarkers(subD@assays$SCT@data,subD@active.ident,block=subD$group,direction="both",pval.type="all",test.type="wilcox")[[4]]
}
for(i in 1:length(markers)){
  print(i)
  markers[[i]]$mean=apply(markers[[i]][,c(4:6)],1,mean)
  markers[[i]]$pct.diff=avgPCT[rownames(markers[[i]]),i+3]-apply(avgPCT[rownames(markers[[i]]),1:3],1,max)
  markers[[i]]$pct=avgPCT[rownames(markers[[i]]),i+3]
}
for(i in 1:4){print(table(markers[[i]]$mean>0.8 & markers[[i]]$FDR<0.05, markers[[i]]$pct.diff>0.3))}
index=list();for(i in 1:4){index[[i]]=which(markers[[i]]$mean>0.8 & markers[[i]]$FDR<0.05 & markers[[i]]$pct.diff>0.3)}
g=unique(names(unlist(index)))
DotPlot2(object = HF6, features = g,assay="SCT",scale=FALSE)
###GSEA enrichment###
gmtfile <- system.file("extdata", "mouse.HF.symbols.gmt", package="clusterProfiler")
mm10<-read.gmt(gmtfile)]
mm10.GSEA=list()
gmtfile <- system.file("extdata", "h.all.v7.4.symbols.gmt", package="clusterProfiler")
hall<-read.gmt(gmtfile)
hall.GSEA=list()
for(i in 1:length(markers)){
  k=markers[[i]][markers[[i]]$pct>0.25,]
  k2=k[order(k$mean,decreasing=T),"mean"]
  mm10.GSEA[[i]]=GSEA(scale(k2)[,1], TERM2GENE=mm10,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(mm10[,1]),]
  hall.GSEA[[i]]=GSEA(scale(k2)[,1], TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
}

########################################################################
#############SCENIC for regulon analysis################################
########################################################################
library(SCENIC)
library(Seurat)
library(RcisTarget)
library(visNetwork)
library(data.table)
library(AUCell)
org="hgnc"
dbDir="/picb/dermatogenomics/software/refgenome"
myDatasetTitle="HF"
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20)
exprMat=as.matrix(single.DF@assays$RNA@counts)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=0.03*ncol(exprMat),minSamples=ncol(exprMat)*0.03)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulons=readRDS("int/3.1_regulons_forAUCell.Rds")
regulons=regulons[rownames(regulonAUC)]
mAUC <- getAUC(regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),])
DF[["regulon"]] <- CreateAssayObject(data = mAUC)
DefaultAssay(DF)="regulon"
TFcor=read.table("TFcor.txt",header=T,stringsAsFactors = F)
TF=rownames(TFcor[TFcor$cor_median>0.25 & TFcor$TFmaxExpr>0.6,])
regulon.markers=findMarkers(DF@assays$regulon@data[TF,],DF@active.ident,block=DF$group,direction="up",pval.type="all")
TF.markers=findMarkers(DF@assays$SCT@data[TF,],DF@active.ident,block=DF$group,direction="up",pval.type="all")
for(i in 1:8){regulon.markers[[i]]$TF.pvalue=TF.markers[[i]][rownames(regulon.markers[[i]]),"p.value"]}
for(i in 1:8){regulon.markers[[i]]$TF.expr=avgExpr[rownames(regulon.markers[[i]]),i]}
for(i in 1:8){regulon.markers[[i]]$TF.pct=avgPCT.TF[rownames(regulon.markers[[i]]),i]}
for(i in 1:8){regulon.markers[[i]]$TF.pct.diff=regulon.markers[[i]]$TF.pct-apply(avgPCT.TF[rownames(regulon.markers[[i]]),-i],1,mean)}
for(i in 1:8){regulon.markers[[i]]$TFlogFCmean=apply(TF.markers[[i]][rownames(regulon.markers[[i]]),c(4:10)],1,mean)}
for(i in 1:8){index=which(regulon.markers[[i]]$FDR<0.01 & regulon.markers[[i]]$TF.pvalue<0.05 & regulon.markers[[i]]$TF.expr>0.6 & regulon.markers[[i]]$TF.pct.diff>0.05);print(levels(DF@active.ident)[i]);sig.regulon.markers[[i]]=regulon.markers[[i]][index,]}

########################################################################
######################Black vs White gene###############################
########################################################################
HF6=readRDS("HF6.merge.Rdata")
DefaultAssay(HF6)="SCT"
F31=HF6[,colnames(HF6)[HF6$group=="F31B" | HF6$group=="F31W"]]
F62=HF6[,colnames(HF6)[HF6$group=="F62B" | HF6$group=="F62W"]]
F31.WvB.gene=list()
F62.WvB.gene=list()
Both.WvB.gene=list()
F31.gsea.hall=list()
F62.gsea.hall=list()
Both.gsea.hall=list()
for(i in 1:12){
  print(i)
  DF.F31=subset(F31,ident=i-1)
  DF.F31@active.ident=as.factor(as.numeric(as.factor(DF.F31$group)))
  names(DF.F31@active.ident)=colnames(DF.F31)
  DF.F31$project="F31"
  DF.F62=subset(F62,ident=i-1)
  DF.F62@active.ident=as.factor(as.numeric(as.factor(DF.F62$group)))
  names(DF.F62@active.ident)=colnames(DF.F62)
  DF.F62$project="F62"
  DF=merge(DF.F31,DF.F62)
  fc=0.05;pct=0.2;pct_diff=0;
  F31.WvB.gene[[i]]<-FindMarkers(DF.F31, ident.1 = 2,ident.2=1, logfc.threshold=fc,min.pct=pct,min.diff.pct=pct_diff,min.cells.group = 1)
  g1=F31.WvB.gene[[i]][,2];names(g1)=rownames(F31.WvB.gene[[i]]);g1=sort(g1,decreasing=T);print(length(g1))
  F31.gsea.hall[[i]] <- GSEA(g1, TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
  F62.WvB.gene[[i]]<-FindMarkers(DF.F62, ident.1 = 2,ident.2=1, logfc.threshold=fc,min.pct=pct,min.diff.pct=pct_diff,min.cells.group = 1)
  g2=F62.WvB.gene[[i]][,2];names(g2)=rownames(F62.WvB.gene[[i]]);g2=sort(g2,decreasing=T);print(length(g2))
  F62.gsea.hall[[i]] <- GSEA(g2, TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
}

##compare to aging
WvB<-F31.gsea.hall<-F62.gsea.hall<-OY1.gsea.hall<-OY2.gsea.hall<-F31.gsea.WvB<-OY1.gsea.WvB<-OY2.gsea.WvB<-OY1.gene<-OY2.gene<-F62.WvB.DEG.up<-F62.WvB.DEG.dn<-list()
OY1=HF[,colnames(HF)[HF$group=="F18" | HF$group=="F59"]];OY1$project="OY1"
OY2=HF[,colnames(HF)[HF$group=="F31B" | HF$group=="F62B"]];OY2$project="OY2"
for(i in 1:12){
  print(i)
  DF.F31=subset(F31,ident=i-1)
  DF.F31@active.ident=as.factor(as.numeric(as.factor(DF.F31$group)))
  names(DF.F31@active.ident)=colnames(DF.F31)
  DF.F31$project="F31"
  DF.F62=subset(F62,ident=i-1)
  DF.F62@active.ident=as.factor(as.numeric(as.factor(DF.F62$group)))
  names(DF.F62@active.ident)=colnames(DF.F62)
  DF.F62$project="F62"
  DF=merge(DF.F31,DF.F62)
  DF.OY1=subset(OY1,ident=i-1);
  DF.OY1@active.ident=as.factor(as.numeric(as.factor(DF.OY1$group)))
  names(DF.OY1@active.ident)=colnames(DF.OY1)
  DF.OY2=subset(OY2,ident=i-1);
  DF.OY2@active.ident=as.factor(as.numeric(as.factor(DF.OY2$group)))
  names(DF.OY2@active.ident)=colnames(DF.OY2)
  
  fc=0.05;pct=0.2;pct_diff=0;
  F31.WvB.gene[[i]]<-FindMarkers(DF.F31, ident.1 = 2,ident.2=1, logfc.threshold=fc,min.pct=pct,min.diff.pct=pct_diff,min.cells.group = 1)
  g1=F31.WvB.gene[[i]][,2];names(g1)=rownames(F31.WvB.gene[[i]]);g1=sort(g1,decreasing=T);print(length(g1))
  F31.gsea.hall[[i]] <- GSEA(g1, TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
  F62.WvB.gene[[i]]<-FindMarkers(DF.F62, ident.1 = 2,ident.2=1, logfc.threshold=fc,min.pct=pct,min.diff.pct=pct_diff,min.cells.group = 1)
  g2=F62.WvB.gene[[i]][,2];names(g2)=rownames(F62.WvB.gene[[i]]);g2=sort(g2,decreasing=T);print(length(g2))
  F62.gsea.hall[[i]] <- GSEA(g2, TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
  
  F62.WvB.DEG.up[[i]]<-FindMarkers(DF.F62, ident.1 = 2,ident.2=1,only.pos=TRUE ,logfc.threshold=0.25,min.pct=0.25,min.diff.pct=0)
  F62.WvB.DEG.dn[[i]]<-FindMarkers(DF.F62, ident.1 = 1,ident.2=2,only.pos=TRUE ,logfc.threshold=0.25,min.pct=0.25,min.diff.pct=0)
  k1=rownames(F62.WvB.DEG.up[[i]])[1:50]
  k2=rownames(F62.WvB.DEG.dn[[i]])[1:50]
  k1=data.frame(term=rep("WvB_UP",length(k1)),gene=k1)
  k2=data.frame(term=rep("WvB_DN",length(k2)),gene=k2)
  WvB[[i]]=rbind(k1,k2)
  print(table(WvB[[i]][,1]))
  F31.gsea.WvB[[i]]=GSEA(g1, TERM2GENE=WvB[[i]],pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[c("WvB_UP","WvB_DN"),]
  print(head(F31.gsea.WvB[[i]]))
  OY1.gene[[i]]<-FindMarkers(DF.OY1, ident.1 = 2,ident.2=1, logfc.threshold=fc,min.pct=pct,min.diff.pct=pct_diff,min.cells.group = 1)
  g4=OY1.gene[[i]][,2];names(g4)=rownames(OY1.gene[[i]]);g4=sort(g4,decreasing=T);print(length(g4))
  OY1.gsea.WvB[[i]] <- GSEA(g4, TERM2GENE=WvB[[i]],pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[c("WvB_UP","WvB_DN"),]
  print(head(OY1.gsea.WvB[[i]]))
  OY1.gsea.hall[[i]] <- GSEA(g4, TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
  OY2.gene[[i]]<-FindMarkers(DF.OY2, ident.1 = 2,ident.2=1, logfc.threshold=fc,min.pct=pct,min.diff.pct=pct_diff,min.cells.group = 1)
  g5=OY2.gene[[i]][,2];names(g5)=rownames(OY2.gene[[i]]);g5=sort(g5,decreasing=T);print(length(g5))
  OY2.gsea.WvB[[i]] <- GSEA(g5, TERM2GENE=WvB[[i]],pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[c("WvB_UP","WvB_DN"),]
  print(head(OY2.gsea.WvB[[i]]))
  OY2.gsea.hall[[i]] <- GSEA(g5, TERM2GENE=hall,pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)[levels(hall[,1]),]
  print("#######")
}

###############################################################
#########trajectory analysis using monocle3####################
###############################################################
#pwd /local/data/wusijie/SingleCellProject/HF_BamCR31/ZL_HF_F31pair/outs/filtered_feature_bc_matrix
library("monocle3")
DF=readRDS("HF6.merge.Rdata")
DF=subset(DF,ident=c(0:9))
DimPlot(object = DF, reduction = "umap",label = TRUE,label.size=5,pt.size = 0.7)
pd2=DF@meta.data[,1:10]   #data frame
fd2 <- data.frame(gene_short_name = row.names(DF), row.names = row.names(DF))  #dataframe
cds <- new_cell_data_set(as(DF@assays$SCT@counts, "sparseMatrix"),cell_metadata = pd2,gene_metadata  = fd2)
sam<-FetchData(DF,vars=c("ident","group","UMAP_1","UMAP_2"))
reducedDims(cds)$UMAP<-sam[,c("UMAP_1","UMAP_2")]
colnames(reducedDims(cds)$UMAP)=c("V1","V2")
cds <- cluster_cells(cds,partition_qval=0.95,k=60)
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds,learn_graph_control=list(minimal_branch_len=2,euclidean_distance_ratio=10))
cds <- order_cells(cds)
plot_cells(cds,color_cells_by = "pseudotime",label_leaves=FALSE,label_roots=FALSE,trajectory_graph_color="black",label_branch_points=FALSE, trajectory_graph_segment_size = 0.5,cell_size=1.5)
save.image("analysis_HF6.monocle3.Rdata")

###############################################################
#########trajectory analysis using destiny#####################
###############################################################
library("destiny")
DM <- DiffusionMap(MxTA_pc)
plot3d(eigenvectors(DM)[,1:3],col=colors[as.numeric(as.character(sam$ident))+1],pch=19,box=F,axes=T,size =4)
dpt_DM <- DPT(DM)

###############################################################
#########trajectory analysis using palantir#####################
###############################################################
pwd /local/data/wusijie/SingleCellProject/HF_BamCR31/ZL_HF_F31pair/outs/filtered_feature_bc_matrix/palantir_dir
source activate py38
import palantir
import scanpy as sc
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
pca_projections=pd.read_csv('MxTA.pca.csv',index_col=0)
umap=pd.read_csv('MxTA.umap.csv',index_col=0)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=8)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
ms_data
terminal_states = pd.Series(['C1_Cortex','C6_IRS','C8_Cuticle','C9_Medulla'], index=['CATGACAAGCTAGGCA-4', 'ACACCCTAGAATGTTG-2', 'TGCGTGGCACATGTGT-4','TTGGCAATCTAACCGA-4'])
start_cell = 'TTAATCCAGGTAGATT-1'
pr_res = palantir.core.run_palantir(ms_data, start_cell, terminal_states=terminal_states.index)
palantir.plot.plot_palantir_results(pr_res, umap)
plt.show()
pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]
plt.figure(figsize=(4,12))
fig1=palantir.plot.plot_palantir_results(pr_res, umap)

###############################################################
###########integrate skin data from Liuguanghui################
###############################################################
#/local/data/wusijie/SingleCellProject/Epi_LiuGuanghui
#integrating the skin data of LiuGuanghui 
DF.list=list()
for(i in 1:9){
  print(i)
  DF.list[[i]]=Read10X_h5(f[i], use.names = TRUE, unique.features = TRUE)
  DF.list[[i]] <- CreateSeuratObject(counts = DF.list[[i]], project = "Skin", min.cells = 10, min.features = 200)
  DF.list[[i]]$group=g[i]
  DF.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object = DF.list[[i]], pattern = "^MT-")
  DF.list[[i]] <- subset(x = DF.list[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
  DF.list[[i]] <- NormalizeData(DF.list[[i]], verbose = FALSE)
  DF.list[[i]] <- SCTransform(DF.list[[i]])
}
###Fast integration Skin and HF using reciprocal PCA (RPCA)#### 
DF.list <- lapply(X = DF.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1200)
})
features <- SelectIntegrationFeatures(object.list = DF.list,nfeatures = 1200)
DF.list <- lapply(X = DF.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
DF.anchors <- FindIntegrationAnchors(object.list = DF.list, anchor.features = features, reduction = "rpca")
DF.combined <- IntegrateData(anchorset = DF.anchors)
DefaultAssay(DF.combined) <- "integrated"
DF.combined <- ScaleData(DF.combined, verbose = FALSE)
DF.combined <- RunPCA(DF.combined, npcs = 50, verbose = FALSE)
npc=25
DF.combined <- RunUMAP(DF.combined, dims = 1:npc,return.model=FALSE)
DF.combined <- FindNeighbors(DF.combined, reduction = "pca", dims = 1:npc)
DF.combined <- FindClusters(DF.combined, resolution = 0.25)
DimPlot(object = DF.combined, reduction = "umap",label = TRUE,label.size=5,pt.size = 0.2,group.by="project")

#########################################################################################
#########################calculate the PCT (expression rate)#############################
#########################################################################################
DF.list <- SplitObject(DF, split.by = "group")
PercentAbove <- function(x, threshold) {return(length(x = x[x > 0]) / length(x = x))}
AVG=list()
PCT=list()
for(i in 1:length(DF.list)){
  print(i)
  #data.features=data.frame(t(as.matrix(DF.list[[i]]@assays$regulon@data)))
  data.features=data.frame(t(as.matrix(DF.list[[i]]@assays$SCT@data)))
  data.features$id=DF.list[[i]]@active.ident
  print(dim(data.features))
  data.plot <- lapply(
    X = levels(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          #return(mean(x = expm1(x = x)))
          return(mean(x = (x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  AVG[[i]]<-PCT[[i]]<-matrix(0,ncol(data.features)-1,length(unique(x = data.features$id)))
  for(j in 1:length(levels(DF.list[[i]]@active.ident))){
    AVG[[i]][,j]=data.plot[[j]][[1]]
    PCT[[i]][,j]=data.plot[[j]][[2]]
  }
}
avgPCT=(PCT[[1]]+PCT[[2]]+PCT[[3]]+PCT[[4]]+PCT[[5]]+PCT[[6]])/6
rownames(avgPCT)<-rownames(PCT[[1]])<-rownames(PCT[[2]])<-rownames(PCT[[3]])<-rownames(PCT[[4]])<-rownames(PCT[[5]])<-rownames(PCT[[6]])<-rownames(HF6)
colnames(avgPCT)=levels(HF6@active.ident)
avgAVG=(AVG[[1]]+AVG[[2]]+AVG[[3]]+AVG[[4]]+AVG[[5]]+AVG[[6]])/6
rownames(avgAVG)<-rownames(AVG[[1]])<-rownames(AVG[[2]])<-rownames(AVG[[3]])<-rownames(AVG[[4]])<-rownames(AVG[[5]])<-rownames(AVG[[6]])<-rownames(HF6)
colnames(avgAVG)=levels(HF6@active.ident)

