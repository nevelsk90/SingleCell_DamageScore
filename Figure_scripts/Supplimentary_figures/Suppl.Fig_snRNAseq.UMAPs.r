

#### load podocyte genes
allPodoGenes <-  readRDS("PROJECTS/PDS/SCSN_allPodoGenes.rda")

#### analyse pdss2 ####

pdss2 <- readRDS ( "/media/tpadvits/WD_tim/PROJECTS/PODOCYTES/RNAseq/Cem_data/Seurat/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt.rda" )
pdss2 <- FindClusters( pdss2, resolution = 0.05 )
pdss2 <- FindClusters( pdss2, resolution = 0.2 )

### prepare markers
{
  # ### identify cluster markers
  pdss2_Marks.r0.05 <-  FindAllMarkers( pdss2, min.pct = 0,logfc.threshold = 0,return.thresh = 1,
                  max.cells.per.ident = 2000,features = allPodoGenes)
  # saveRDS( pdss2_Marks, "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt_Marks.r0.1.rda")
  # saveRDS( pdss2_Marks.r0.05, "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt_Marks.r0.05.rda")
  # saveRDS( pdss2_Marks, "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt_Marks.r0.1.rda")
   saveRDS( pdss2_Marks.r0.05, "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt_Marks.r0.05.rda")
  
  ###  prepare cellßtype markers markers
  marks.hmphr.sn <- read.table(header = T, sep = ";","ANNOTATIOS/Mouse/Humphreys2020_mouse_snRNA.csv")
  # marks.hmphr.sc <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2022_mouse_scRNA.csv")
  marks.sztk.sc <- readxl::read_xlsx(sheet = 1, skip = 1,  "ANNOTATIOS/Mouse/Susztak2021_mouse.scRNA.xlsx")
  marks.sztk.sc <- lapply( as.list(marks.sztk.sc) , function(X) X <- X[!is.na(X)] )
  names(marks.sztk.sc) <- sub( " ","_" , names(marks.sztk.sc))
  # convert long df to a list
  marks.hmphr.sn <- split(marks.hmphr.sn$gene, marks.hmphr.sn$celltype)
  # marks.hmphr.sc <- split(marks.hmphr.sc$gene, marks.hmphr.sc$celltype)
  podoMarks_db <- c(marks.hmphr.sn,  marks.sztk.sc)
  podoMarks_db <- lapply(podoMarks_db, function(X) X[X %in% allPodoGenes])
  # podoMarks_db <- lapply( seq(podoMarks_db), function(ii) podoMarks_db[[ii]][1:100])
  names(podoMarks_db) <- c( paste(names(marks.hmphr.sn),"hmphr.sn",sep = "_"),
                            paste(names(marks.sztk.sc),"sztk.sc",sep = "_")
  )
  
  
  
}


### do gsea to identify enrichment of cluster markers in ctype markers
{
  pdss2_Marks.r0.05 <- readRDS( file = "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt_Marks.r0.05.rda" )
  pdss2_Marks.r0.1 <- readRDS( file = "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt_Marks.r0.1.rda" )
  
  library(fgsea)
  markers <- pdss2_Marks.r0.1 
  pdss2_clustCtype.gsea.r0.05 <- lapply( levels( markers$cluster) ,
                                   function(clID)
                                   {
                                     print(clID)
                                     clXmark_LFC <- markers[
                                       markers$cluster==clID ,]
                                     
                                     clXmark_LFC <- setNames(clXmark_LFC$avg_log2FC , nm = clXmark_LFC$gene)
                                     clXmark_LFC <- clXmark_LFC[ order(-clXmark_LFC)]
                                     X <- fgseaMultilevel( stats= clXmark_LFC , scoreType="pos",
                                                           pathways = podoMarks_db ,nPermSimple = 10000,
                                                           nproc=4 , eps = 1e-100)
                                     return(X)
                                   } )
  
  # select top3
  gsea_res <-  pdss2_clustCtype.gsea.r0.1
  pdss2_clustCtype.gsea_top2 <- Reduce( union,  lapply( seq( gsea_res ), function(ii){
    print(ii)
    X <- gsea_res[[ ii ]]
    XX <- dplyr::slice_max(X, n = 3, order_by = NES )
    XX <- XX$pathway[XX$padj<0.05]
    print(XX)
    return(XX)
  }) )
  
  ### plot
  pdss2_clustCtype.gsea_p <- Reduce( cbind, lapply( seq(gsea_res), function(ii)
  {
    X <- gsea_res[[ii]]
    X$logP <- scale( -log10(X$pval), center = F)
    
    X <- X[  match( pdss2_clustCtype.gsea_top2,  X$pathway),]
    X$logP
    # X$clust <- levels(sce.decontX_subsMarks$cluster)[ii]
    # return(X)
  }))
  colnames(pdss2_clustCtype.gsea_p) <- levels(markers$cluster)
  rownames(pdss2_clustCtype.gsea_p)<- pdss2_clustCtype.gsea_top2
  pdss2_clustCtype.gsea_p[is.na(pdss2_clustCtype.gsea_p)]<-0
  
  toPlot <- t(pdss2_clustCtype.gsea_p)
  gg0 <-  pheatmap::pheatmap(toPlot,
                             color=viridis::viridis(9),
                             cluster_rows = F,
                             fontsize = 16)
  
  ### assign cell-types
  ctypes <- setNames( c( "PT.S1/2", "CNT", "EC","Podo", "PT.S2/3", "Per",
                         "TAL", "Immune", "Podo", "PT", "IC", "Pod?PT", "PT" ),
                      nm = levels( pdss2$seurat_clusters ) )
  pdss2@meta.data$ctype <- ctypes[ match( pdss2@meta.data$seurat_clusters ,
                                              names(ctypes))]
  ## add metadata
  pdss2$nCount.log_RNA <- log( pdss2$nCount_RNA +1 )
  pdss2$percent_mito <- PercentageFeatureSet(object = pdss2, pattern = "^mt-")
  pdss2$percent_ribo <- PercentageFeatureSet(object = pdss2, pattern = "^mt-")
  
  saveRDS( pdss2 , "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt.rda")
  
}


#### plot snRNAseq ####
library(Seurat)
### load data
nphs2 <- readRDS ( "PROJECTS/PDS/WT1hetdel_decontX.allcells_Seur.4w12w.scDblFilt.rda")

wt1 <- readRDS ("PROJECTS/PDS/Nphs2mut_decontX.allcells_Seur.scDblFilt.rda" )
wt1_no4w <- subset( wt1 , subset= age!= 4 )

pdss2 <- readRDS ( "PROJECTS/PDS/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt.rda" )
pdss2_noCoQ2 <- subset( pdss2 , subset= group!= "CoQ2_6")

### UMAPs
gg01 <- DimPlot(nphs2, label = T, group.by = "ctype", shuffle = T)
gg02 <-DimPlot(nphs2, group.by = "group", shuffle = T)

gg11 <- DimPlot(wt1_no4w, label = T,  group.by = "ctype", shuffle = T)
gg12 <-DimPlot(wt1_no4w, group.by = "group", shuffle = T)

gg21 <- DimPlot(pdss2_noCoQ2, label = T, group.by = "ctype", shuffle = T)
gg22 <-DimPlot(pdss2_noCoQ2, group.by = "group", shuffle = T)

ggl<-  cowplot::plot_grid(plotlist = list(gg01,gg02,
                                          gg11,gg12,
                                          gg21, gg22 ), ncol = 2)

pdf( width = 12 ,height = 12, file = "PROJECTS/PDS/KFO_UMAPs_mnscrpt.pdf")
ggl
dev.off()

png( width = 1000 ,height = 1000, file = "PROJECTS/PDS/KFO_UMAPs_mnscrpt.png")
ggl
dev.off()

### Vlnplot
vl01 <- VlnPlot(nphs2,features ="percent_mito" , group.by = "sample", pt.size = 0 , y.max = 1)
vl02 <-VlnPlot(nphs2,features ="nCount.log_RNA" , group.by = "sample", pt.size = 0 )

vl11 <- VlnPlot(wt1_no4w,features ="percent_mito" , group.by = "sample", pt.size = 0 , y.max = 1)
vl12 <-VlnPlot(wt1_no4w,features ="nCount.log_RNA" , group.by = "sample", pt.size = 0 )

vl21 <- VlnPlot(pdss2_noCoQ2,features ="percent_mito" , group.by = "sample", pt.size = 0 , y.max = 1)
vl22 <-VlnPlot(pdss2_noCoQ2,features ="nCount.log_RNA" , group.by = "sample", pt.size = 0 )

ggl2<-  cowplot::plot_grid(plotlist = list(vl01,vl02,
                                           vl11,vl12,
                                           vl21, vl22 ), ncol = 2)


pdf( width = 12 ,height = 10, file = "PROJECTS/PDS/KFO.snRNAseq_QC.vln_mnscrpt.pdf")
ggl2
dev.off()

png( width = 1000 ,height = 800, file = "PROJECTS/PDS/KFO.snRNAseq_QC.vln_mnscrpt.png")
ggl2
dev.off()

