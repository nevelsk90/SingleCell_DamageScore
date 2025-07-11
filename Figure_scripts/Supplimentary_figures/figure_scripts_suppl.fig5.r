###=====# Single-Cell Resolution of Cellular Damage Illuminates Disease Progression #=======###
###=====# supplimentary figure 5 PDS code #=======###

library( Seurat )
library( ggplot2 )
require( cowplot )

#### load podocyte genes
allPodoGenes <-  readRDS("PROJECTS/PDS/SCSN_allPodoGenes.rda")


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