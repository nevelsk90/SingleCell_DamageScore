# ###################################################### #
### circadian gene analysis ###
# release memory
mallinfo::mallinfo()
mallinfo::malloc.trim()
gc()

setwd("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian")

.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))
options( connectionObserver = NULL )
library( MetaCycle)
library( org.Mm.eg.db )
library( BSgenome.Mmusculus.UCSC.mm10 )
library( ggplot2 )
library( cowplot)
library( reshape2 )
library( plyr )
library( Seurat )
library( viridis )
library( ggthemes )
library( AUCell )
library( ggpubr )
library( biomaRt )
library( GSEABase )
library( RColorBrewer )

# mart_homo <- useMart( "ensembl",dataset="hsapiens_gene_ensembl" , host="www.ensembl.org")
mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl" )
tx2gene <- getBM( attributes=c( 'ensembl_gene_id', 'external_gene_name',"entrezgene_id"),  mart = mart_mouse)

# tx2prot <-  getBM( attributes=c( "uniprotswissprot","uniprot_gn_id", 'external_gene_name',"mgi_symbol"),  mart = mart_mouse)
#### load functions
# generate damage signature
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
# calculate damage score
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

### TF of interest 
TFpodo.Nr <- c(  "Rora", "Nr2f2",
                 "Rorc", "Nr1d1", "Thra", "Nr1h2", "Esrra")
cc.core <- c("Rora", "Nr1d1", "Arntl", "Clock", "Tef","Dbp","Npas2","Arntl2",
             "Cry1", "Cry2", "Per1", "Per2")

ccTFs <- c( "Rora", "Nr1d1", "Arntl", "Clock", "Tef","Dbp","Npas2","Arntl2",
            "Rorc","Nr1d2" )

# Nfil3 - inhibitory molecule

#### load data #### 
countMat <- read.table( row.names = 1, header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circadian_counts.tsv")
countMat <- countMat[,-1]

### load core circadian genes
core_clock_genes <-read.table( header = F,  "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/tempo/core_clock_genes.txt")$V1
core_clock_genes <- c( core_clock_genes , "Bmal1","Bmal2", "Bhlhe40")

### load 
# listSCSN <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_22.12.23.rda")
listSCSN.1K <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_1K.29.09.23.rda")
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")


# load GRN
ATACseq_tgenes <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/ATACseq/atac.podo_tobias.fimo77TF.p5e4_TFtc.rda")
ATACseq_tgenes.S <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "score" ))
ATACseq_tgenes.Q <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "qvalue" ))
colnames(ATACseq_tgenes.S) <- colnames(ATACseq_tgenes.Q) <- names(ATACseq_tgenes)
ATACseq_tgenesM <- ATACseq_tgenes.S * ( ATACseq_tgenes.Q <0.1)
rownames(ATACseq_tgenesM) <- rownames(ATACseq_tgenes[[1]])
ATACseq_tgenesM <- ATACseq_tgenesM[ rowSums(ATACseq_tgenesM)>0 , colSums(ATACseq_tgenesM)>0  ]

# load podocyte gene sets
PodoPathGSet <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.05.01.24.7plus.rda")

# podocyte exprsd genes
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

# all podo genes based on bulk
podoGenesFACS <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/wtPodo_bulkRNAseq/bulkRNAseq_podo_ExInt.rda")
podoGenesFACS <- podoGenesFACS$gene[rowMeans(podoGenesFACS[,1:3])>10]

## 
DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")


#### analyse bulk circadian RNAseq ####  
# visual inspection
DESeq2::plotPCA( countMat )
# remove one sample with very few counts
countMat <- countMat[, colnames(countMat)!="SN8580185_21056_806Aligned.sortedByCoord.out.bam"]
annot <- gsub(".*_|Aligned.*", "", colnames(countMat))
annot <-data.frame(  ID= colnames(countMat) ,
                     condition= annot,
                     age=c(rep(10,16), rep(80,15)),
                     time= sub("^..","",annot))
annot$age <- as.factor( annot$age )
annot$time <- as.numeric( annot$time )

# # normalised counts
# dds_rlog.norm.y <- read.table(row.names = 1, header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/dds_rlog.norm_young.tsv")
# dds_rlog.norm.o <- read.table(row.names = 1, header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/dds_rlog.norm_old.tsv")
# dds_rlog.norm <- cbind(dds_rlog.norm.y, dds_rlog.norm.o)

### get DEseq2 normalise valuse 
{
 
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix( countData = countMat ,
                                 colData = annot ,
                                 design= ~ age )
   # 
  dds_rlog <- rlog( dds )
  dds_rlog.norm <- assay( dds_rlog )
  
  #
  saveRDS( dds_rlog.norm , file= "bulkRNAseq.Kidney.cirk_rlog.rda")
  # # write data for 2 ages in seperte tabs
  # write.table(dds_rlog.norm[,annot$age==10], col.names = NA, 
  #             file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/dds_rlog.norm_young.tsv")
  # write.table(dds_rlog.norm[,annot$age!=10], col.names = NA, 
  #             file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/dds_rlog.norm_old.tsv")
  # 
}

### cluster genes
{
  toPlot <- dds_rlog.norm[ tx2gene$ensembl_gene_id[ 
    tx2gene$external_gene_name %in% 
      c( "Bmal1","Bmal2" , core_clock_genes.xtnd) ], ]
  rownames(toPlot)<- tx2gene$external_gene_name[ 
    match( rownames(toPlot), tx2gene$ensembl_gene_id)]
  library(dtwclust)
  # Partitional
  pc <- tsclust( toPlot , "hierarchical", 
                 distance = "sbd", trace = TRUE,
                 control = hierarchical_control(method = "ward.D2"))
  plot(pc)
  
}

### find cycling genes with MetaCycle
{
  MetaCycle_young <- meta2d(infile = "dds_rlog.norm_young.tsv", cycMethod =  "JTK",
                            filestyle="txt",outputFile = F,
                            timepoints	= as.numeric( annot$time[annot$age==10]))
  MetaCycle_old <- meta2d(infile = "dds_rlog.norm_old.tsv", cycMethod =  "JTK",
                          filestyle="txt",outputFile = F,
                          timepoints	= as.numeric( annot$time[annot$age!=10]))
  
  # saveRDS( MetaCycle_young , "MetaCycle_bulk.kidney_young.rda")
  # saveRDS( MetaCycle_old ,  "MetaCycle_bulk.kidney_old.rda")
  
  ### extract results
  MetaCycle_young <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/MetaCycle_bulk.kidney_young.rda")
  MetaCycle_old <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/MetaCycle_bulk.kidney_old.rda")
  
  datt.y <- MetaCycle_young$meta
  datt.o <- MetaCycle_old$meta
  
  library(org.Mm.eg.db)
  gene.symbols <- select( org.Mm.eg.db, keys = union(datt.y$CycID, datt.o$CycID) , keytype = 'ENSEMBL', columns = 'SYMBOL')
  
  datt.y$gName <- gene.symbols$SYMBOL[ match( datt.y$CycID, gene.symbols$ENSEMBL)]
  datt.y.podo <- datt.y[ datt.y$gName %in% allPodoGenes , ]
  
  datt.o$gName <-gene.symbols$SYMBOL[ match( datt.o$CycID, gene.symbols$ENSEMBL)]
  datt.o.podo <- datt.o[ datt.o$gName %in% allPodoGenes , ]
  # 
  # TFpodo.Nr.circ <- TFpodo.Nr[TFpodo.Nr %in% datt.y$gName[datt.y$JTK_BH.Q<0.1] | 
  #                               TFpodo.Nr %in% datt.o$gName[datt.o$JTK_BH.Q<0.1]]
  
  
  # datt.sig <- datt[datt$JTK_BH.Q<0.01,]
  # datt.sig <- datt.sig[order(datt.sig$JTK_BH.Q),]
  
}

### GLMMcosinor, calculate acrophase
{
  library(GLMMcosinor)
  
  dds_rlog.norm <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian")
  
#   # all genes identified in podocytes
#    datt <- cbind.data.frame( Time=annot$time, Age= annot$age , 
#                  t(dds_rlog.norm[ rownames( dds_rlog.norm) %in%
#       tx2gene$ensembl_gene_id[ tx2gene$external_gene_name %in% 
#                                 c( "Bmal1","Bmal2", 
#                                    union( allPodoGenes, podoGenesFACS) )
#                                ] , ]) )
#   
#   # estimate acrophase uncertainty, run as a background
#   GLMMcosinor.all <- lapply( 3:ncol( datt), function(ii)
#   {
#     print(ii)
#     dattTest <- datt[, c(1:2,ii)]
#     colnames(dattTest) <- c("Time","Age","Y")
#     GLMMfit <- tryCatch( cglmm(
#       Y ~ Age + amp_acro( time_col=Time, period = 24, group = "Age"),
#       data = dattTest ) , error = function(e) NA)
#     
#     if(!is.na(GLMMfit)) {
#       XX <- summary(GLMMfit)$transformed.table[
#         grep(".*acr.*", ignore.case = T, rownames( summary(GLMMfit)$transformed.table)),] 
#     } else XX <- NA
#     return(XX)
#   } )
#   names(GLMMcosinor.all) <- tx2gene$external_gene_name[
#     match( colnames(datt)[-c(1:2)],tx2gene$ensembl_gene_id)]
# # save
#   saveRDS(GLMMcosinor.all, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/bulkKidney.all_GLMMfit.rda")
#   
  ### top scsn and core clock genes
  # select top scsn genes
  scsnTop0.95 <- Reduce( intersect , lapply( listSCSN , function(XX){
    rownames(XX@assays$RNA@data)[ rowMeans(XX@assays$RNA@data) > quantile( rowMeans(XX@assays$RNA@data), 0.95 )]
  }))
  
  ddd <-   t(dds_rlog.norm[ rownames( dds_rlog.norm) %in%
                                    tx2gene$ensembl_gene_id[ 
                                      tx2gene$external_gene_name %in% 
                                        core_clock_genes
                                        # union( core_clock_genes , scsnTop0.95 )
                                      ] , ]) 
  # ddd <- scale( ddd )
  datt.circ.core <- cbind.data.frame( Time=annot$time, Age= annot$age , ddd )
  
  Z <- rowSums( ddd[,c("ENSMUSG00000040998","ENSMUSG00000055116")])

  # estimate acrophase uncertainty, run as a background
  GLMMcosinor.circ <- lapply( 3:ncol( datt.circ.core), function(ii)
    {
    print(ii)
    dattTest <- datt.circ.core[, c(1,2, ii)]
    colnames(dattTest) <- c("Time","Age","Y")
    GLMMfit <- tryCatch( cglmm(
      Y ~ amp_acro( time_col=Time, period = 24 ),
      data = dattTest ) , error = function(e) NA)
    
    if(!is.na(GLMMfit)) {
      XX <- summary(GLMMfit)$transformed.table[
        grep(".*acr.*", ignore.case = T,
             rownames( summary(GLMMfit)$transformed.table)),]
    } else XX <- NA
    return(list(GLMMfit, XX))
  } )
  
  
  GLMMcosinor.circ.model <- lapply( GLMMcosinor.circ, "[[",1)
  GLMMcosinor.circ.Tab <- lapply( GLMMcosinor.circ, "[[",2)
  names(GLMMcosinor.circ) <- names(GLMMcosinor.circ.Tab) <-  
    names(GLMMcosinor.circ.model) <- tx2gene$external_gene_name[
      match( colnames( ddd) ,tx2gene$ensembl_gene_id)]
  GLMMcosinor.circ.Tab <- Reduce( rbind, GLMMcosinor.circ.Tab)
  GLMMcosinor.circ.Tab$gName <- names(GLMMcosinor.circ)
  GLMMcosinor.circ.Tab$gName <- plyr::revalue(GLMMcosinor.circ.Tab$gName , 
                                              c("Bmal1" = "Arntl", "Bmal2" = "Arntl2"))
  
  # plot
  autoplot(GLMMcosinor.circ.model$Bmal1, superimpose.data = TRUE)
  polar_plot(GLMMcosinor.circ.model$Bmal1)
  autoplot(GLMMcosinor.circ.model$Bmal2, superimpose.data = TRUE)
  polar_plot(GLMMcosinor.circ.model$Bmal2)
  autoplot(GLMMcosinor.circ.model$Clock, superimpose.data = TRUE)
  polar_plot(GLMMcosinor.circ.model$Clock)
 
  # save 
  saveRDS( GLMMcosinor.circ , file="GLMMcosinor__scsnTop0.95_AND_core.clock.rda")
  saveRDS( GLMMcosinor.circ.Tab , file="bulkKidney_GLMMcosinor.AcrPrior__scsnTop0.95.AND.core.clock.rda")
  
  
  # plot coefficients
  GLMMcosinor.circ.Tab$CI.95<- GLMMcosinor.circ.Tab$standard.error*2*1.96
  toPlot <- GLMMcosinor.circ.Tab[ ( GLMMcosinor.circ.Tab$CI.95 < pi & 
                                      GLMMcosinor.circ.Tab$gName %in% core_clock_genes )| 
                                    GLMMcosinor.circ.Tab$CI.95 < 0.5  ,]
  toPlot$estimate <-ifelse( toPlot$estimate<0, 2*pi+toPlot$estimate, toPlot$estimate)
  toPlot$lower.CI <-ifelse( toPlot$lower.CI<0, 2*pi+toPlot$lower.CI, toPlot$lower.CI)
  toPlot$upper.CI <-ifelse( toPlot$upper.CI<0, 2*pi+toPlot$upper.CI, toPlot$upper.CI)
  
  # toPlot <- toPlot[! toPlot$gName %in% c("Arntl2", "Csnk1a1"),]
  gg<- ggplot(data= toPlot , aes(x= reorder(gName, estimate)))+ 
    geom_col( aes( y= estimate) ,position = position_dodge()) +
    theme_bw()+ scale_fill_colorblind()+
    geom_errorbar(aes(x=gName, ymin= lower.CI, ymax= upper.CI ),
                  position = position_dodge(0.9),show.legend = T)+  
    xlab("core circadian genes")+ ggtitle("acrophase estimated by GLMMcosinor\nfrom kidney bulkRNAseq, error bars show 95CI")+
    theme( text = element_text(size=18))+
    # geom_text(aes( y= estimate, label = ifelse(p.value<0.01, "*", "")), 
    #           position = position_dodge(width =  .9), size = 9 , color="red") +
    coord_flip()
  gg
  # save
  saveRDS( toPlot$gName , file="core_clock_genes.xtnd2.rda")
  
  ### plot expression of the cycling genes in snsc
  core_clock_genes.scsnMean <- Reduce( cbind.data.frame, lapply( seq(listSCSN), function(ii){
    datt <- listSCSN[[ii]]
    XX <- rowMeans( datt@assays$RNA@data[ rownames(datt@assays$RNA@data) %in% toPlot$gName , ])
    XX <- XX[ match( toPlot$gName , names(XX))]
  }) )
  colnames(core_clock_genes.scsnMean)<- names(listSCSN)
  core_clock_genes.scsnMean$gName= rownames(core_clock_genes.scsnMean)
  
  toPlot <- melt(core_clock_genes.scsnMean)
  gg2 <- ggplot( data=toPlot , 
          aes(y=reorder( gName, value),x=value, color=variable))+
    geom_point( alpha=0.75, size=5)+theme_bw()+theme( text=element_text( size=18))+
    ylab("average Expr.lvl")+ggtitle("snRNAseq Nphs2mut. podocytes")+
    scale_color_tableau()
  
  # combine plots 
  ggl <- cowplot::plot_grid( plotlist = list( gg, gg2),align = T)
  
  pdf( width = 16, height = 8, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/tempo/GLMMcosinor.acrophase_core.clock_barplot.pdf")
    print( ggl)
  dev.off()
   
   

}

### plot cycling gene curves
{  
  dds_rlog.norm <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/bulkRNAseq.Kidney.cirk_rlog.rda")

  getAMP <- function(expr, per, pha, tim= as.numeric(annot$time)[1:16])
  { 
    trendt <- tim - mean(tim[!is.na(tim) & !is.nan(tim)])
    cost <- cos(2*pi/per*(tim - pha))
    fit <- lm(expr~trendt + cost)
    fitcoef <- fit$coefficients
    basev <- fitcoef[1]
    trendv <- fitcoef[2]
    ampv <- fitcoef[3]
    fitexp <- basev + trendv*trendt + ampv*cost
    outL <- list("base"=basev, "trend"=trendv, "amp"=ampv, "fit"=fitexp)
    return(outL)
  }
  
  ggenes <- TFpodo.Nr
  ggenes <- c("Bmal1","Nr1d1","Cry2","Per1")
  # core circadian machinery
  # genes that sig correlate with PDS in many studies and cycling
  # ggenes <- datt.y.podo$gName[ datt.y.podo$gName %in% 
  #                                Reduce( union, PodoPathGSet[1:4]) &
  #                                (datt.y.podo$JTK_BH.Q < 0.1 |  datt.o.podo$JTK_BH.Q < 0.1) ]
  
  pdf(width = 10, height = 12,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/XX.pdf")
  par( mfrow= c(2, 2))
  
  lapply( seq(ggenes), function(ii)
  {
    ggeneID <- ggenes[[ii]]
    ggeneID <- tx2gene$ensembl_gene_id[ tx2gene$external_gene_name == ggeneID ]
    # ampL.y <- getAMP( expr=as.numeric(dds_rlog.norm[ggeneID,annot$age==10]), 
    #                   per=datt.y[datt.y$CycID==ggeneID, "JTK_period"], 
    #                   pha=datt.y[datt.y$CycID==ggeneID, "JTK_adjphase"],
    #                   tim=annot$time[annot$age==10] )
    # ampL.o <- getAMP( expr=as.numeric(dds_rlog.norm[ggeneID,annot$age==80]), 
    #                   per=datt.o[datt.o$CycID==ggeneID, "JTK_period"], 
    #                   pha=datt.o[datt.o$CycID==ggeneID, "JTK_adjphase"],
    #                   tim=annot$time[annot$age==80] )
    # lay<-layout(cbind(1, 2), widths=c( lcm(cm(4.5)), lcm(cm(1.5)) ), heights=lcm(cm(4.5)) )
    # par(mai=c(0.65,0.6,0.4,0.05),mgp=c(2,0.5,0),tck=-0.01)
    # xrange <- c(18, 65)
    # yrange <- c(200, 2350)
    # 
    # plot(annot$time[annot$age==10] , as.numeric(dds_rlog.norm[ ggeneID ,annot$age==10]), 
    #      type="b", xlab="Circadian time(CT)", 
    #      ylab="Expression value", main=datt.sig$gName[
    #        datt.sig$CycID==ggeneID], cex.main=1.2 , col="blue",)
    # par(new=T)
    # 
    # plot(annot$time[annot$age==80] , as.numeric(dds_rlog.norm[ ggeneID ,annot$age==80]), 
    #      type="b", xlab="Circadian time(CT)", 
    #      ylab="Expression value", main=datt.sig$gName[
    #        datt.sig$CycID==ggeneID], cex.main=1.2 , col="red",)
    
    
    loessD.y <- data.frame( expd= as.numeric(dds_rlog.norm[ ggeneID ,annot$age==10]),
                            tp=annot$time[annot$age==10] )
    # loessD.y$expd <- loessD.y$expd/max( loessD.y$expd)
    exploess.y <- loess( expd~tp, loessD.y, span = 0.4)
    expsmooth.y <- predict(exploess.y, data.frame(tp=annot$time[annot$age==10]))
    
    loessD.o <- data.frame(expd= as.numeric(dds_rlog.norm[ ggeneID ,annot$age==80]),
                           tp=annot$time[annot$age==80] )
    # loessD.o$expd <-loessD.o$expd/max( loessD.o$expd)
    exploess.o <- loess( expd~tp, loessD.o, span = 0.4 )
    expsmooth.o <- predict( exploess.o, data.frame(tp=annot$time[annot$age==80]))
    
    plot( annot$time[ annot$age == 10 ] , 
          ylim=c(min(expsmooth.y,expsmooth.o), 
                 max(expsmooth.y,expsmooth.o)),
          xlab="circadian time", ylab="rlogNorm expr.lvl",
          expsmooth.y, main= paste0(ggenes[[ii]],"\n","young qval = ",
                                    round( datt.y$JTK_BH.Q[ datt.y$CycID==ggeneID ], 3),
                                    "\n","old qval = ",
                                    round( datt.o$JTK_BH.Q[ datt.o$CycID==ggeneID ], 3)),
          lwd=2 , col = "blue" , type = "l" )
    lines( annot$time[annot$age==80] , expsmooth.o , 
           lwd=2, col="red" , type="l" )
    
    # # par(new=T)
    # plot( annot$time[ annot$age==10], ampL.y[[4]],
    #       type="b", col="lightblue", xlab="", ylab="", main="" )
    # par(new=T)
    # plot( annot$time[ annot$age==80], ampL.o[[4]],
    #       type="b", col="lightsalmon", xlab="", ylab="", main="" )
  })
  
  
  dev.off()
  
  
}

### plot sample PCA for cycling genes only
{
  
  ### do PCA
  library(factoextra)
  library(viridis)
  
  
  PCAdat_sel <- dds_rlog.norm
  PCAdat_sel <- PCAdat_sel[ rownames(PCAdat_sel) %in% union(
    datt.o.podo$CycID[datt.o.podo$JTK_BH.Q<0.05] ,  datt.y.podo$CycID[datt.y.podo$JTK_BH.Q<0.05]
  ),]
  colnames(PCAdat_sel) <- paste0( annot$age ,"_", annot$time,  "hrs")
  ## run PCA  
  res.pca <- prcomp( t(PCAdat_sel), scale = F)
  
  # plot results
  fviz_eig(res.pca)
  biplot(res.pca, choices = )
  gg <-  factoextra::fviz_pca_ind(res.pca, 
                           habillage = ifelse( annot$age==10, "young","old" ))+ 
    scale_color_colorblind()
  
  
  
}

### analyse cycling of PDS, MDS and AGS
{
  library(GLMMcosinor)
  
  
  # ddd <- scale( ddd )
  datt.circ.core <- cbind.data.frame( Time=annot$time, Age= annot$age , 
                                      PDS=PDS ,MDS.all=MDS.all , AGS.all = AGS.all)
  # datt.circ.core_melt <- reshape2::melt(data =   datt.circ.core, id.var = c("Time","Age"))

  # estimate acrophase uncertainty, run as a background
    dattTest <- datt.circ.core[, c(1,2, 3)]
    colnames(dattTest) <- c("Time","Age","PDS")
    GLMMfit.PDS <- tryCatch( cglmm(
      PDS ~ amp_acro( time_col=Time, period = 24, group = "Age" ),
      data = dattTest ) , error = function(e) NA)
    
    dattTest<- datt.circ.core[, c(1,2, 4)]
    colnames(dattTest) <- c("Time","Age","PDS")
    GLMMfit.MDS <- tryCatch( cglmm(
      PDS ~ amp_acro( time_col=Time, period = 24, group = "Age" ),
      data = dattTest ) , error = function(e) NA)
    
    dattTest<- datt.circ.core[, c(1,2, 5)]
    colnames(dattTest) <- c("Time","Age","PDS")
    GLMMfit.AGS <- tryCatch( cglmm(
      PDS ~ amp_acro( time_col=Time, period = 24, group = "Age" ),
      data = dattTest ) , error = function(e) NA)
    
    autoplot(GLMMfit.PDS)
    autoplot(GLMMfit.MDS)
    autoplot(GLMMfit.AGS)
    
 
}

#### check overrepresentation of cycling genes in other gsets ####

### significant overrrepresentation in podo paths
{
  thrsh <- 0.1
  datt <- Reduce( rbind.data.frame, lapply( seq( PodoPathGSet ) , 
                                            function(jj)  {
                                              print( jj )
                                              ggnames <- PodoPathGSet[[jj]]
                                              
                                              # datt<- sapply( seq(mmotifs), function(ii){
                                              #   grn  <- rownames( ATACseq_tgenesM)[ ATACseq_tgenesM[ ,mmotifs[ii]] >0]
                                              # x1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% intersect( grn, ggnames) & 
                                              #                                 datt.y$gName %in% allPodoGenes ] < thrsh ) 
                                              x1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% ggnames &
                                                                              datt.y$gName %in% allPodoGenes ] < thrsh )
                                              x2 <- table( datt.o$JTK_BH.Q[ datt.o$gName %in% ggnames & 
                                                                              datt.o$gName %in% allPodoGenes ] < thrsh )
                                              X1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% allPodoGenes ] < thrsh )
                                              X2 <- table( datt.o$JTK_BH.Q[ datt.o$gName %in% allPodoGenes ] < thrsh )
                                              p.y <- phyper( x1[2]-1,  X1[2], X1[1],sum(x1), lower.tail = F, log.p = FALSE)
                                              p.o <-phyper( x2[2]-1,X2[2], X2[1],sum(x2), lower.tail = F, log.p = FALSE)
                                              
                                              # print( paste0(
                                              #   # "hyper geometric test for enrichment of ",
                                              #   mmotifs[ii],
                                              #   " tgenes in circadian genes, ", "young pval=", 
                                              #   round( p.y,3) , ", old pval=", round( p.o,3)))
                                              # return(c(p.y,p.o))
                                              
                                              # colnames(datt) <- mmotifs
                                              # rownames(datt) <-  c( paste0("young_" ,plotrix::color.id(levels(cclust)[jj]) ),
                                              #                      paste0("old_" , plotrix::color.id(levels(cclust)[jj])))
                                              return(c(p.y,p.o))
                                            } ) )
  
  colnames(datt) <- c( "young_phyper" , "old_phyper")
  rownames(datt) <-  names(PodoPathGSet)
  
  datt$young_phyper_padj <- p.adjust(datt$young_phyper, method = "fdr")
  datt$old_phyper_padj <- p.adjust(datt$old_phyper, method = "fdr")
  siglbl <- round( datt[,3:4], 3)
  # siglbl <- ifelse( datt[,3:4] < 0.1,"","X")
  # siglbl <- siglbl[ rowSums(!is.na(toplot))>0, ]
  
  pheatmap::pheatmap( -log10(datt[,c(3,4)]), color =(viridis::viridis(12)), 
                      display_numbers = siglbl , fontsize_number = 14 ,
                      number_color = "red") 
  
  ### check overrepresenttion of circadian in DS
  thrsh <- 0.05
  ggnames <- DS_all$gene_symbol[1:42]
  x1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% ggnames & 
                                  datt.y$gName %in% allPodoGenes ] < thrsh )
  X1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% allPodoGenes ] < thrsh )
  phyper( x1[2]-1,  X1[2], X1[1],sum(x1), lower.tail = F, log.p = FALSE)
  
  
}

### use GRN - check how enriched are predicted targets in cycling genes
thrsh <- 0.01
GRNatac <- ATACseq_tgenesM.podo

GRNatac_circ.test <- Reduce( cbind, lapply( seq( from=0.01 , to=0.21,  0.025 ), 
                                            function(thrsh){
                                              print(thrsh)
                                              
                                              GRNatac<- ATACseq_tgenesM.podo
                                              
                                              circ.tgenes.phyper <- sapply( seq(colnames(GRNatac)), function(ii)
                                              {
                                                
                                                
                                                datt <- GRNatac[ , ii ]
                                                datt <- datt[datt>0 & names(datt)%in% allPodoGenes]
                                                x1 <- table( datt.y.podo$JTK_BH.Q[ datt.y.podo$gName %in% names(datt) ] < thrsh ) 
                                                # x2 <- table( datt.o$JTK_BH.Q[ datt.o$gName %in% names(datt) & 
                                                # datt.o$gName %in% allPodoGenes ] < thrsh )
                                                X1 <- table( datt.y.podo$JTK_BH.Q < thrsh )
                                                # X2 <- table( datt.o$JTK_BH.Q[ datt.o$gName %in% allPodoGenes ] < thrsh )
                                                p.circ <- phyper( x1[2]-1,  X1[2], X1[1],sum(x1), lower.tail = F, log.p = FALSE)
                                                # p.o <-phyper( x2[2]-1,X2[2], X2[1],sum(x2), lower.tail = F, log.p = FALSE)
                                                
                                                # print( paste0(
                                                #   # "hyper geometric test for enrichment of ",
                                                #   colnames(GRNatac)[ii],
                                                #   " tgenes in circadian genes, ", "young pval=", 
                                                #   round( p.y,3) , ", old pval=", round( p.o,3)))
                                                # return(c(p.y,p.o))
                                              })
                                              
                                              # circ.tgenes.phyper <- t(circ.tgenes.phyper )
                                              # rownames( circ.tgenes.phyper) <- colnames(ATACseq_tgenesM)
                                              # colnames( circ.tgenes.phyper )<- c("young_phyper", "old_phyper")
                                              # p.adjust( unlist(circ.tgenes.phyper), method="fdr")
                                              
                                              # datt <- as.data.frame( circ.tgenes.phyper )
                                              # datt$young_phyper_padj <- p.adjust(datt$young_phyper, method = "fdr")
                                              # datt$old_phyper_padj <- p.adjust(datt$old_phyper, method = "fdr")
                                              
                                              # return( circ.tgenes.phyper[1,])
                                            }) )


rownames( GRNatac_circ.test) <- colnames(ATACseq_tgenesM)
colnames( GRNatac_circ.test )<-seq( from=0.01 , to=0.21,  0.025 )

toPlot <- apply( GRNatac_circ.test, 2, p.adjust, method="fdr")
siglbl <- round( GRNatac_circ.test, 3)
# siglbl <- ifelse( datt[,3:4] < 0.1,"","X")
# siglbl <- siglbl[ rowSums(!is.na(toplot))>0, ]

pheatmap::pheatmap( -log10( toPlot ), color =(viridis::viridis(12)), 
                    display_numbers = siglbl , fontsize_number = 14 ,
                    number_color = "red") 


#### analyse scsnRNAseq studies study and relate to PDS #### 

nphs2.podo_seu <- subset( listSCSN$Nphs2 , subset= sample %in% c(
  "140739", "140738", "140740", "140741", "139919",
  "139917", "139921", "139913", "139915", "139911" ))
# nphs2.podo_seu <- listSCSN$Nphs2
Idents(nphs2.podo_seu) <- "sample"
set.seed(42)
nphs2.podo_seu <- subset( nphs2.podo_seu , downsample=500)
nphs2.podo_seu$sampleDiscr <- paste0( nphs2.podo_seu$sample,"_",
                                      nphs2.podo_seu$gtypeDE)

### perform DE
{
  listSCSN_gtypeDE<- lapply( seq(listSCSN), function(ii){
    podo_DE <- listSCSN[[ii]]
    if(ii %in% c(7,9)) 
      podo_DE<-  merge( podo_DE, subset( listSCSN[[8]] , subset=gtypeDE=="control"))
    Idents(podo_DE) <- "gtypeDE"
    
    podo_DE <- subset( podo_DE, downsample=1000)
    podo_DE <-  FindMarkers( podo_DE,  ident.1 = "experimental", ident.2 = "control", verbose = T)
  })
  names(listSCSN_gtypeDE) <- names(listSCSN)
  
  
  ### test over-representation of cycling genes in DE genes
  thrsh <- 0.05
  
  gtypeDE.circTest <- Reduce( rbind.data.frame, 
                              lapply( seq(listSCSN_gtypeDE),
                                      function(ii){
                                        
                                        datt <- listSCSN_gtypeDE[[ii]]
                                        ggnames <- rownames(datt)[datt$p_val_adj< thrsh]
                                        x1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% ggnames &
                                                                        datt.y$gName %in% allPodoGenes ] < thrsh )
                                        x2 <- table( datt.o$JTK_BH.Q[ datt.o$gName %in% ggnames & 
                                                                        datt.o$gName %in% allPodoGenes ] < thrsh )
                                        X1 <- table( datt.y$JTK_BH.Q[ datt.y$gName %in% allPodoGenes ] < thrsh )
                                        X2 <- table( datt.o$JTK_BH.Q[ datt.o$gName %in% allPodoGenes ] < thrsh )
                                        p.y <- phyper( x1[2]-1,  X1[2], X1[1],sum(x1), lower.tail = F, log.p = FALSE)
                                        p.o <-phyper( x2[2]-1,X2[2], X2[1],sum(x2), lower.tail = F, log.p = FALSE)
                                        
                                        print( paste0(
                                          # "hyper geometric test for enrichment of ",
                                          names(listSCSN)[ii],
                                          " tgenes in circadian genes, ", "young pval=", 
                                          round( p.y,3) , ", old pval=", round( p.o,3)))
                                        return(c(p.y,p.o))
                                        # # 
                                        # nphs2.podo_DE.circ <- rownames( nphs2.podo_DE)[ (rownames(nphs2.podo_DE) %in% c(
                                        #   datt.y$gName[ datt.y$JTK_BH.Q < thrsh], datt.o$Name[ datt.y$JTK_BH.Q < thrsh] )) &
                                        #     nphs2.podo_DE$p_val_adj<thrsh]
                                      }))
  
  colnames(gtypeDE.circTest) <- c( "young_phyper" , "old_phyper")
  rownames(gtypeDE.circTest) <- names(listSCSN)
  
  gtypeDE.circTest$young_phyper_padj <- p.adjust(gtypeDE.circTest$young_phyper, method = "fdr")
  gtypeDE.circTest$old_phyper_padj <- p.adjust(gtypeDE.circTest$old_phyper, method = "fdr")
  siglbl <- round( gtypeDE.circTest[,3:4], 3)
  # siglbl <- ifelse( datt[,3:4] < 0.1,"","X")
  # siglbl <- siglbl[ rowSums(!is.na(toplot))>0, ]
  
  pheatmap::pheatmap( -log10(gtypeDE.circTest[,3:4]), cluster_rows = F,
                      color =(viridis::viridis(12)), 
                      display_numbers = siglbl , fontsize_number = 14 ,
                      number_color = "red") 
  
}

### calculate corr for each sample
{
  listSCSN_PDScorr.smpl<- lapply( seq(listSCSN), function(jj, datt=listSCSN)
  {
    
    podo.seu <- datt[[jj]]
    Idents(podo.seu) <- "sample"
    set.seed(42)
    podo.seu <- subset( podo.seu , downsample=500)
    podo.seu$sampleDiscr <- paste0( podo.seu$sample,"_",
                                    podo.seu$gtypeDE)
    
    snames=names(table( podo.seu$sample))[table(podo.seu$sample)> 50 ]
    snamesDE=names(table( podo.seu$sampleDiscr))[table(podo.seu$sample)> 50 ]
    
    PDS.SpCor_persample <- lapply( seq(snames),
                                   function(jj){
                                     print(snames[jj])
                                     # extract normalised data to avoid spurios corr, 
                                     # due to FSGS related change the lib.size
                                     datt.smpl <- podo.seu@assays$RNA@data[
                                       , podo.seu$sample==snames[jj] ]
                                     datt.smpl <- datt.smpl[ rowSums(datt.smpl>0) > 
                                                               sqrt(ncol(datt.smpl)) , ]
                                     datt.smpl <- t( as.matrix( datt.smpl ))
                                     
                                     
                                     # datt.smpl <- datt.smpl[ datt.smpl$PDS.42>  mean(datt.smpl$PDS.42)-3*sd(datt.smpl$PDS.42) & 
                                     #                           datt.smpl$PDS.42<  mean(datt.smpl$PDS.42)+3*sd(datt.smpl$PDS.42), ]
                                     PDS.42 <- podo.seu$PDS[
                                       podo.seu$sample== snames[jj] ]
                                     
                                     cor.smpl <- psych::corr.test(   y= PDS.42 ,
                                                                     x = datt.smpl,
                                                                     method = "spearman" )
                                     
                                     
                                     # calculate q-values
                                     cor.qval.smpl <- qvalue::qvalue( cor.smpl$p )$qvalues
                                     
                                     return( cbind.data.frame( cor.r=cor.smpl$r , 
                                                               cor.p=cor.smpl$p  , 
                                                               cor.qval = cor.qval.smpl ,
                                                               sample=snames[jj] ,
                                                               gtypeDE=unique( podo.seu$gtypeDE[
                                                                 podo.seu$sample== snames[jj] ]
                                                               )) )
                                   })
    
    names(PDS.SpCor_persample) <- snames
    return(PDS.SpCor_persample)
  })
  
  listSCSN_PDScorr.smpl.indv <- Reduce( c, listSCSN_PDScorr.smpl)
  # nphs2.podo_PDS.SpCor_persample <- Reduce( rbind, nphs2.podo_PDS.SpCor_persample)
  allgenes <- Reduce( union, lapply( seq(listSCSN_PDScorr.smpl.indv) , 
                                     function(ii) rownames(
                                       listSCSN_PDScorr.smpl.indv[[ii]] )))
  
  ### extract r and p
  {
    podo_PDS.SpCor_r <- Reduce( cbind.data.frame,  lapply( seq(listSCSN_PDScorr.smpl.indv) , 
                                                           function(ii) listSCSN_PDScorr.smpl.indv[[ii]]$cor.r[
                                                             match( allgenes, rownames( listSCSN_PDScorr.smpl.indv[[ii]]))] ) )
    # centrd
    podo_PDS.SpCor_r.cntrd <- apply( podo_PDS.SpCor_r, 2, scale, scale = F)
    
    # pval
    podo_PDS.SpCor_p <- Reduce( cbind.data.frame,  lapply( seq(listSCSN_PDScorr.smpl.indv) , 
                                                           function(ii) listSCSN_PDScorr.smpl.indv[[ii]]$cor.p[
                                                             match( allgenes, rownames( listSCSN_PDScorr.smpl.indv[[ii]]))] ) )
    # qval
    podo_PDS.SpCor_q <- Reduce( cbind.data.frame,  lapply( seq(listSCSN_PDScorr.smpl.indv) , 
                                                           function(ii) listSCSN_PDScorr.smpl.indv[[ii]]$cor.qval[
                                                             match( allgenes, rownames( listSCSN_PDScorr.smpl.indv[[ii]]))] ) )
    
  }
  
  colnames( podo_PDS.SpCor_r) <- colnames( podo_PDS.SpCor_p) <- 
    colnames( podo_PDS.SpCor_r.cntrd) <- 
    colnames( podo_PDS.SpCor_q) <- names( listSCSN_PDScorr.smpl.indv )
  
  
  
  rownames( podo_PDS.SpCor_r) <- rownames( podo_PDS.SpCor_r.cntrd) <- 
    rownames( podo_PDS.SpCor_p) <-rownames( podo_PDS.SpCor_q) <- allgenes
  
}

### plot circadian genes PDScorr heatmap for individual samples 
{
  thrsh <- 0.1
  ggenes1 <- c("Rora", "Nr1d1", "Arntl", "Clock", "Cry1", "Cry2", "Per1", "Per2")
  ggenes2 <- datt.y.podo$gName[ datt.y.podo$gName %in% 
                                  Reduce( union, PodoPathGSet[1:4]) &
                                  (datt.y.podo$JTK_BH.Q < 0.1 |  datt.o.podo$JTK_BH.Q < 0.1) ]
  ggenes <- union(ggenes1, ggenes2)
  
  toplot <- podo_PDS.SpCor_r.cntrd[ 
    match( ggenes,  rownames(podo_PDS.SpCor_r.cntrd))  ,]
  # toplot[ nphs2.podo_PDS.SpCor_q[
  #  match( ggenes,  rownames(nphs2.podo_PDS.SpCor_q)),] > thrsh] <- 0
  rownames(toplot) <- ifelse( rownames(toplot)%in% DS_all$gene_symbol[1:42],
                              paste0( rownames(toplot),"_DS"),
                              rownames(toplot) )
  rownames(toplot) <- ifelse( rownames(toplot)%in% ggenes1,
                              paste0( rownames(toplot),"_CC"),
                              rownames(toplot) )
  
  toplotR <- toplot[ rowSums(!is.na(toplot))>0, ]
  siglbl <- ifelse( podo_PDS.SpCor_q[ match( 
    ggenes,  rownames( podo_PDS.SpCor_q)), ] < thrsh,"","X")
  siglbl <- siglbl[ rowSums(!is.na(toplot))>0, ]
  
  # toplot[is.na(toplot)] <- "0"
  rg <- max(abs(toplotR), na.rm = T);
  gg <- pheatmap::pheatmap( toplotR , cluster_cols = F, cluster_rows =  T, 
                            clustering_method = "ward.D2",
                            color=colorRampPalette(c("navy", "white", "red"))(9),
                            fontsize_number = 14 ,
                            display_numbers = siglbl, breaks = seq(-rg, rg, length.out = 10))
  
}


### plot circadian genes PDScorr heatmap for studies
{
  PDS42.SpCor_barcode.combo<-  read.table( sep = "\t",header = T,
                                           file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/PDS42.sprmCorr.barcode_centered.05.01.24.tsv")
  
  toplot <- genecorrPDS.Sprmn_r.cntrd[ 
    match( ggenes,  rownames(genecorrPDS.Sprmn_r.cntrd))  , 10:18]
  
  siglbl <- ifelse( genecorrPDS.Sprmn_q[ match( ggenes,  rownames( genecorrPDS.Sprmn_q)),
                                         10:18] < thrsh ,"","X")
  rownames(toplot) <- ifelse( rownames(toplot)%in% DS_all$gene_symbol[1:42],
                              paste0( rownames(toplot),"_DS"),
                              rownames(toplot) )
  rownames(toplot) <- ifelse( rownames(toplot)%in% ggenes1,
                              paste0( rownames(toplot),"_CC"),
                              rownames(toplot) )
  rg <- max(abs(toplot), na.rm = T);
  pheatmap::pheatmap( toplot, clustering_method = "ward.D2", 
                      color=colorRampPalette(c("navy", "white", "red"))(15),
                      cluster_cols = F,cluster_rows =  T, display_numbers = siglbl,
                      breaks = seq(-rg, rg, length.out = 16), fontsize_number = 14)
}



#### analyse circadian (dis)regulation CRD ####
# circ.genes.bulk0.01 <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circ.genes.bulk0.01.rda")
podo.circ <-   readxl::read_xlsx("PROJECTS/PODOCYTE/DiseaseScore/circadian/Podocyte circadian transcriptomics meta2d BHQ _0.1 .xlsx")

# podo.circ.0.05 <- podo.circ$Gene_Symbol[ 
#   podo.circ$Gene_Symbol%in% allPodoGenes & 
#     podo.circ$JTK_BH.Q < 0.05 ] 
podo.circ.0.01 <- podo.circ$Gene_Symbol[
  podo.circ$Gene_Symbol%in% allPodoGenes &
    podo.circ$JTK_BH.Q < 0.01 ]
# circ.genes.bulk0.05 <-  readRDS("circ.genes.bulk0.05.rda")
# listSCSN.1K <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_1K.22.12.23.rda")

### CRD score https://github.com/leihe2021/CRDscore/ 
  {
  library( CRDscore )

  ### retrieve gene length for TPM calculation for all genes, just in case
  gene.symbols <- select( org.Mm.eg.db, keys = rownames(countMat) ,
                          keytype = 'ENSEMBL', columns = 'SYMBOL')
  gene.symbols <- gene.symbols[!is.na(gene.symbols$SYMBOL),]

  GeneLengthAndGCContent <- EDASeq::getGeneLengthAndGCContent(
    org= "mm10" , mode="org.db" , gene.symbols$ENSEMBL )

  GeneLengthAndGCContent <- as.data.frame(GeneLengthAndGCContent)
  GeneLengthAndGCContent$gName <- gene.symbols$SYMBOL
  GeneLengthAndGCContent <- aggregate( .~ gName, data=GeneLengthAndGCContent, FUN=mean)
  GeneLengthAndGCContent <- GeneLengthAndGCContent[-1,]


  # calculate CRD score
  listSCSN.CRD <- lapply( seq(listSCSN.1K.sampl), function(ii)
      {
    print(ii)

    seu <- listSCSN.1K.sampl[[ii]]

    seu$sampleDiscr <- paste0( seu$sample,"_", seu$gtypeDE)


    # select genes with estimated gene length (for TPM calc)
    ccounts <-  seu@assays$RNA@counts
    ccounts <- ccounts[ rownames(ccounts) %in%
                          GeneLengthAndGCContent$gName[
                            !is.na(GeneLengthAndGCContent$length)
                          ], ]




    ## calculate tpm
    seu_tpm <- counts_to_tpm(  ccounts ,
                               featureLength=GeneLengthAndGCContent$length[
                                 match(rownames(ccounts), GeneLengthAndGCContent$gName)],
                               meanFragmentLength= rep( 89,ncol(ccounts) ))


    # glomerular CRD (legacy)
    seu$CRD.bulk0.01 <- CRDscore::cal_CRDscore(
      seu_tpm ,
      circadians = circ.genes.bulk0.01 ,
      study.type= "scRNAseq" )

    # podocyte specific CRD
    seu$CRD <- CRDscore::cal_CRDscore(
      seu_tpm ,
      circadians = podo.circ.0.01 ,
      study.type= "scRNAseq" )

    return(seu)
  })
  names(listSCSN.CRD) <- names(listSCSN.1K.sampl)

  ### combine nephritis days in one study
  listSCSN.CRD.mod <- listSCSN.CRD
  listSCSN.CRD.mod[[8]] <- merge( listSCSN.CRD.mod[[8]] , listSCSN.CRD.mod[[9]])
  listSCSN.CRD.mod[[7]] <- merge( listSCSN.CRD.mod[[7]] ,
                                  subset( listSCSN.CRD.mod[[8]] , subset= gtypeDE=="control" ))
  listSCSN.CRD.mod <- listSCSN.CRD.mod[1:8]
  names(listSCSN.CRD.mod) <- c("Nphs2","Wt1","Pdss2","Lmx1b",
                               "btbr","cd2ap","doxo","nephritis")
  saveRDS( listSCSN.CRD.mod , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/listSCSN.CRDs_08.05.25.rda")

  listSCSN.CRD.mod <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/listSCSN.CRDs_08.05.25.rda")
  
  listSCSN.CRD_table <- Reduce( rbind, lapply( seq(listSCSN.CRD.mod), function(ii){
    XX <- listSCSN.CRD.mod[[ii]]@meta.data[ , c( 
      "age" ,  "gtype"  , "sample", "group", "gtypeDE","sampleDiscr",  "PDS" , 
      grep("CRD",colnames( listSCSN.CRD.mod[[ii]]@meta.data ),ignore.case = T,value = T) )]
    XX$dataSet <- names(listSCSN.CRD.mod)[ii]
    return(XX)
  }))
  saveRDS( listSCSN.CRD_table , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/listSCSN.CRDs_table_08.05.25.rda")
  # circ.genes <- circ.genes[ circ.genes%in% allPodoGenes]
  }

### plot CRD
{
  listSCSN.CRD_table <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/listSCSN.CRDs_table_08.05.25.rda" )
    toPlot <- listSCSN.CRD_table
    vvc <-  grep("CRD",colnames( toPlot ),  ignore.case = T,value = T) 
    ggl <- lapply( seq(vvc), function(ii)
      {
      
      vv <- vvc[ii]
      gg1 <- ggplot( data= toPlot,
                     aes( y = toPlot[[vv]] , 
                          x = PDS,color=gtypeDE ) ) +
        geom_point( alpha=0.3) + 
        # geom_smooth(method = lm) + 
        geom_smooth(method = lm,  se = F) +
        theme_bw() + theme( text = element_text(size=20), legend.position = "bottom") +
        coord_cartesian(ylim = c( quantile( toPlot[[ vv]], 0.025),
                                  quantile( toPlot[[ vv ]], 1)))+
        stat_cor( aes(color=gtypeDE) ,r.accuracy = 0.01,
                  method = "spearman", size=6,cor.coef.name = "rho")+
        scale_color_colorblind()+ labs(y =  vv ) + 
        facet_grid(rows = vars(dataSet))
      
      gg2 <- ggplot( data= toPlot, aes( 
        y = toPlot[[ vv ]] , 
        x = PDS, 
        color = sampleDiscr)  ) +
        geom_smooth(method = lm, aes(  linetype=gtypeDE ), se = F, lwd=3) +
        theme_bw() + theme( text = element_text(size=20), 
                            legend.position = "none") + labs(y = vv )+
        stat_cor(  aes(color=sampleDiscr) , r.accuracy = 0.01,
                   method = "spearman", 
                   label.y = 0 ,size=6, cor.coef.name = "rho")  + 
        facet_grid( rows = vars(dataSet), scales = "free_y") 
      
      return(list(gg1,gg2))
    })
   

   gglgrid <- cowplot::plot_grid(plotlist = unlist(ggl, recursive = F) , ncol =  4, rel_widths = c(0.4,0.6))
  
  pdf( height = 20, width = 30 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/CRDs.vs.PDS_listSCSN_scatter.pdf")
  gglgrid
  dev.off(  )
  
  png( height = 1500, width = 2000 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/CRDs.vs.PDS_listSCSN_scatter.png")
  gglgrid
  dev.off(  )
  
  
  ### make a dot plot that summarises results across studies
  CRDvsPDS.spcor <- Reduce( rbind , lapply( seq( unique(listSCSN.CRD_table$sample)), function(ii) {
    samp <- unique(listSCSN.CRD_table$sample)[ii]
    datt <- listSCSN.CRD_table[ listSCSN.CRD_table$sample == samp, ]
    res <- tryCatch( cor.test( datt$PDS , datt$CRD, method = "spearman"),
                     error = function(e) NA)
    if(is.na(res)) c(NA,NA) else c(res$p.value , res$estimate)
  }))
  
  CRDvsPDS.spcor_tab <- data.frame( SpCor = CRDvsPDS.spcor[,2],
                                    pval.SpCor = CRDvsPDS.spcor[,1],
                                    listSCSN.CRD_table[ match(
                                      unique(listSCSN.CRD_table$sample),
                                      listSCSN.CRD_table$sample),])
  CRDvsPDS.spcor_tab$sig <- ifelse(CRDvsPDS.spcor_tab$pval.SpCor < 0.05 , 
                                   "significant", "non-significant")
  CRDvsPDS.spcor_tab$seqType <- ifelse( CRDvsPDS.spcor_tab$dataSet %in%c("Nphs2","Wt1","Pdss2","Lmx1b") , "single nucleus","single cell")
  CRDvsPDS.spcor_tab <- CRDvsPDS.spcor_tab[ !is.na(CRDvsPDS.spcor_tab$pval.SpCor),]
  gg <- ggplot(CRDvsPDS.spcor_tab , aes(x=gtypeDE, y=SpCor)) +
    geom_point( aes(color=sig, shape=seqType), alpha = 0.5, size=3,
                position=position_jitter(width =0.1) )+
    theme_minimal()+theme( text=element_text( size = 20))+
    scale_color_colorblind()
  
  pdf( height = 5, width = 5 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/CRDs.vs.PDS_listSCSN_dotPlot.sum.pdf")
  gg
  dev.off(  )
}



#### use tempo to estimate circadian time ####
### prepare data for tempo
  {
  ## save podocine data 
  nphs2.podo_seu <-listSCSN$Nphs2 
  # nphs2.podo_seu <- listSCSN$Nphs2
  Idents(nphs2.podo_seu) <- "sample"
  set.seed(42)
  nphs2.podo_seu <- subset( nphs2.podo_seu , downsample=500)
  nphs2.podo_seu$sampleDiscr <- paste0( nphs2.podo_seu$sample,"_",
                                        nphs2.podo_seu$gtypeDE)
  
  
  # save samples in seperate files
  lapply( names(table(nphs2.podo_seu$sample)) , function(smplN){
    datt <- subset( nphs2.podo_seu, subset= sample== smplN )
    datt <-  CreateSeuratObject( datt@assays$RNA@counts[ 
      rownames(datt) %in% allPodoGenes,  ])
    
    ffname <- paste0( "KFO.snRNAseq.Nphs2_",smplN,".h5Seurat")
    SeuratDisk::SaveH5Seurat( datt , filename = ffname )
    SeuratDisk::Convert( ffname , dest = "h5ad" , verbose = T )
    
  })
  
  # datt <- CreateSeuratObject( listSCSN$Lmx1b@assays$RNA@counts[ allPodoGenes,  ])
  # SeuratDisk::SaveH5Seurat( datt , filename = "KFO.snRNAseq.Lmx1b.podoExpr.h5Seurat")
  # SeuratDisk::Convert( "KFO.snRNAseq.Lmx1b.podoExpr.h5Seurat", dest = "h5ad" )
  
}

### read in tempo results
# tempo.out.core.(~20 genes) performs reasonably, 
# while using more genes gives strange artifacts
lldir <- list.dirs( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/tempo.out/Nphs2/KFO.snRNAseq.Nphs2_tempo.out.core.xtnd/",
                    full.names = T , recursive = F )
tempo.out.list <- lapply( seq(lldir) , function( ii ){
  tempo.out <-  tryCatch( read.table( header = T , row.names = 1 , sep = "\t" ,
                            file=paste0( lldir[[ii]] , "/tempo_results/opt/cell_posterior.tsv") ),
                          error = function(e) NA )
  
})
names(tempo.out.list) <- names( table(nphs2.podo_seu$sampleDiscr) )
tempo.out.list <- tempo.out.list[!is.na(tempo.out.list)]

tempo.out.maxP <- Reduce( rbind, lapply( seq(tempo.out.list), function(ii){
  # toPlot <- reshape2::melt( tempo.out.list[[ii]] )
  toPlot <- data.frame ( value= rowMaxs( as.matrix( tempo.out.list[[ii]] )) )
  toPlot$sample <- names(tempo.out.list)[ii]
  toPlot$gtypeCol <- sub(".*_","",names(tempo.out.list)[ii])
  return(toPlot)
  
}))

# # plot max P distribution
# gg <- ggplot( data=tempo.out.maxP , aes(value, color=gtypeCol) ) +
#   geom_density( lwd=1.5 )+ 
#   # geom_histogram(bins = 10) + 
#   # coord_cartesian(xlim = c(min(gglist$value), 0.4))+
#   theme_bw()+ facet_wrap(~sample ,ncol = 1) + 
#   scale_color_colorblind()+
#   stat_summary(aes(xintercept = ..x.., y = 0), fun = median, 
#                geom = "vline", orientation = "y", lwd=1.5 , linetype="dashed")
# gg

# pdf(width = 5, height = 15, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/tempo/Nphs2.samples_tempo.out_posterior.max.Density.pdf")
# gg
# dev.off()

### plot max P correlation with PDS
tempo.out.maxP$PDS <- nphs2.podo_seu$PDS[ rownames(tempo.out.maxP)]
toPlot <- tempo.out.maxP
# toPlot <- tempo.out.maxP[tempo.out.maxP$value<0.5,]
ggplot( data=toPlot , aes(y=value, x=PDS, color=sample) ) +
  geom_smooth(  method ="lm" ,lwd=1.5 ) + 
   # coord_cartesian(ylim = c(0, 0.3))+
  # geom_histogram(bins = 10) + 
  theme_bw()+ facet_wrap(~gtypeCol ,ncol = 3) + 
  # scale_color_colorblind()+
  stat_cor( 
    # p.accuracy = 0.0001,
             r.accuracy = 0.01,
            method = "spearman", size=6, cor.coef.name = "rho")


### plot time uncertainty VS PDS
tempo.out.maxP.bin <- Reduce( rbind, lapply( seq(tempo.out.list), function(ii){
  # toPlot <- reshape2::melt( tempo.out.list[[ii]] )
  toPlot <- data.frame ( value= apply( tempo.out.list[[ii]] , 1, function(X){
    which( X == max(X))} ) )
  toPlot$sample <- names(tempo.out.list)[ii]
  toPlot$gtypeCol <- sub(".*_","",names(tempo.out.list)[ii])
  return(toPlot)
}))
tempo.out.maxP.bin$PDS <- nphs2.podo_seu$PDS[ rownames(tempo.out.maxP)]

### plot time of maximum P
toPlot <- tempo.out.maxP.bin[ tempo.out.maxP.bin$sample=="140740_experimental",]
ggplot( data= toPlot , aes(y=value, x=PDS) ) +
  # geom_smooth(  method = "glm",lwd=1.5 ) + 
  geom_point(alpha = 0.5) + 
  theme_bw() 
  # scale_color_colorblind()+
  # stat_cor( p.accuracy = 0.0001, r.accuracy = 0.1,
            # method = "spearman", size=6, cor.coef.name = "rho")

toPlot2 <- tempo.out.maxP[ tempo.out.maxP$sample=="140740_experimental",]

ggplot( data= toPlot2 , aes(y=value, x=PDS) ) +
  # geom_smooth(  method = "glm",lwd=1.5 ) + 
  geom_point(alpha = 0.5) + 
  theme_bw() 
  # scale_color_colorblind()+
  # stat_cor( p.accuracy = 0.0001, r.accuracy = 0.1,
  #           method = "spearman", size=6, cor.coef.name = "rho")

