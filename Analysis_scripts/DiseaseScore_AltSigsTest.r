# ====================================================================================== #
# ### script for testing signatures of various origins as podocyte damage estimators ### #
# ====================================================================================== #

mallinfo::malloc.trim()
gc()

options( connectionObserver = NULL )
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))

library( AUCell)
library( ggplot2)
library( ggrepel)
library( ggthemes)
library( plyr)
library( biomaRt)
library( Seurat)
library( ggpubr)
library( org.Mm.eg.db)

#### load functions and data
# generate damage signature
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
# calculate damage score
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")

DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")

allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

### load subsampled sc data
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
listSCSN.1K.sampl.AltSig <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23_AltSigs.rda")


# #### prepare MDS and AGS, ctype  signatures #### 
# 
#   # load mesenchimal drift signature
#   MTsig <- readxl::read_xlsx( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/Table_S2_MT.signature.xlsx")
#   gene_symbol.mouse <- fun_homoTO.FROMmouse( MTsig$gene_symbol, TO=T )
#   MTsig.mouse <- MTsig
#   MTsig.mouse<- cbind(  MTsig.mouse , gene_symbol.mouse[
#     match( MTsig.mouse$gene_symbol , rownames(gene_symbol.mouse )),] )
#   MTsig.mouse <- MTsig.mouse[MTsig.mouse$MUS %in% allPodoGenes , ]
#   MTsig.mouse$mean_rank <- sample(1:nrow(MTsig.mouse) )
#   MTsig.mouse$MUS <- unlist( MTsig.mouse$MUS)
#   
#   ### podo aging signature
#   podo.ags <- readxl::read_xlsx("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/NIHMS1607271-supplement-2.xlsx")
#   podo.ags <- as.data.frame( podo.ags[ podo.ags$gene.symbol %in% allPodoGenes , ] )
#   
#   podo.ags$mean_rank <- rank(podo.ags$pval)
#   podo.ags$direction_foldchange <- ifelse( podo.ags$log2FoldChange<0, -1, 1)
#   
#   ### podo ctype signature
#   # prepare gene sets
#   marks.hmphr.sn <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2020_mouse_snRNA.csv")
#   marks.hmphr.sn <- split(marks.hmphr.sn$gene, marks.hmphr.sn$celltype)
#   marks.hmphr.sn <-  data.frame( gName= marks.hmphr.sn$Pod,
#                                  rank= 1:length( marks.hmphr.sn$Pod),
#                                  row.names = marks.hmphr.sn$Pod )
#   marks.hmphr.sc <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2022_mouse_scRNA.csv")
#   marks.hmphr.sc <- split(marks.hmphr.sc$gene, marks.hmphr.sc$celltype)
#   marks.hmphr.sc <-  data.frame( gName= marks.hmphr.sc$Pod,
#                                  rank= 1:length( marks.hmphr.sc$Pod),
#                                  row.names = marks.hmphr.sc$Pod )
#   marks.szt.sc <- read.table(header = T, sep = "\t","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Susztak2021_mouse_scRNA.csv")
#   marks.szt.sc <- split(marks.szt.sc$gene, marks.szt.sc$celltype)
#   marks.szt.sc <-  data.frame( gName= marks.szt.sc$Pod,
#                                  rank= 1:length( marks.szt.sc$Pod),
#                                row.names = marks.szt.sc$Pod )
#   
#   # create a damage signature from a podocyte marker list
  # gg <- union(  union( marks.hmphr.sn$gName , marks.hmphr.sc$gName), marks.szt.sc$gName)
  # CTS <- cbind( marks.hmphr.sn$rank[ match( gg , rownames(marks.hmphr.sn))],
  #               marks.hmphr.sc$rank[ match( gg , rownames(marks.hmphr.sc))],
  #               marks.szt.sc$rank[ match( gg , rownames(marks.szt.sc))])
  # 
  # CTS <- data.frame( gene_symbol= gg ,
  #                    mean_rank= rowMeans(CTS, na.rm = T),
  #                    direction_foldchange = -1 )
  # CTS <- CTS[CTS$gene_symbol %in% allPodoGenes, ]
   # 
  # save( podo.ags ,  MTsig.mouse , CTS , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/MDS.PAGS.CTS.sigs.Rdata")

### podo age  GSE136138
{
  library(GEOquery)
  gsm <- GEOquery::getGEO("GSE136138")
  annot <- (pData( phenoData(gsm$GSE136138_series_matrix.txt.gz) ))
  annot$sample <- gsub("-", ".", annot$title)
  annot$age <- annot$`age:ch1`
  
  age.podo.counts <- read.table( header = T, sep = ",",gzfile("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/GSE136138_Counts.csv.gz"))
  rownames(age.podo.counts) <- make.unique( age.podo.counts$Gene)
  age.podo.counts <- age.podo.counts[, -c(1,2)]
  
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix( countData = age.podo.counts ,
                                 colData = annot ,
                                 design= ~ age )
  dds.r <- DESeq(dds)
  dds.tab <- results(dds.r)
    
  podo.ags2 <- dds.tab[ dds.tab$padj < 0.1 & !is.na(dds.tab$padj) &
                          rownames(dds.tab) %in% allPodoGenes , ]
    podo.ags2$mean_rank <- rank(podo.ags2$padj)
    podo.ags2$gene_symbol <- rownames(podo.ags2)
    podo.ags2$direction_foldchange <- ifelse( podo.ags2$log2FoldChange<0, 1, -1)
    podo.ags2 <- podo.ags2[ order(podo.ags2$mean_rank),]
}


load( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/MDS.PAGS.CTS.sigs.Rdata")

#### test overlap of signatures #### 
  # load circadian genes 
  # glom.circ <- readxl::read_xlsx("PROJECTS/PODOCYTE/DiseaseScore/circadian/Glomerular circadian transcriptomics meta2d BHQ _0.1.xlsx")
  # glom.circ <- glom.circ[1:375, ]
  podo.circ <- readxl::read_xlsx("PROJECTS/PODOCYTE/DiseaseScore/circadian/Podocyte circadian transcriptomics meta2d BHQ _0.1 .xlsx")
  
  # 
  
  geneset0 <- DS_all$gene_symbol[1:42]
  geneset1 <- MTsig.mouse$MUS
  geneset2 <- CTS$gene_symbol
  geneset33 <- podo.ags$gene.symbol
  geneset3 <- podo.ags$gene.symbol[ podo.ags$log2FoldChange <0 ]
  geneset4 <-  podo.ags$gene.symbol[ podo.ags$log2FoldChange >0 ]
  geneset5 <-  podo.circ$Gene_Symbol[ podo.circ$Gene_Symbol %in% 
                                        allPodoGenes & podo.circ$JTK_BH.Q< 0.01] 
  
  
  library(GeneOverlap)
  ggsets <- list(geneset0 , geneset1 , geneset2,
                 geneset3, geneset4 ,
                 geneset5)
  
  names(ggsets) <- c( "PDS" , "MDS" , "CTS" , 
                      "PAGS.down" , "PAGS.up" ,
                      "circ.Podo.01" )

  gom.obj <- newGOM( ggsets, genome.size = length(allPodoGenes))
  # drawHeatmap( XX, cutoff = 1, adj.p = F)
   drawHeatmap( XX , what = "Jaccard" , log.scale= T ,cutoff = 1)
  drawHeatmap( XX , what = "odds.ratio" , log.scale= T  ,cutoff = 1)
  gom.p <- getMatrix(gom.obj, name="pval")
  gom.or <- getMatrix(gom.obj, name="odds.ratio")
  gom.or[ gom.or==0] <- NA
  options("scipen"=-100, "digits"=4)
  
  library(RColorBrewer)
  
  
  gg0 <- pheatmap::pheatmap( log10( gom.or ), cluster_rows = F, cluster_cols = F,
                        color =(RColorBrewer::brewer.pal(name = "Greens", 9)), 
                      display_numbers = gom.p ,
                      fontsize_number = 14 ,
                      number_color = "red") 
  
  pdf( width = 5, height = 4,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/AltSigs_gSet.ovlp.OR.Heatmap.pdf")
  gg0
  
  dev.off() 
  



#### calculate scores for SCSN  #### 
podo.ags.noPDS <-  podo.ags[ !( podo.ags$gene.symbol %in% DS_all$gene_symbol), ] 
CTS.noPDS <-  CTS[ !( CTS$gene_symbol %in% DS_all$gene_symbol), ] 

podo.ags2.noPDS  <- podo.ags2[ !( podo.ags2$gene_symbol %in% DS_all$gene_symbol), ] 

# gg2 <- intersect(  intersect( marks.hmphr.sn$gName , marks.hmphr.sc$gName), marks.szt.sc$gName)
# CTSmin <- CTS[CTS$gene_symbol %in% gg2, ]
# CTSmin.noPDS <-  CTSmin[ !( CTSmin$gene_symbol %in% DS_all$gene_symbol), ] 

listSCSN.1K.sampl.AltSig <- lapply( seq(listSCSN.1K.sampl.AltSig), function(ii)
  {
    print(ii)
    datt<- listSCSN.1K.sampl.AltSig[[ii]]
    
    exprMatrices <- datt@assays$RNA@counts
    exprMatrices <- exprMatrices[ rowSums(round(exprMatrices) > 0) > 0.01*ncol(exprMatrices) , ]
    
    # 
    # datt$MDS.all <- DS_calc.func( exprMatrices , MTsig.mouse ,
    #                               geneIDname = "MUS",
    #                               ntop= nrow(MTsig.mouse) , ceilThrsh=0.05 , progStat =F ,
    #                               wghtd=T )
    # 
    # datt$AGS.all <- DS_calc.func( exprMatrices , podo.ags , 
    #                               geneIDname = "gene.symbol", 
    #                               ntop= nrow(podo.ags) , ceilThrsh=0.05 , progStat =F ,
    #                               wghtd=T )
    # 
    # 
    # datt$AGS.noPDS <- DS_calc.func( exprMatrices , podo.ags.noPDS, 
    #                                 geneIDname = "gene.symbol", 
    #                                 ntop= nrow(podo.ags.noPDS) , ceilThrsh=0.05 , progStat =F ,
    #                                 wghtd=T )
    # 
    # datt$AGS.noPDS.42 <- DS_calc.func( exprMatrices , podo.ags.noPDS, 
    #                                    geneIDname = "gene.symbol", 
    #                                    ntop= 42 , ceilThrsh=0.05 , progStat =F ,
    #                                    wghtd=T )
    
    # datt$CTS <- DS_calc.func( exprMatrices , CTS, 
    #                           geneIDname = "gene_symbol", 
    #                           ntop= nrow(CTS) , ceilThrsh=0.05 , progStat =F ,
    #                           wghtd=T )
    # 
    # datt$CTS.noPDS <- DS_calc.func( exprMatrices , CTS.noPDS, 
    #                           geneIDname = "gene_symbol", 
    #                           ntop= nrow(CTS.noPDS) , ceilThrsh=0.05 , progStat =F ,
    #                           wghtd=T )
    
    datt$PAGS2 <- DS_calc.func( exprMatrices , podo.ags2,
                                    geneIDname = "gene_symbol",
                                    ntop= nrow(podo.ags2) , ceilThrsh=0.05 , progStat =F ,
                                    wghtd=T )

    datt$PAGS2.noPDS <- DS_calc.func( exprMatrices , podo.ags2.noPDS,
                                       geneIDname = "gene_symbol",
                                       ntop= nrow(podo.ags2.noPDS) , ceilThrsh=0.05 , progStat =F ,
                                       wghtd=T )
    
    # print( names(listSCSN.1K.sampl)[ii] )
    
    # print("PDS VS MDS")  
    # print(cor.test( datt$PDS , datt$MDS, method = "spearman"), sep = " ")
    # print("PDS VS MDS all")  
    # print(cor.test( datt$PDS , datt$MDS.all, method = "spearman"), sep = " ")
    # print("PDS VS AGS")
    # print(cor.test( datt$PDS , datt$AGS, method = "spearman"), sep = " ")
    # print("PDS VS AGS all" )
    # print( cor.test( datt$MDS , datt$AGS.all, method = "spearman"), sep = " ")
    
    
    return(datt)
  })
  
  names(listSCSN.1K.sampl.AltSig)<- names(listSCSN.1K.sampl)
  names(listSCSN.1K.sampl.AltSig)<- c("Nphs2","Wt1","Pdss2","Lmx1b",
                                      "btbr","cd2ap","doxo","nephr.D1",
                                      "nephr.D5")
  saveRDS(listSCSN.1K.sampl.AltSig , "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23_AltSigs.rda")
  
  
#### correlations between scores  #### 
  AltSigCor <- Reduce( rbind, lapply( seq(listSCSN.1K.sampl.AltSig), function(ii)
  {
    print(ii)
    datt<- listSCSN.1K.sampl.AltSig[[ii]]
    
    Reduce( rbind , lapply( seq(unique(datt$sample)), function(jj){
      print(jj)
      
      sampleDat <- datt@meta.data[ datt@meta.data$sample==unique(datt$sample)[jj], ]
      
      tryCatch( cbind.data.frame( Reduce( rbind, list( 
        
       unlist( cor.test( sampleDat$PDS , sampleDat$MDS.all, method = "spearman")[c("p.value","estimate")] ),
        unlist(cor.test( sampleDat$PDS , sampleDat$AGS.all, method = "spearman")[c("p.value","estimate")] ),
        unlist(cor.test( sampleDat$PDS , sampleDat$AGS.noPDS, method = "spearman")[c("p.value","estimate")] ),
        unlist(cor.test( sampleDat$PDS , sampleDat$CTS, method = "spearman")[c("p.value","estimate")] ),
        unlist(cor.test( sampleDat$PDS , sampleDat$CTS.noPDS, method = "spearman")[c("p.value","estimate")] )
        
      )) ,
      compare=c("PDSvsMDS" ,"PDSvsAGS","PDSvsAGS.noPDS",
                "PDSvsCTS" , "PDSvsCTS.noPDS"),
      sample=unique(sampleDat$sample), 
      gtypeDE=unique(sampleDat$gtypeDE) ,
      dataset=names(listSCSN.1K.sampl.AGS)[ii]) , error = function(e) 
        as.data.frame( matrix(NA, nrow = 0, ncol = 4, 
                              dimnames =  list(NULL, c("p.value", "estimate.rho", "sample", "dataset"))))
      ) 
    }) ) 
    
    
  })) 
  
  AltSigCor$padj <- p.adjust(AltSigCor$p.value, method = "fdr")  
  
  AltSigCor$sigg <- ifelse(AltSigCor$p.value< 0.01, TRUE , F)  
  toPlot <- AltSigCor[AltSigCor$compare%in% 
                            c("PDSvsMDS" ,"PDSvsAGS","PDSvsAGS.noPDS",
                              "PDSvsCTS" , "PDSvsCTS.noPDS") ,]
  
  gg <- ggplot( toPlot, aes( x= dataset, y = estimate.rho ,
                             color=sigg, 
                             shape=gtypeDE))+ 
    geom_point(size=5, alpha=0.5)+
    geom_hline(yintercept = 0, color="red")+
    facet_wrap(vars(compare))+ 
    theme_bw() + scale_color_colorblind()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf( width = 8, height = 4,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGS.SpCor_SCSN.dotplot.pdf")
  print(gg)
  dev.off()
  
  png( width = 700, height = 300,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGS.SpCor_SCSN.dotplot.png")
  print(gg)
  dev.off() 


#### test scores on aging  data sets #### 
### test signature on aging DR data
  {
    OldYoungArDr <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Kif3a/Seurat/OldYoungArDr_Seurat.SCT.rds")
    OldYoungArDr_sub <- subset( OldYoungArDr,  downsample=2000,  subset= diet %in% c("AL","DR"))
    
    # DSs <-  lapply( seq(unique(OldYoungArDr_sub$seurat_clusters)) , function(ii) {
    #   
    #   print(ii)
      # datt <- subset( OldYoungArDr_sub , seurat_clusters==unique(OldYoungArDr_sub$seurat_clusters)[ii])
      datt <- subset( OldYoungArDr_sub , seurat_clusters=="Podo")
      
      exprMatrices <- datt@assays$SCT@data
      
      exprMatrices <- exprMatrices[ rowSums(round(exprMatrices) > 0) > 0.01*ncol(exprMatrices) , ]
      
      PDS.dr <- DS_calc.func( exprMatrices , DSignature = DS_all , 
                           ntop=  42 , ceilThrsh=0.1 , progStat =F ,
                           wghtd = T )
      
      MDS.dr <- DS_calc.func( exprMatrices , MTsig.mouse , 
                               geneIDname = "MUS", 
                               ntop= nrow(MTsig.mouse) , ceilThrsh=0.1 , progStat =F ,
                               wghtd=T )
      PAGS.dr <- DS_calc.func( exprMatrices , podo.ags , 
                               geneIDname = "gene.symbol", 
                               ntop= nrow(podo.ags) , ceilThrsh=0.1 , progStat =F ,
                               wghtd=T )
      PAGS.noPDS.dr <- DS_calc.func( exprMatrices , podo.ags.noPDS , 
                               geneIDname = "gene.symbol", 
                               ntop= nrow(podo.ags.noPDS) , ceilThrsh=0.1 , progStat =F ,
                               wghtd=T )
      
      llist <- list(PDS,MDS.all,AGS.all)
    #   names(llist) <- c("PDS","MDS.all","AGS.all")
    #   return(llist)
    # })
    
    DSs_df <- data.frame( PDS=  unlist( lapply( DSs, "[[", 1) ),
                          MDS.all= unlist( lapply( DSs, "[[", 2)),
                          AGS.all=  unlist( lapply( DSs, "[[", 3)))
    
    OldYoungArDr_sub@meta.data <- cbind(OldYoungArDr_sub@meta.data , 
                                        DSs_df[ rownames(OldYoungArDr_sub@meta.data),] )
    
    datt@meta.data<- cbind(datt@meta.data , 
                           PDS= PDS.dr , 
                           MDS=MDS.dr,
                           PAGS=PAGS.dr,
                           PAGS.noPDS=PAGS.noPDS.dr)
    
    # lapply( seq(unique(OldYoungArDr_sub$seurat_clusters)), function(ii){
    #   datt <- OldYoungArDr_sub@meta.data[ OldYoungArDr_sub@meta.data$seurat_clusters==
    #                                         unique(OldYoungArDr_sub$seurat_clusters)[ii]  , ]
    #   
    #   print("PDS VS MDS all")  
    #   print(cor.test( datt$PDS , datt$MDS.all, method = "spearman"), sep = " ")
    #   print("PDS VS AGS all" )
    #   print( cor.test( datt$PDS , datt$AGS.all, method = "spearman"), sep = " ")
    #   print( cor.test( datt$MDS.all , datt$AGS.all, method = "spearman"), sep = " ")
    #   
    # })
    
    toPlot <- OldYoungArDr_sub@meta.data
    toPlot <- datt@meta.data
    
    # toPlot <- toPlot[ !(toPlot$diet %in% c("TX227a","TX227b")), ]
    # toPlot <- toPlot[ (toPlot$diet %in% c("AL","DR")), ]
    # toPlot <- toPlot[ (toPlot$seurat_clusters %in% c("Podo")), ]
    
    gg0 <- ggplot( toPlot, aes( x= condition , y = PDS , color=age))+
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    
    gg2 <- ggplot( toPlot, aes( x= condition , y = MDS , color=age))+ 
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    gg3 <- ggplot( toPlot, aes( x= condition , y = PAGS , color=age))+ 
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    gg4 <- ggplot( toPlot, aes( x= condition , y = PAGS.noPDS , color=age))+ 
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    
    cowplot::plot_grid( plotlist = list( gg0, gg2, gg3,gg4),nrow = 2)
    # geom_hline(yintercept = 0, color="red")+
    # facet_wrap(vars(dataset))+ theme_bw() + scale_color_colorblind()
    
    ggg1 <- ggplot( toPlot, aes( x= condition , y = PDS , color=batch))+ 
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    
    ggg2 <- ggplot( toPlot, aes( x= condition , y = MDS , color=batch))+
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    ggg3 <- ggplot( toPlot, aes( x= condition , y = PAGS , color=batch))+
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
    ggg4 <- ggplot( toPlot, aes( x= condition , y = PAGS.noPDS , color=batch))+
      geom_boxplot()+ facet_wrap(vars(seurat_clusters))+theme_clean()
     cowplot::plot_grid( plotlist = list( ggg1, ggg2, ggg3,ggg4),nrow = 2)
    
    # X1 <- toPlot[toPlot$SID == "SID136207",]
    # X2 <- toPlot[toPlot$SID == "SID118987",]
    
    pdf( width = 12, height = 9,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGSvsMDS_scDR.Podo.plot.pdf")
    cowplot::plot_grid( plotlist = list( ggg1, ggg2, ggg3),nrow = 1)
    dev.off() 
  }
  
#### test PDS and MDS and AGS on aging kidney data, bulk RNAseq 
{
  
  
  # remove one sample with very few counts
  countMat <- read.table( row.names = 1, header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circadian_counts.tsv")
  countMat <- countMat[,-1]
  
  countMat <- countMat[, colnames(countMat)!="SN8580185_21056_806Aligned.sortedByCoord.out.bam"]
  annot <- gsub(".*_|Aligned.*", "", colnames(countMat))
  annot <-data.frame(  ID= colnames(countMat) ,
                       condition= annot,
                       age=c(rep(10,16), rep(80,15)),
                       time= sub("^..","",annot))
  annot$age <- as.factor( annot$age )
  annot$time <- as.numeric( annot$time )
  
  ### calculate damage scores 
  gene.symbols <- select( org.Mm.eg.db, keys = rownames(countMat) , keytype = 'ENSEMBL', columns = 'SYMBOL')
  
  countMat$gNames <- gene.symbols$SYMBOL[ match( rownames(countMat),gene.symbols$ENSEMBL ) ]
  countMat <- countMat[ !is.na(countMat$gNames),]
  countMat <- aggregate(.~gNames, data=countMat, FUN=mean)
  rownames(countMat) <- countMat$gNames
  countMat <- countMat[ , -1]
  exprMatrices <- countMat[ rowSums(round(countMat) > 0) > 0.01*ncol(countMat) , ]
  
  PDS.k <- DS_calc.func( exprMatrices , DSignature = DS_all , 
                       ntop=  42 , ceilThrsh=0.1 , progStat =F ,
                       wghtd = T )
  
  MDS.k <- DS_calc.func( exprMatrices , MTsig.mouse , 
                           geneIDname = "MUS", 
                           ntop= nrow(MTsig.mouse) , ceilThrsh=0.1 , progStat =F ,
                           wghtd=T )
  PAGS.k <- DS_calc.func( exprMatrices , podo.ags , 
                           geneIDname = "gene.symbol", 
                           ntop= nrow(podo.ags) , ceilThrsh=0.1 , progStat =F ,
                           wghtd=T )
  PAGS.noPDS.k <- DS_calc.func( exprMatrices , podo.ags.noPDS , 
                          geneIDname = "gene.symbol", 
                          ntop= nrow(podo.ags.noPDS) , ceilThrsh=0.1 , progStat =F ,
                          wghtd=T )
  annot_DS.k <- cbind( annot ,PDS=PDS.k ,MDS=MDS.k, 
                     PAGS = PAGS.k , PAGS.noPDS = PAGS.noPDS.k)
  
  gg0 <- ggplot( annot_DS.k, aes( x= age , y = PDS , color=age))+
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+ 
    geom_text_repel(aes(label=time))+theme_clean()
  gg2 <- ggplot( annot_DS.k, aes( x= age , y = MDS , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+ 
    geom_text_repel(aes(label=time))+theme_clean()
  gg3 <- ggplot( annot_DS.k, aes( x= age , y = PAGS , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+ 
    geom_text_repel(aes(label=time))+theme_clean()
  gg4 <- ggplot( annot_DS.k, aes( x= age , y = PAGS.noPDS , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+ 
    geom_text_repel(aes(label=time))+theme_clean()
  
   cowplot::plot_grid( plotlist = list( gg0, gg2, gg3,gg4),nrow = 1)
  
  pdf( width = 9, height = 9,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/AltSigs_bulkCircadian.plot.pdf")
  cowplot::plot_grid( plotlist = list( gg0, gg2, gg3,gg4),nrow = 1)
  dev.off()  
  
  
}

### test PDS and MDS and AGS on aging glomerular data, bulk RNAseq
{
  ll <-  list.files(path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/GSE240375_RAW",
                    pattern = "counts" , full.names = T )
  
  
  glomAge_bulk <- Reduce( cbind , lapply( seq(ll), function(ii) {
    read.table( gzfile(ll[[ii]] ), header = T, row.names = 1)[,1, drop=F]
  }))
  colnames(glomAge_bulk) <- sub( "_.*","" , basename( ll ))
  
  gene.symbols <- select( org.Mm.eg.db, keys = rownames(glomAge_bulk) , keytype = 'ENSEMBL', columns = 'SYMBOL')
  
  glomAge_bulk$gNames <- gene.symbols$SYMBOL[ match( rownames(glomAge_bulk),gene.symbols$ENSEMBL ) ]
  glomAge_bulk <- glomAge_bulk[ !is.na(glomAge_bulk$gNames),]
  glomAge_bulk <- aggregate(.~gNames, data=glomAge_bulk, FUN=mean)
  rownames(glomAge_bulk) <- glomAge_bulk$gNames
  glomAge_bulk <- glomAge_bulk[ , -1]
  exprMatrices <- glomAge_bulk[ rowSums(round(glomAge_bulk) > 0) > 0.01*ncol(glomAge_bulk) , ]
  
  PDS <- DS_calc.func( exprMatrices , DSignature = DS_all , 
                       ntop=  42 , ceilThrsh=0.1 , progStat =F ,
                       wghtd = T )
  
  MDS.all <- DS_calc.func( exprMatrices , MTsig.mouse , 
                           geneIDname = "MUS", 
                           ntop= nrow(MTsig.mouse) , ceilThrsh=0.1 , progStat =F ,
                           wghtd=T )
  PAGS <- DS_calc.func( exprMatrices , podo.ags , 
                           geneIDname = "gene.symbol", 
                           ntop= nrow(podo.ags) , ceilThrsh=0.1 , progStat =F ,
                           wghtd=T )
  PAGS.noPDS <- DS_calc.func( exprMatrices , podo.ags.noPDS , 
                        geneIDname = "gene.symbol", 
                        ntop= nrow(podo.ags.noPDS) , ceilThrsh=0.1 , progStat =F ,
                        wghtd=T )
  
  glomAge_DS <- cbind.data.frame( PDS=PDS ,MDS=MDS.all , 
                                  PAGS = PAGS , 
                                  PAGS.noPDS= PAGS.noPDS,
                                  age=c(rep("young",3),rep("old",3) ) )
  
  zz0 <- ggplot( glomAge_DS, aes( x= age , y = PDS , color=age))+
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+   
    theme_clean()
  zz2 <- ggplot( glomAge_DS, aes( x= age , y = MDS , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+theme_clean()
  zz3 <- ggplot( glomAge_DS, aes( x= age , y = PAGS , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+theme_clean()
  zz4 <- ggplot( glomAge_DS, aes( x= age , y = PAGS.noPDS , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.5)+ stat_compare_means()+theme_clean()
  
   cowplot::plot_grid( plotlist = list( zz0, zz2, zz3, zz4),nrow = 1)
  
  pdf( width = 9, height = 6,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGSvsMDS_bulkGlom.plot.pdf")
  cowplot::plot_grid( plotlist = list( zz0, zz2, zz3),nrow = 1)
  dev.off() 
}

### test PDS and MDS and AGS on aging glomerular data, sc RNAseq
{
  library( annotation)
  
  
  ll <-  list.files(path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/GSE240375_RAW",
                    pattern = "matrix.mtx" , full.names = T )
  XX <- lapply( seq(ll), function(ii){
    ReadMtx( mtx = ll[[ii]] ,cells= sub( "matrix.mtx", "barcodes.tsv",ll[[ii]]), 
             features =  sub( "matrix.mtx", "features.tsv",ll[[ii]]) , feature.column = 2)
  })
  
  glomAge_sc <- CreateSeuratObject(XX, )
  glomAge_sc$sample <- c( rep( "GSM7697106", ncol(XX[[1]])) ,
                          rep( "GSM7697107", ncol(XX[[2]])) ,
                          rep( "GSM7697108", ncol(XX[[3]])) , 
                          rep( "GSM7697109", ncol(XX[[4]]))  )
  # VlnPlot(glomAge_sc, features = "nCount_RNA", group.by = "sample")
  glomAge_sc <- NormalizeData(glomAge_sc, normalization.method = "LogNormalize", scale.factor = 10000)
  
  glomAge_sc <- FindVariableFeatures(glomAge_sc, selection.method = "vst", nfeatures = 2000)
  glomAge_sc <- ScaleData( glomAge_sc )
  glomAge_sc <- RunPCA(glomAge_sc, features = VariableFeatures(object = glomAge_sc))
  ElbowPlot(glomAge_sc, ndims = 50)
  
  glomAge_sc <- FindNeighbors(glomAge_sc, dims = 1:30)
  glomAge_sc <- FindClusters(glomAge_sc, resolution = 0.1)
  glomAge_sc <- RunUMAP(glomAge_sc, dims = 1:30)
  glomAge_sc$nCount_RNA.log10 <- log10(glomAge_sc$nCount_RNA)
  glomAge_sc[["percent.mt"]] <- PercentageFeatureSet(glomAge_sc, pattern = "^mt-")
  
  
  FeaturePlot( glomAge_sc , features = c("nCount_RNA.log10", "Wt1"))
  DimPlot( glomAge_sc , label = T)
  DimPlot( glomAge_sc, group.by = "sample")
  
  VlnPlot(glomAge_sc,features = c("nCount_RNA.log10", "Wt1","Nphs2",'percent.mt'), pt.size = 0)
  
  glomAge_sc_pod <- subset(glomAge_sc, subset= seurat_clusters %in% c(2,9) )
  
  VlnPlot(glomAge_sc_pod,features = c("nCount_RNA.log10", "Wt1","Nphs2",'percent.mt'), 
          pt.size = 0, group.by="sample")
  
  glomAge_sc_pod.jnd <- JoinLayers(glomAge_sc_pod)
  
  exprMatrices <-glomAge_sc_pod.jnd[["RNA"]]$counts
  exprMatrices <- exprMatrices[ rowSums(round(exprMatrices) > 0) > 0.01*ncol(exprMatrices) , ]
  
  PDS.sc <- DS_calc.func( exprMatrices , DSignature = DS_all , 
                       ntop=  42 , ceilThrsh=0.05 , progStat =F ,
                       wghtd = T )
  
  
  
  MDS.sc <- DS_calc.func( exprMatrices , MTsig.mouse , 
                           geneIDname = "MUS", 
                           ntop= nrow(MTsig.mouse) , ceilThrsh=0.05 , progStat =F ,
                           wghtd=T )
  
  PAGS.sc <- DS_calc.func( exprMatrices , podo.ags , 
                        geneIDname = "gene.symbol", 
                        ntop= nrow(podo.ags) , ceilThrsh=0.1 , progStat =F ,
                        wghtd=T )
  PAGS.noPDS.sc <- DS_calc.func( exprMatrices , podo.ags.noPDS , 
                              geneIDname = "gene.symbol", 
                              ntop= nrow(podo.ags.noPDS) , ceilThrsh=0.1 , progStat =F ,
                              wghtd=T )
  
  glomAge_scDS <- cbind.data.frame( PDS.sc=PDS.sc ,MDS.sc=MDS.sc , 
                                    PAGS.sc = PAGS.sc , PAGS.noPDS.sc=PAGS.noPDS.sc,
                                    age= ifelse( glomAge_sc_pod.jnd$sample %in% 
                                                   c("GSM7697106","GSM7697107"),
                                                 "young", "old"),
                                    sample= glomAge_sc_pod.jnd$sample)
  
  yy0 <- ggplot( glomAge_scDS, aes( x= sample , y = PDS.sc , color=age))+
    geom_boxplot()+ geom_point(size=3, alpha=0.2)+ stat_compare_means()+theme_clean()
  yy2 <- ggplot( glomAge_scDS, aes( x= sample , y = MDS.sc , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.2)+ stat_compare_means()+theme_clean()
  yy3 <- ggplot( glomAge_scDS, aes( x= sample , y = PAGS.sc , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.2)+ stat_compare_means()+theme_clean()
  yy4 <- ggplot( glomAge_scDS, aes( x= sample , y = PAGS.noPDS.sc , color=age))+ 
    geom_boxplot()+ geom_point(size=3, alpha=0.2)+ stat_compare_means()+theme_clean()
  
  cowplot::plot_grid( plotlist = list( yy0, yy2, yy3,yy4),nrow = 1)
  
  pdf( width = 9, height = 6,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGSvsMDS_scGlom.plot.pdf")
  cowplot::plot_grid( plotlist = list( zz0, zz2, zz3),nrow = 1)
  dev.off() 
}
# to do: check aging signature young vs old  in different cell tzpes, the idea is that in podocyte generic aging signature may show the opposite behavior
# to do: check aging signature in bulk circadian and bulk aging data

#### analyse cycling of scores #### 
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
  colnames(dattTest) <- c("Time","Age","MDS.all")
  GLMMfit.MDS <- tryCatch( cglmm(
    MDS.all ~ amp_acro( time_col=Time, period = 24, group = "Age" ),
    data = dattTest ) , error = function(e) NA)
  
  dattTest<- datt.circ.core[, c(1,2, 5)]
  colnames(dattTest) <- c("Time","Age","AGS.all")
  GLMMfit.AGS <- tryCatch( cglmm(
    AGS.all ~ amp_acro( time_col=Time, period = 24, group = "Age" ),
    data = dattTest ) , error = function(e) NA)
  
  p1 <- autoplot(GLMMfit.PDS)
  p2 <-autoplot(GLMMfit.MDS)
  p3 <-autoplot(GLMMfit.AGS)
  
  cowplot::plot_grid( plotlist = list( p1, p2, p3), nrow = 1)
  
  pdf( width = 12, height = 5,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGSvsMDS_bulkKidney.circ.pdf")
  cowplot::plot_grid( plotlist = list( p1, p2, p3), nrow = 1)
  dev.off() 
  
  png( width = 1000, height = 400,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDSvsAGSvsMDS_bulkKidney.circ.png")
  cowplot::plot_grid( plotlist = list( p1, p2, p3), nrow = 1)
  dev.off() 


#### correlate scores with morphology/physiology #### 

  ### SD length
  {
    PDS_tab <- Reduce( rbind, lapply( c(1,3), function(ii){
      datt<- listSCSN.1K.sampl.AltSig[[ii]]@meta.data[ , c(
        "PDS", "MDS.all",  "CTS","CTS.noPDS",
        "AGS.all","AGS.noPDS", "PAGS2","PAGS2.noPDS",
        "gtype", "age", "sample", "group")]
      datt$dataset <- names(listSCSN.1K.sampl.AltSig)[ii]
      return(datt)
    }))
    
    PDS_tab <- PDS_tab[ PDS_tab$sample %in% 
                          names(table(PDS_tab$sample)[
                            table(PDS_tab$sample)>10]), ]
    aggPDS <- aggregate( .~sample , FUN = mean , na.rm=T,
                         data= PDS_tab[, c(1:8, 11)] )
    aggPDS <- cbind(aggPDS , PDS_tab[ match(aggPDS$sample , PDS_tab$sample),
                                      c("gtype", "age", "group","dataset")])
    aggPDS$age_simple <- aggPDS$age
    
    aggPDS$age_simple[aggPDS$age_simple==8] <- "8_9"
    aggPDS$age_simple[aggPDS$age_simple==6] <- "6_7"
    aggPDS$age_simple[aggPDS$age_simple%in%c(12,14)] <- "12_14"
    
    aggPDS$age_simple <- as.factor(aggPDS$age_simple)
    
    gg3 <- ggplot(aggPDS, aes(y=PDS, x=gtype))+ geom_boxplot(outlier.size = 0)+
      geom_point( position = position_dodge(0.1), alpha = 0.5, 
                  aes( color=age_simple, shape=dataset) ,
                  size=3)+
      # stat_compare_means()+
      theme_bw() + theme( 
        text = element_text( size = 14))
    
    ### combine SD and PDS
    aggPDS.SD <- merge( aggregate( .~age_simple+gtype+dataset , FUN = mean , na.rm=T,
                                   data= aggPDS[,c("age_simple","gtype","PDS","MDS.all",
                                                   "CTS","CTS.noPDS", 
                                                   "AGS.all" , "AGS.noPDS",
                                                   "PAGS2","PAGS2.noPDS", "dataset")] ) , 
                        aggregate( .~age_simple+gtype+dataset , FUN = mean , na.rm=T,
                                   data= SD.both[,c("age_simple","gtype","average_SD","dataset")] ),  )
    
    scoresS <- c("PDS","MDS.all","CTS","CTS.noPDS", "AGS.all" , "AGS.noPDS" )
    scoresS <- c("PDS", "AGS.all" , "AGS.noPDS","PAGS2","PAGS2.noPDS" )
    
    aggPDS.SD$age <- as.numeric( sub( "_" ,"." , aggPDS.SD$age_simple)) 
    aggPDS$ageS <- as.numeric( sub( "_" ,"." , aggPDS$age_simple)) 
    
    ggl <-  lapply( seq(scoresS), function(ii){
      ggplot(aggPDS.SD, aes(y=average_SD, x=aggPDS.SD[,scoresS[ii]] ))+
        geom_point(  alpha = 0.5,
                     aes( color=gtype , shape=dataset) ,
                     size=3)+stat_cor()+
        geom_smooth(method = "lm")+ xlab(scoresS[ii])+
        theme_bw() + theme(  text = element_text( size = 14))
      
      # ggplot(aggPDS, aes(y= age, x=aggPDS[,scoresS[ii]] ))+
      #   geom_point(  alpha = 0.5,
      #                aes( color=gtype , shape=dataset) ,
      #                size=3)+stat_cor()+
      #   geom_smooth(method = "lm")+ xlab(scoresS[ii])+
      #   theme_bw() + theme(  text = element_text( size = 14))
      
    })
    
    cowplot::plot_grid( plotlist = ggl , nrow = 1)
    
    gll <- cowplot::plot_grid( plotlist = ggl , nrow = 2)
    
    pdf( height = 8, width = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsSDlength/PDS.MDS.AGSvsSD_Nphs2Pdss2.pdf")
    print(gll)
    dev.off()
    
    png( height = 600, width = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsSDlength/PDS.MDS.AGSvsSD_Nphs2Pdss2.png")
    print(gll)
    dev.off()
  }
  
  ### AlbCr
  {
    #### bulk ####
    
    stud <- c( "GSE117571", "GSE108629", "GSE17709" ,
               "GSE117987", "GSE126217" ,"GSE154955",
               "GSE112116", "GSE131266", "GSE134327", 
               "GSE110092" ,"GSE77717" , "KFO.Wt1"   )
    
    annot_bulkStage <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/annot_bulkStage.rda")
    
    # calculate PDS
    blkData <- explist[ names(explist) %in% stud]
    PDSsize.bulk <-  lapply( seq(blkData) , function(ii , ceilThrsh = 0.05  )
    {
      DS_calc.func( exprMatrices = blkData[[ii]], 
                    ceilThrsh = ceilThrsh , 
                    wghtd = F, useOrder = "mean_rank",
                    DSignature= DS_all, ntop = 42)
    })
    names(PDSsize.bulk) <- names(blkData)
    # order like in the stud vector
    PDSsize.bulk <- PDSsize.bulk[stud]
    # GSE117571 use only gloms
    PDSsize.bulk[["GSE117571"]] <- PDSsize.bulk[["GSE117571"]][1:4]
    
    # combine in one df
    ll<- lapply( seq(PDSsize.bulk), function(jj)
    {
      score <- PDSsize.bulk[[jj]]
      id <- names(PDSsize.bulk)[jj]
      # combine score and annotation, then aggregate the score by annotation
      # if a multistaged study - use the prepared annotation list,
      # otherwise extract annotation from sample names
      if( id %in% names(annot_bulkStage)) {
        annot <- annot_bulkStage[[id]]$groups[ 
          match( sub( "__.*", "", names(score)), 
                 rownames( annot_bulkStage[[id]] ))]
      } else {
        annot <- sub( ".*__", "", names(score))
      }
      
      XX <- data.frame(score =score, study  = id,
                       groups= as.factor(annot) )
      XX <- aggregate( .~groups+study, data=XX , FUN=median)
      
      return(XX)
    })
    names(ll) <- stud 
    datMean  <- Reduce( rbind, ll)
    
    # ad proteinuria
    mgmg.Bulk <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/mgmg.Bulk.rda")
    datMean$mgmg <-  unlist(mgmg.Bulk) # group  Wt1h.d. KFO
    
    # add platform type annotation
    datMean$platform <- c( rep( "MA" , 2) , rep( "MA" ,3 ) ,  rep( "MA" , 2) , 
                           rep(  "bulk" , 2) , rep( "bulk", 2) ,  rep( "bulk" , 3),
                           rep(  "MA" , 2) , rep( "MA", 2) , rep( "bulk", 2) ,  
                           rep( "bulk", 2) , rep(  "bulk", 2), rep(  "bulk", 4))
    ## add stage
    datMean$stage <- ifelse( datMean$groups %in% c("experiment","mutant", "ko_4w","D9 after ADR injection","LMB2day4"),
                             "stage.1", ifelse(datMean$groups %in% c( "ko_12w", "D14 after ADR injection","LMB2day7"), 
                                               "stage.2", "control"))
    
    # make a plot 
    gg <- ggplot2::ggplot( data = datMean, aes( x=score , y=log(mgmg) )) +
      geom_point( aes( color=study , shape=stage), size=6 )+ theme_bw() +  
      ggtitle(paste("42 genes damage signature",sep = ""))+
      theme( text = element_text(size = 22) )  + 
      # geom_text(hjust=0, vjust=0)+
      geom_smooth(method='lm', se = FALSE) + stat_cor(size=7, method = "spearman") 
    gg
    
    # save the plot
    pdf(height = 6, width = 9, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/PDS.42vsAlbCr_bul.pdf")
    print(gg)
    dev.off()
    # save the plot
    png(height = 500, width = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/PDS.42vsAlbCr_bul.png")
    print(gg)
    dev.off()
    
    
    
    #### sn KFO ####
    
    # laod annotation
    library(ggpubr)
    annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
    annot_tab$group <- paste(annot_tab$Genotype,annot_tab$Age_weeks,sep = "_")
    
    # combine metadata from 3 experiments 
    datt <- Reduce( rbind , lapply( listSCSN.1K.sampl.AltSig[1:3], function(X) {
      XX <- X@meta.data 
      return(XX[,c( "group","sample","gtype", "PDS",  "MDS.all", 
                    "CTS" , "CTS.noPDS","AGS.all" , "AGS.noPDS",
                    "PAGS2","PAGS2.noPDS"
      )])
    }))
    
    # treat carefully 21 week Pdss2 samples since they have only per group measurements
    datt1 <- datt[ datt$sample %in% c( "146985", "146986", "143485" , "143486") ,]
    # aggregate
    aggPDS1 <- aggregate( .~group , FUN = mean , 
                          data= datt1[ , c("group", "PDS",  "MDS.all", 
                                           "CTS" , "CTS.noPDS","AGS.all" , 
                                           "AGS.noPDS" ,"PAGS2","PAGS2.noPDS") ] )
    aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
    aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
    colnames(aggPDS1)[1] <- "sample"
    # the rest of samples
    datt2 <- datt[ !(datt$sample %in%  c("146985", "146986", "143485" , "143486")), ]
    aggPDS2 <- aggregate( .~sample , FUN=mean,
                          data= datt2[ , c( "sample","PDS", "MDS.all", 
                                            "CTS" , "CTS.noPDS","AGS.all" , 
                                            "AGS.noPDS","PAGS2","PAGS2.noPDS" )] ) 
    
    aggPDS <- rbind(aggPDS2, aggPDS1)
    
    aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ match( sub("SID","" ,aggPDS$sample) , 
                                                      annot_tab$CCG_Sample_ID)]
    aggPDS$group <- annot_tab$group[ match( sub("SID","" ,aggPDS$sample ), 
                                            annot_tab$CCG_Sample_ID)]
    aggPDS$gtype <- as.factor(annot_tab$Genotype[ match( sub("SID","" ,aggPDS$sample ), 
                                                         annot_tab$CCG_Sample_ID)])
    
    PDSvec <- c("PDS",  "AGS.all" , "AGS.noPDS","PAGS2","PAGS2.noPDS")
    
    aggPDS$AlbCrRatio <- as.numeric(aggPDS$AlbCrRatio)
    
    gglist <- Reduce(c,  lapply( seq( PDSvec ), 
                                 function(ii){
                                   # plot lm and correlation
                                   gg1 <- ggplot2::ggplot( data = aggPDS, 
                                                           aes( x=aggPDS[,PDSvec[ii]], 
                                                                y=log(AlbCrRatio))) +
                                     geom_point(  size=6, aes(col=gtype)) +
                                     theme_bw() +  theme( text = element_text(size = 22)) + 
                                     geom_smooth(method='lm', se = FALSE) + ggtitle( PDSvec[ii])+ 
                                     stat_cor( size=7, method = "spearman" ) 
                                   # geom_text(aes(label = sample  ), size=6, position = "dodge")
                                   
                                   # density plots 
                                   gg2 <- ggplot2::ggplot( data = datt, 
                                                           aes( x=datt[,PDSvec[ii]], 
                                                                color=gtype)) +
                                     geom_density(size=1.5) + ggtitle( PDSvec[ii])+ 
                                     theme_bw() +  theme( text = element_text(size = 22)) 
                                   
                                   # dotplot for samples of Nphs2mut
                                   aggPDS.nphs2 <- aggPDS[aggPDS$sample%in% listSCSN.1K.sampl.AltSig$Nphs2$sample,]
                                   gg3 <- ggplot2::ggplot( data = aggPDS.nphs2, 
                                                           aes( y=aggPDS.nphs2[,PDSvec[ii]], 
                                                                x=group, 
                                                                color=group)) +
                                     geom_jitter(size=1.5) + ggtitle( PDSvec[ii])+ 
                                     theme_bw() +  theme( text = element_text(size = 22)) +
                                     geom_label(aes(label=sample))
                                   
                                   # density plots for nphs2
                                   datt.nphs2 <- datt[datt$sample%in% listSCSN.1K.sampl.AltSig$Nphs2$sample,]
                                   gg4 <- ggplot2::ggplot( data = datt.nphs2, 
                                                           aes( x=datt.nphs2[,PDSvec[ii]], 
                                                                color=group)) +
                                     geom_density(size=1.5) + ggtitle( PDSvec[ii])+ 
                                     theme_bw() +  theme( text = element_text(size = 22)) 
                                   
                                   return( list(gg1, gg2
                                                # , gg3 ,gg4
                                   ))
                                 }))
    gll0 <-  cowplot::plot_grid(plotlist = gglist , ncol=2)
    
    pdf( height = 32, width = 24, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDS.PAGSvsAlbCr_KFO.pdf")
    print(gll0)
    dev.off()
    
    png( height = 2000, width = 2000, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsMTS/PDS.MDS.AGSvsAlbCr_KFO.png")
    print(gll0)
    dev.off()
    
  }


