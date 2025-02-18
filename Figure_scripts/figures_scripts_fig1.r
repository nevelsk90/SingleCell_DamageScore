###=====# PDS publication figure 1 code #=======###
# an input for each figure should be a data, 
# which requires no significant computations for plotting

# release memory
mallinfo::malloc.trim()

#### load DS and code ####
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))
options( connectionObserver = NULL )
library( org.Mm.eg.db )
library( ggplot2 )
library(cowplot)
library( reshape2 )
library( plyr )
library( Seurat )
library( viridis )
library( ggthemes )
library( AUCell )
library( ggpubr )
library( biomaRt )
library( RColorBrewer )

# setwd
setwd( "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1" )
inputdir <-  "/media/tim_nevelsk/WD_tim/PROJECTS/WRITING/PDS_manuscript/Figure_input" 

# load necessary code
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# load expression data
listSCSN.1K <- readRDS( paste0( inputdir, "/listSCSN_1K.22.12.23.rda") )


#### 1. Scheme of the damage score pipeline, made in inkscape ####

#### 2. Dim.Red. plots showing PDS gradient in the population of damaged podocytes ####
### UMAP of podocyte colored by PDS, from Wt1het.del snRNAseq study
  {
    toPlot <- listSCSN.1K$Wt1
    
   gg1 <- FeaturePlot( toPlot , features = "PDS", 
                        cells = sample(colnames(toPlot)),
                        cols = brewer.pal( n = 9 , name = "YlOrBr"),
                       pt.size =1.5)+ 
      theme(legend.position = "bottom", 
            legend.key.height=unit(5,"mm"),
            text = element_text(size=20),
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) 
      # coord_cartesian(xlim = c(-10,-6), ylim = c(-2,3))+
      # theme(legend.position = "none")
    
    gg2 <- DimPlot( toPlot , group.by = "gtype",
                   pt.size = 1.5 , shuffle = T)+
      # coord_cartesian(xlim = c(-10,-6), ylim = c(-2,3))+
      scale_colour_colorblind()+
      theme(legend.position = "bottom",
            legend.key.height=unit(15,"mm"),
            text = element_text(size=20),
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())
   
    
    ## densito plot for PDS
    gg3 <- ggplot(data=toPlot@meta.data,
                  aes( x=PDS, color=gtype))+
      geom_density(lwd=3)+ theme_bw()+
      scale_colour_colorblind()+
      theme(legend.position = "none",
            axis.title.y =element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            text = element_text(size=20))
    
    ggpnl <-cowplot::plot_grid(plotlist = list(gg2, gg1, gg3), 
                               rel_heights = c(1,0.25), nrow = 2)
    
    
    pdf(height = 6, width = 8, file = "Wt1hd.snRNAseq_podo_PDS.umap.NOdim.pdf")
    ggpnl
    dev.off()
    png(height = 600, width = 600, file = "Wt1hd.snRNAseq_podo_PDS.umap.NOdim.png")
    ggpnl
    dev.off() 
  }


### suppl. fig UMAPS for 3 KFO studies to show the data
{
 
  nphs2 <- readRDS ( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_Seur.scDblFilt.rda")
  wt1 <- readRDS ("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_WT1hetdel/Seurat/snRNAseq.WT1hetdel_decontX.allcells_Seur.4w12w.scDblFilt.rda" )
  pdss2 <- readRDS ( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Cem_data/Seurat/Pdss2CoQ2_decontX.allcells_Seur.scDblFilt.rda" )

  varsPlot <- c("ctype", "groups")
  gg01 <- DimPlot(nphs2, label = T, group.by = "ctype", shuffle = T)
  gg02 <-DimPlot(nphs2, group.by = "group", shuffle = T)
  
  gg11 <- DimPlot(wt1, label = T,  group.by = "ctype", shuffle = T)
  gg12 <-DimPlot(wt1, group.by = "group", shuffle = T)
  
  gg21 <- DimPlot(nphs2, group.by = "ctype", shuffle = T)
  gg22 <-DimPlot(nphs2, group.by = "group", shuffle = T)
  
  ggl<-  cowplot::plot_grid(plotlist = list(gg01,gg02, 
                                            gg11,gg12), ncol = 2)
  
 ###  only podo
  listSCSN <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_22.12.23.rda")
  
  pp01 <- DimPlot( listSCSN$Nphs2, group.by = "group", shuffle = T)+
    coord_cartesian(xlim = c(-13,-6), ylim = c(-1,7))+theme(legend.position = "bottom")
  pp02 <- DimPlot( listSCSN$Wt1, group.by = "group", shuffle = T,)+
    coord_cartesian(xlim = c(-10,-6), ylim = c(-2,3))+theme(legend.position = "bottom")
  pp03 <- DimPlot( listSCSN.1K$Pdss2, group.by = "group", shuffle = T)+
    coord_cartesian(xlim = c(5,14), ylim = c(-3,7))+theme(legend.position = "bottom")
  
  ppl<-  cowplot::plot_grid(plotlist = list(pp01,pp02, pp03), nrow = 1)
  
  }

#### 3. model cross validation plots ####
###  Disease-model leave-one out cross-validation
  {
  CV_model.list_PDS42.agg <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/CV_model.list_PDS42.sampAgg.rda")
  
  # scale PDS
  toPlot.list <- lapply(seq(CV_model.list_PDS42.agg), function(ii){
    toPlot <- CV_model.list_PDS42.agg[[ii]]
    Reduce( rbind, lapply( seq(toPlot), function(jj){
      datt <- toPlot[[jj]]
      datt$PDS.42 <- scale(datt$PDS.42)
      return(datt)
    }))
  }) 
  names(toPlot.list) <- names(CV_model.list_PDS42.agg)
  
  ### plot individual models
  ppglist <- lapply(seq(toPlot.list), function(ii){
    datt <- toPlot.list[[ii]]
    ggplot2::ggplot( datt , 
                     aes( y= PDS.42, x=gtypeDE, color=gtypeDE)) +
      scale_color_colorblind() +  geom_boxplot(lwd=1.5) + 
      geom_jitter(size=3)+
      theme_bw() + theme( text = element_text(size = 24) , axis.title.x=element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
      # stat_compare_means( size = 6   )+
      facet_grid(cols = vars(dataset),scales = "free")
  })
  ppglist <- cowplot::plot_grid( plotlist = ppglist, ncol = 1)
  
  pdf( height = 20, width = 20, file = "Supl.Fig1/PDS_cv.models_boxplot_Indv.pdf" )
  print( ppglist )
  dev.off()
  png( height = 1800, width = 1500, file = "Supl.Fig1/PDS_cv.models_boxplot_Indv.png" )
  print( ppglist )
  dev.off()
  
  ### make one plot for all models
  toPlot <- Reduce( rbind, lapply(seq(toPlot.list), function(ii){
    datt <- toPlot.list[[ii]]
    datt <- aggregate( PDS.42 ~ gtypeDE + dataset + CV_model, data=datt , FUN=mean )
    names(datt)[4]<- "PDS.42"
    return(datt)
  }))
  
  gg<-  ggplot2::ggplot( toPlot , 
                         aes( y= PDS.42, x=gtypeDE, color=gtypeDE)) +
    scale_color_colorblind() +  geom_boxplot(lwd=1.5) + 
    geom_jitter(size=5)+
    theme_bw() + theme( text = element_text(size = 24) , axis.title.x=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means( size = 6   ) +
    facet_grid(cols = vars(CV_model),scales = "free")
  
  pdf( height = 6, width = 12, file = "Figure1/PDS_cv.models_boxplot_agg.pdf" )
  print( gg )
  dev.off()
  
  # png( height = 600, width = 1000, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1/PDS_cv.models_boxplot_agg.png" )
  # print( gg )
  # dev.off()
  }

### Suppl fig. platform CV leave-one out cross-validation
### Suppl fig. DS randomisation test
### Suppl. Signature size test ####
  {
  library(GSEABase)
  DSignature <- DS_all
  # DSignature <- DSignature[order(DSignature$mean_rankSC),]
  
  
  size_vec <- c( 5, 10, 15, 20, 30, 42 , 50, 100, 200, nrow(DSignature) )
  size_Gsets <- Reduce( c , lapply(size_vec, function(ssize){
    
    genesUP = DSignature$gene_symbol[1:ssize][ 
      which( DSignature$direction_foldchange[1:ssize]==1)]
    genesDOWN =  DSignature$gene_symbol[1:ssize][ 
      which( DSignature$direction_foldchange[1:ssize]== -1)]
    
    list( GeneSet( genesUP , 
                   setName= paste0("UP.",ssize )) , 
          GeneSet( genesDOWN , 
                   setName= paste0("DOWN.", ssize)))
  }) )
  size_Gsets <- GSEABase::GeneSetCollection( size_Gsets  )
  
  ### analyse effect of a set size on individual studies
  ceilThrsh <- 0.05
  
  # iterate over sc studies and calculate PDS
  set.seed(42)
  listSCSN.PDS_DSsizeTest <- lapply( seq(listSCSN.1K.sampl),
                                     function(ii)
                                     {
                                       ### uppercase row. names
                                       print(ii)
                                       newSeu <- listSCSN.1K.sampl[[ii]]
                                       # subset 
                                       Idents(newSeu) <- newSeu$sample
                                       
                                       # get counts
                                       expr <- newSeu@assays$RNA@counts
                                       expr <- expr[ rowSums(round(expr) > 0) > 0.01*ncol(expr) , ]
                                       
                                       ###  1. Build gene-expression rankings for each cell  
                                       cellRanks <- AUCell_buildRankings( expr, nCores=4, plotStats=F, 
                                                                          verbose = T )
                                       
                                       ###  2. Calculate enrichment for the gene signatures (AUC)
                                       # gene sets should include UP and DOWN regulated genes seperately
                                       cells_AUC <- AUCell_calcAUC( size_Gsets, cellRanks , verbose = T,
                                                                    # aucMaxRank parameter is a tricky one,
                                                                    # for single cell 0.05 is a reasonable 
                                                                    # but for bulk a bigger, e.g. 0.2 can be used
                                                                    aucMaxRank = ceiling( ceilThrsh * nrow(cellRanks))
                                       )
                                       cells_AUC <- getAUC( cells_AUC) #this will extract results 
                                       
                                       ## make sure even empty gSets are in the matrix
                                       cells_AUC <- cells_AUC[ match(
                                         names(size_Gsets), rownames(cells_AUC)),]
                                       cells_AUC[is.na(cells_AUC)] <- 0
                                       
                                       # calculate coefficients to balance by N of UP and DOWN
                                       coeff <- sapply( seq(size_Gsets), function(jj) length(
                                         size_Gsets[[jj]]@geneIds)/rep(size_vec, each=2)[jj])
                                       ## adjust 
                                       cells_AUC.wghtd <- cells_AUC*coeff
                                       
                                       # subtract DOWN from UP AUCell scores
                                       DSsizeTest = cells_AUC.wghtd[ seq( 1,length(size_Gsets), by=2) ,] - 
                                         cells_AUC.wghtd[seq( 1,length(size_Gsets), by=2)+1,]
                                       rownames(DSsizeTest) <- paste0("PDS.size",size_vec)
                                       
                                       # add PDS to the metadata
                                       newSeu@meta.data <- cbind(newSeu@meta.data, t(DSsizeTest))
                                       return(newSeu)
                                     })
  # saveRDS(listSCSN.PDS_DSsizeTest, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS42_DSsizeTest.rda")
  # saveRDS(listSCSN.PDS_DSsizeTest, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS_DSsizeTest.rda")
  
  ### plot PDS distributions
  {
    ## prepare the table
    toPlot <- lapply( seq(listSCSN.PDS_DSsizeTest) , function(ii){
      datt <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data
      datt <-   newSeu@meta.data
      datt <- datt[,c( "gtypeDE","gtype","group","sample", 
                       grep("PDS",names(datt), value = T))]
      datt[,grep("PDS.*",names(datt), value = T)] <- scale(
        datt[,grep("PDS.*",names(datt), value = T)], scale = T )
      return(datt)
    })
    
    toPlot.melt <- reshape2::melt( toPlot , id.vars=c("gtypeDE","gtype","group","sample"))
    toPlot.melt$L1 <- factor( toPlot.melt$L1, labels =(names(listSCSN.1K.sampl)))
    # remove pdss2 6 w as having no effect, and old PDS
    toPlot.melt <- toPlot.melt[toPlot.melt$group!="Pdss2_6" &
                                 toPlot.melt$variable!="PDS" ,]
    toPlot.melt$variable <- droplevels(toPlot.melt$variable)
    
    ### make plots
    ppg1 <- ggplot( toPlot.melt , aes( y=value , 
                                       x=gtypeDE , 
                                       color=gtypeDE)) +
      geom_violin(trim = T,scale="width") + theme_bw() + scale_color_colorblind() + 
      stat_summary(
        fun.data = "mean_sdl",  fun.args = list(mult = 1), 
        geom = "pointrange", color = "black"
      )+ theme( text = element_text(size=20), legend.position = "bottom", 
                axis.text.x = element_blank()
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      ) +
      facet_grid( cols = vars(variable))
    # # separate each sc.sn study
    # ppg2 <- ppg1+facet_grid(rows = vars(L1), scales = "free")+
    #   stat_summary(fun=median, geom="line",
    #                aes(group=gtypeDE ,color=gtypeDE))
    
    # ##  add t stats
    # t_stat<- sapply( unique(toPlot.melt$variable), function(ii){
    #   tres <- t.test(x= toPlot.melt$value[toPlot.melt$variable==unique(toPlot.melt$variable)[ii] & 
    #                                         toPlot.melt$gtypeDE=="experimental" ],
    #                  y= toPlot.melt$value[toPlot.melt$variable==unique(toPlot.melt$variable)[ii] &
    #                                         toPlot.melt$gtypeDE=="control"] )
    #   tres$statistic
    # })
    # names(t_stat) <- unique(toPlot.melt$variable)
    # # create df of tests
    # test.means <- compare_means(
    #   value ~ gtypeDE, data = toPlot.melt, group.by = "variable",
    #   method = "t.test", ref.group = "control",p.adjust.method = "bonferroni"
    # )
    # # add tstat to the test df
    # test.means$tstat <- round(t_stat,2)
    # # plot
    # ppg3 <- ppg1 +
    #   stat_pvalue_manual(data = test.means, label="tstat",
    #                      x = "variable", 
    #                      y.position = 4 )+ 
    #   scale_x_discrete(limits = rev(levels(toPlot.melt$variable)))+
    #   stat_summary(fun=mean, geom="line",
    #                aes(group=gtypeDE ,color=gtypeDE))+ 
    #   coord_flip()+ facet_grid( cols = vars(variable))
    
    
    
    # save plots for 2 groups, used for DE
    pdf(height = 6, width = 12, file =  paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/SetSize/",
                                               "PDS_sizeTest_scsn.VlnPlot.pdf", sep = "") )
    print(ppg1)
    dev.off()
    
    
  }
  
  
  ### plot correlation with PDS
  {
    # combine metadata from 3 experiments 
    datt <- Reduce( rbind , lapply( 1:3, function(ii) {
      XX <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data 
      return(XX[,c( "gtypeDE","group","sample","gtype", 
                    grep("PDS",names(XX), value = T)
      )])
    }))
    datt <- datt[, colnames(datt)!="PDS.42"]
    toPlot.melt <- reshape2::melt( datt , id.vars=c("gtypeDE","gtype","group","sample"))
    # toPlot.melt$L1 <- factor( toPlot.melt$L1, labels =(names(listSCSN.PDS)))       
    
    # treat carefully 21 week Pdss2 samples since they have only per group measurements
    datt1 <- datt[ datt$sample %in% c( "146985", "146986", "143485" , "143486") ,]
    # aggregate
    aggPDS1 <- aggregate( .~group , FUN = mean , 
                          data= datt1[ , c("group", 
                                           grep("PDS",names(datt), value = T))] )
    aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
    aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
    colnames(aggPDS1)[1] <- "sample"
    # the rest of samples
    datt2 <- datt[ !(datt$sample %in%  c("146985", "146986", "143485" , "143486")), ]
    aggPDS2 <- aggregate( .~sample , FUN=mean,
                          data= datt2[ , c( "sample",
                                            grep("PDS",names(datt), value = T))] ) 
    
    aggPDS <- rbind(aggPDS2, aggPDS1)
    
    aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ match( sub("SID","" ,aggPDS$sample) , 
                                                      annot_tab$CCG_Sample_ID)]
    aggPDS$group <- annot_tab$group[ match( sub("SID","" ,aggPDS$sample ), 
                                            annot_tab$CCG_Sample_ID)]
    aggPDS$gtype <- as.factor(annot_tab$Genotype[ match( sub("SID","" ,aggPDS$sample ), 
                                                         annot_tab$CCG_Sample_ID)])
    # melt 
    aggPDS.melt <- reshape2::melt( aggPDS ,
                                   id.vars= names(aggPDS)[!names(aggPDS) %in%
                                                            grep("PDS",names(aggPDS), value = T)])
    gg <- ggplot2::ggplot( data = aggPDS.melt, aes( x=value, y=log(AlbCrRatio))) +
      geom_point(  size=6, aes(col=gtype)) +
      theme_bw() +  theme( text = element_text(size = 22)) + 
      geom_smooth(method='lm', se = FALSE) + ggtitle( "Nphs2mut podocytes")+ 
      stat_cor( size=5, method = "spearman" )  +  
      geom_text(aes(label = sample  ), size=6, position = "dodge")
    # facet_grid(cols = vars(variable), scales = "free")
    gglist <-  gg+ ggforce::facet_wrap_paginate( ~variable, scales = "free", 
                                                 ncol = 2 ) 
    
    
    # distribution plot for NPHS2
    dattX <- datt[datt$sample%in% listSCSN.PDS$Nphs2$sample,]
    dattX.melt <- reshape2::melt( dattX ,
                                  id.vars= names(dattX)[!names(dattX) %in%
                                                          grep("PDS",names(dattX), value = T)])
    ggplot2::ggplot( data = dattX.melt, aes( x=value, color=group)) +
      geom_density(size=1.5) +
      theme_bw() +  theme( text = element_text(size = 22))+
      ggforce::facet_wrap_paginate( ~variable, scales = "free", 
                                    ncol = 2 )
    
    XX <- aggPDS[aggPDS$sample%in% listSCSN.PDS$Nphs2$sample,]
    XX.melt <- reshape2::melt(XX,  id.vars= names(XX)[!names(XX) %in%
                                                        grep("PDS",names(XX), value = T)])
    ggplot(XX.melt, aes(x=group,y=value, color=group))+geom_jitter()+theme_bw()+
      geom_label(aes(label=sample))+
      ggforce::facet_wrap_paginate( ~variable, scales = "free", 
                                    ncol = 2 )
  }
  
}
### Suppl fig. sling unsupervised sc/sn trajectory
  {  
  library( scater )
  library(slingshot)
  
  sce.sling <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.1K_sling.ptime_sce.rda" )
  
  embedded <- embedCurves(sce.sling, "UMAP.CCA")
  embedded <- slingCurves(embedded)[[1]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord,])
  
  # plot with pseudotime
  gg1 <- plotReducedDim( sce.sling, colour_by=("slingPseudotime_1"),
                   dimred = "UMAP.CCA") + ggtitle( "trajectory ptime" ) +
    geom_path( data=embedded, aes( x=umapcca_1, y=umapcca_2), size = 1.2 )+
    theme( text = element_text( size =24) , legend.position = "bottom")
  # xlim(-13,-8)
  # plot with PDS
  gg2 <-plotReducedDim( sce_clust_SCE , colour_by=("PDS"),
                  dimred= "UMAP.CCA" ) + ggtitle("aucell.42 PDS")+
    theme( text = element_text( size =24) , legend.position = "bottom")
  # plot conditions
  gg0 <-  plotReducedDim( sce_clust_SCE , colour_by=("gtypeDE") ,
                    text_by= "groupSling" , text_colour="red" , text_size = 8,
                    dimred = "UMAP.CCA" ) + 
    ggtitle( "experimental conditions" )+
    theme( text = element_text( size =24), legend.position = "bottom")
    # xlim(-13,-8)+ggtitle("Wt1het.del. podocytes")
  ggl <- cowplot::plot_grid( gg0 , gg1 , gg2 , nrow = 1 )
  
  # save plots
  pdf(height = 6, width = 15, file = "Supl.Fig1/listSCSN.1K_sling.ptimeDE_dimRed.pdf")
    ggl
  dev.off()
  png(height = 600, width = 1500, file = "Supl.Fig1/listSCSN.1K_sling.ptimeDE_dimRed.png")
    ggl
  dev.off() 
  
  ### DE along the trajectory
  pseudo <- TSCAN::testPseudotime( sce.sling, 
                            pseudotime=sce.sling$slingPseudotime_1)
  pseudo$Gene.Symbol <- rownames( pseudo)
  pseudo <- pseudo[order(pseudo$p.value),]
  colnames(pseudo)[1:2] <-  c("log2FoldChange", "pvalue" )
  
  saveRDS( pseudo, file = "listSCSN.1K_sling.ptimeDE.rda")
  
}
### Suppl fig. DE LFC heatmap of damage markers
  {
  ### load DE results
  FSGS_MA_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/MA_DElist.04.04.22.rda")
  FSGS_bulk_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/bulk_DElist.04.04.22.rda")
  FSGS_sc_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/sc_DElist.04.04.22.rda")
  
  names(FSGS_sc_DE) <- paste0(c("Coq2","Pdss2","GSE127235" ,"GSE139506",
                                "btbr", "cd2ap" ,"doxo", "nephr.D1", "nephr.D5", 
                                "GSE164273","GSE174013","GSE174102", "Nphs2","Wt1" ),"_sc")
  DE_list <- c( FSGS_MA_DE , FSGS_bulk_DE , FSGS_sc_DE)
  
  DEall_LFC <- Reduce( cbind , lapply( seq(DE_list), function(ii){
    # select 50 genes and column with LFC values
    X <- DE_list[[ii]]
    X$log2FoldChange[!is.finite(X$log2FoldChange)] <- max(X$log2FoldChange[is.finite(X$log2FoldChange)])
    X$log2FoldChange <- scale( X$log2FoldChange, center = F )
    X <- X[ match( DS_all.42$gene_symbol , 
                   X$Gene.Symbol ) , "log2FoldChange" ]
    return(X)
  }))
  colnames (DEall_LFC ) <- names(DE_list)
  rownames( DEall_LFC ) <- DS_all.42$gene_symbol
  
  # tame LFC outliers
  toPlot <- DEall_LFC
  toPlot[toPlot< -3] <- -3
  toPlot[toPlot> 3] <- 3
  # chose color scheme
  paletteLength <- 20
  myColor <-  colorspace::divergingx_hcl(palette = 'RdBu', 
                                         rev = T, n=20)
  #  save plots 
  library(viridis)
  setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise")
  heatmaply::heatmaply(toPlot, colors =myColor , na.value	=  "grey25", 
                       file = c("heatmaply_plot.pdf", "heatmaply_plot.png"),
                       width=1500, height= 1000)
  
}


#### 4. validate PDS in mouse spatial transcriptomics #### 
### read glom coordinates
TimGlom_coord <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/TimGlom_coord.rda" )

### GSM5713367 Spatial plot, with glom morphology and PDS
  {
    
    ### read glom coordinates
    TimGlom_coord <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/TimGlom_coord.rda" )
    
    ### read cell-type annotation
    ll <-list.files( pattern = "glom_annot" , full.names = T,
                     path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/BeadLocationsPodoTim/" )
    TimGlom_ctypes <- lapply( ll,  read.csv , row.names=1)
    names(TimGlom_ctypes) <- sub( "_.*","",basename(ll))
    
    
    ### Spatial plot: highlight gloms with circles colored by PDS
    library(ggforce)
    toPlotID <- "GSM5713367"
    toPlot <- TimGlom_coord[TimGlom_coord$geoID== "GSM5713367",] 
    
    ### read seurat object
    seur <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/filtSeur/GSM5713367_filtSeur.rda")
    # annotate all glomerular cells
    seur$glom.KNN <-  TimGlom_ctypes$GSM5713367$cell_type[ match( colnames(seur),
                                                                  TimGlom_ctypes$GSM5713367$barcode)]
    seur$glom.KNN[!is.na(seur$glom.KNN)] <- "glom" 
    
   
    # spatial plot of Slide'seq.v2 data from kidney-tissue slice
    pp<- SpatialDimPlot(seur, group.by = "glom.KNN",stroke = 0.1)+
      scale_fill_colorblind()+ theme( text = element_text( size=16))

    # plot circles coloured by PDS
    cc <- ggplot(toPlot, aes(x0=y, y0=x, r=radii, color=PDS, fill=PDS)) + 
      geom_circle(alpha=0.4) + coord_equal() + 
      ggeasy::easy_remove_axes()+theme_void()+  
      theme( text = element_text( size=16))+
      theme( legend.position = "bottom") +  scale_y_reverse() +
      scale_color_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
      scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))
    
      # combine on 2 panels, merge panels in the illustrator
   ppl  <- cowplot::plot_grid( cc, pp , rel_heights = c(1, 1), rows = 2)
      
   
   # save plots
   pdf(height = 15, width = 15, file = "Dimplot_glom.PDS_GSM5713367.pdf")
    ppl
   dev.off()
   
  }

### Suppl fig. D)
### Spearman correlation Heatmap of PDS and glom. attributes
  {
  toPlot <- TimGlom_coord[ ,c(
    "radii","PDS","Podocyte","gtypeBin","nFeatures","PodoFrct","normPodo")]
  corMatVal <- psych::corr.test( toPlot , method = "spearman", adjust = "fdr")
  corrplot::corrplot(corMatVal$r, p.mat  = corMatVal$p , 
                     col=rev(COL2('RdBu', 200)),method = 'color',
                     # addCoef.col = "black",
                     tl.cex = 2 ,
                     sig.level = c(0.01, 0.05),insig = 'label_sig')
  
  toPlot1 <- TimGlom_coord[TimGlom_coord$gtype!="BTBR-ob/ob",c(
    "radii","PDS","Podocyte","PodoFrct","nFeatures","normPodo")]
  corMatVal1 <- psych::corr.test( toPlot1 , method = "spearman", adjust = "fdr")
  corrplot::corrplot(corMatVal1$r, p.mat = corMatVal1$p , type = "upper",
                     col=rev(COL2('RdBu', 200)),method = 'color',
                     # addCoef.col = "black",  
                     tl.cex = 2, 
                     tl.pos = "full",title = "control",
                     sig.level = c(0.01, 0.05),insig = 'label_sig')
  
  toPlot2 <- TimGlom_coord[TimGlom_coord$gtype=="BTBR-ob/ob",c(
    "radii","PDS","Podocyte","PodoFrct","nFeatures","normPodo")]
  corMatVal2 <- psych::corr.test( toPlot2 , method = "spearman", adjust = "fdr")
  corrplot::corrplot( corMatVal2$r, p.mat = corMatVal2$p , 
                      # addCoef.col = "black", 
                      tl.cex = 2 , 
                      col=rev(COL2('RdBu', 200)),
                      method = 'color',
                      add = T , type = "lower", tl.pos = "n" ,
                      sig.level = c(0.01, 0.05),insig = 'label_sig',
                      title = "BTBR-ob/ob")
}

### Suppl fig. E) 
### regression and boxplots of PDS and podo N (normalised by the glom area)
  {
  library(rstatix)
  toPlot <- TimGlom_coord[, !colnames(TimGlom_coord)%in%
                            c("x","y","gtypeBin","EC","MC","PDS","oldPodo","PodoFrct",
                              "nFeatures")]
  toPlot <- reshape2::melt( toPlot ,
                            id.vars=c("glomID","geoID","gtype") )
  stat.test <- toPlot %>%
    group_by(variable) %>%
    wilcox_test(value ~ gtype) %>%
    adjust_pvalue( method = "fdr") 
  
  ## boxplot
  gg1 <- ggplot( data = toPlot , aes(x=gtype, y=value , color=gtype))+
    geom_boxplot(lwd=1.5)+ scale_color_colorblind()+
    # geom_dotplot(alpha=0.5)+ 
    theme_bw()+ 
    facet_grid( rows = vars(variable), scales="free")+
    theme( axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
           text=element_text( size=20),legend.position="bottom")+
    stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(350,40,0.003),size = 6)
  # facet_grid(gtype~gtype, scales = "free_y")
  
  ### scatterplots
  toPlot <- reshape2::melt( TimGlom_coord[,c("glomID","PDS","gtype","normPodo")], 
                            id.vars=c("glomID","gtype","PDS"))
  gg2 <- ggplot2::ggplot( data=toPlot, aes(y=PDS, x=value, color=gtype))+
    geom_point( size=5, alpha=0.5)+ scale_color_colorblind()+
    geom_smooth( method=lm ,se = FALSE, lwd=2) + theme_bw()+
    theme( text = element_text( size=20),
           legend.position = "none")+
    stat_cor(method="spearman",size=8 )+
    facet_grid( cols = vars(variable), scales = "free" )
  
  gg <- cowplot::plot_grid( plotlist = list(gg2,gg1), nrow = 1, rel_widths = c(2,1))
  
  pdf( width = 8, height = 6 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/Glom.metrics_scatter.bxplts.test.pdf")
  print(gg)
  dev.off()
  
}








