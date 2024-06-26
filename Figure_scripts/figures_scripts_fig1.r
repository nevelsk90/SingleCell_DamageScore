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








