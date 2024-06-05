###=====# PDS publication Figure 2 code #=====###
# an input for each figure should be a data, 
# which requires no significant computations for plotting

#### load libraries, code and DS ####
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

# load necessary code
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# convert gene names to human orthologs
DS.HOMO <-  fun_homoTO.FROMmouse( gns = DS_all$gene_symbol, TO = F)
DS_all.HOMO <-  DS_all
DS_all.HOMO$HOMO <- unlist(DS.HOMO$HOMO)

setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure2")


#### 1. functional validation with mouse data #### 
### Main Fig A) PDS vs Proteinura in Mice
  {
  
}
#### 2. validate PDS in Humans RNAseq data #### 
### relate PDS and Clinical traits in podocytes from sc and snRNAseq KPMP data

# load precomputed data, filtered for CKD and AKI samples with at least 3 podocytes
  library(corrplot)
  toPlot.filt <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/HumanValidation/KPMP.CKD.AKI_scsnRNAseq_podo.PDS.aggregated.rda")

  ### Main Fig. B)
  ### box/dotplots of PDS in KPMP samples with different status
    {
    
    
  }

  ### Suppl Fig. B) 
  ### box/dotplots of PDS in KPMP samples with different status
    {         
    
    
  }  
  
  ### Main Fig. C)  
  ### spearman correlation heatmap for PDS vs selected clinical trait
    {
      toPlot.filt.red <- toPlot.filt[,c("PDS" ,"Age", "Sex", "Race", "Diabetes.Duration..Years.",
                                        "Hypertension.Duration..Years.", "eGFR.CrCys",
                                        "ACR.mean")]
      cor.test(toPlot.filt.red$PDS, toPlot.filt.red$ACR.mean, method = "spearman")
      
      corMat.red <- Reduce( cbind.data.frame, lapply( 1:ncol(toPlot.filt.red),  function(ii){
        X <- toPlot.filt.red[,ii]
        if( !is.numeric(X)) { 
          X <-  as.numeric( as.factor(X))
          return(X)} else return(X)
      }))
      rownames(corMat.red) <- toPlot.filt.red$partIDsnsc
      colnames(corMat.red) <- c("PDS" ,"Age", "Sex", "Race", "Diabetes Duration",
                                "Hypertension Duration", "eGFR Cr.Cys",
                                "ACR mg/g")
      corMat.red <- as.data.frame( corMat.red[,apply(corMat.red, 2, sd,na.rm=T)>0] )
      
      corMat.test <-  psych::corr.test( corMat.red,  method = "spearman", 
                                        use = "pairwise.complete.obs" , adjust = "none")
      # make a plot
      corrplot::corrplot(corMat.test$r,order="hclust", col= rev(COL2('RdBu', 200)),
                         method = 'color',p.mat =corMat.test$p,
                         sig.level = 0.1,
                         title = "KPMP sn.sc.RNAseq",insig = 'label_sig',
                         tl.cex = 1.5, cl.cex = 1, pch.cex = 2)
      
    }

  ### Suppl Fig. C) 
  ### spearman correlation heatmap for PDS vs all clinical traits
    {         
    
    corMat <- Reduce( cbind.data.frame, lapply( 1:ncol(toPlot.filt),  function(ii){
      X <- toPlot.filt[,ii]
      if( !is.numeric(X)) { 
        X <-  as.numeric( as.factor(X))
        return(X)} else return(X)
    }))
    rownames(corMat) <- toPlot.filt$partIDsnsc
    colnames(corMat) <- colnames(toPlot.filt)
    corMat <- as.data.frame( corMat[,apply(corMat, 2, sd,na.rm=T)>0] )
    
    corMatVal.sc.test <-  psych::corr.test( corMat[ corMat$seqType==1 , 
                                                    !colnames(corMat)%in% c("partIDsnsc","partID","Participant.ID","seqType")],
                                            # "seqType","Tissue.Type")],
                                            method = "spearman", 
                                            use = "pairwise.complete.obs" )
    corMatVal.sc <- corMatVal.sc.test$r
    corMatVal.sc[is.na(corMatVal.sc)]<- 0
    
    corMatVal.sn.test <-  psych::corr.test( corMat[ corMat$seqType==2 , 
                                                    !colnames(corMat)%in% c("partIDsnsc","partID","Participant.ID","seqType")],
                                            # "seqType","Tissue.Type")], 
                                            method = "spearman", 
                                            use = "pairwise.complete.obs" )
    corMatVal.sn <- corMatVal.sn.test$r
    corMatVal.sn[is.na(corMatVal.sn)]<- 0
    
    corMat.test <-  psych::corr.test( corMat[ ,!colnames(corMat)%in% c("partIDsnsc","partID",
                                                                       "Participant.ID")], 
                                      method = "spearman", adjust = "none",
                                      use = "pairwise.complete.obs" )
    corMat.test$r[is.na(corMat.test$r)]<- 0
    
    # make a plot
    corrplot::corrplot(corMatVal.sc, type = "upper", 
                       tl.pos = "full", 
                       order="original",
                       method = 'color',  col=rev(COL2('RdBu', 200)),
                       sig.level = 0.05,insig = 'label_sig',
                       p.mat =  corMatVal.sc.test$p, 
                       title = "KPMP scRNAseq",
                       tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
    corrplot::corrplot(corMatVal.sn, type = "lower", tl.pos = "n", 
                       order="original",  col=rev(COL2('RdBu', 200)),
                       method = 'color',insig = 'label_sig',
                       sig.level = 0.05,p.mat =  corMatVal.sn.test$p, 
                       add = T,   title = "KPMP snRNAseq",
                       tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
    
    corrplot::corrplot(corMat.test$r,order="hclust",
                       col=rev(COL2('RdBu', 200)),
                       method = 'color',p.mat =corMat.test$p,
                       sig.level = 0.05,
                       # title = "KPMP sn.sc.RNAseq",
                       insig = 'label_sig',
                       tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
    
    
  }

  


#### 3. validate PDS in Human spatial transcriptomics #### 
### visualise PDS in the highQ histo-image from KPMP patient 29.10282 
  
  # load Seurat object with annotated gloms
  kpmp.sptl_29.10282.seur <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/kpmp.sptl_29.10282.seur.rda")
  
  ### Main fig. E) 
  ### Spatial plot, PDS visualised 
  {
    kpmp.sptl_29.10282.seur$glom.PDS <- ifelse(kpmp.sptl_29.10282.seur$glom.annot=="",
                                               NA, kpmp.sptl_29.10282.seur$PDS.42)
    

    gg1<- SpatialFeaturePlot( kpmp.sptl_29.10282.seur, features = "glom.PDS",
                              pt.size.factor = 2 ) + 
      scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
      theme(legend.position = "bottom", text = element_text(size=20)) +
      scale_x_reverse()
    
    gg2<- SpatialDimPlot( kpmp.sptl_29.10282.seur, group.by = "glom.annot",
                          pt.size.factor = 2) + 
      theme(legend.position = "bottom", text = element_text(size=20)) +
      guides(fill = guide_legend(override.aes = list(size=8)))+
      scale_x_reverse()
    gg3<- SpatialFeaturePlot( kpmp.sptl_29.10282.seur, features =  "WT1NPHS12",
                              pt.size.factor = 2) + 
      scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
      theme(legend.position = "bottom", text = element_text(size=20))+
      scale_x_reverse()
    # scale_fill_viridis_b(option = "A",direction = -1,n.breaks=7,)
    
    cowplot::plot_grid(plotlist = list( gg2 , gg1 , gg3), ncol = 1)
    
        }
  
  
  ### Suppl. fig D)
  ### summary boxplot: PDS by glom type
  ggplot(kpmp.sptl_29.10282.seur@meta.data, 
         aes(x=glomType, y=PDS.42, fill=glomType))+ 
    geom_boxplot()+theme_bw()+
    geom_point(size=3,alpha = 0.3)+
    theme( text=element_text(size = 20),
           axis.text.x = element_blank())+
    annotate(geom="text",label= "AUCell.thrsh=0.3", x=2, y=0)