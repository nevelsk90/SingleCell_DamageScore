###=====# PDS publication Figure 2 code #=====###
# an input for each figure should be a data, 
# which requires no significant computations for plotting

# release memory
gc()
mallinfo::malloc.trim()
gc()

#### load libraries, code and DS ####
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", 
            "/media/tim_nevelsk/WD_tim/SOFT/R"))
options( connectionObserver = NULL )
# library( org.Mm.eg.db )
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
library( RColorBrewer )

# setwd
setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure2")
inputdir <-  "/media/tim_nevelsk/WD_tim/PROJECTS/WRITING/PDS_manuscript/Figure_input" 

# load necessary code
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# convert gene names to human orthologs
DS.HOMO <-  fun_homoTO.FROMmouse( gns = DS_all$gene_symbol, TO = F )
DS_all.HOMO <-  DS_all
DS_all.HOMO$HOMO <- unlist(DS.HOMO$HOMO)
# load scsn podo genes
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

#### 1. functional validation with mouse data #### 
### PDS vs Proteinura in KFO snRNAseq datasets
  {
    # read expression data
    listSCSN <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_22.12.23.rda" )
    listSCSN.1K.sampl <- lapply( seq( listSCSN ), function(ii)
      {
      print(ii)
      # if( names(listSCSN)[ii]=="Nphs2") {
      #   newSeu <- subset( listSCSN[[ii]] , subset= sample %in% c(
      #     "140739", "140738", "140740", "140741", "139919",
      #     "139917", "139921", "139913", "139915", "139911" ))
      # } else 
      newSeu <- listSCSN[[ii]]
      
      
      # balance samples
      Idents(newSeu)<- newSeu$sample
      newSeu <- subset( newSeu , downsample=1000 )
      
      return(newSeu)
    })
    names(listSCSN.1K.sampl ) <- names(listSCSN)
    saveRDS(listSCSN.1K.sampl , "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
    
    
    # laod annotation
    annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
    annot_tab$group <- paste(annot_tab$Genotype,annot_tab$Age_weeks,sep = "_")
    
    # combine metadata from 3 experiments 
    datt <- Reduce( rbind , lapply( seq(listSCSN.1K.sampl), function(ii) {
      XX <- listSCSN.1K.sampl[[ii]]@meta.data 
      XX <- XX[,c( "group","sample","gtype", "PDS")]
      XX$dataSet <- names(listSCSN.1K.sampl)[ii]
      return(XX)
    }))
    
    # treat carefully 21 week Pdss2 samples since they have only per group measurements
    datt1 <- datt[ datt$sample %in% c( "146985", "146986", "143485" , "143486") ,]
    # aggregate  Pdss2 samples
    aggPDS1 <- aggregate( .~group , FUN = mean , 
                          data= datt1[ , c("group", "PDS")] )
    aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
    aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
    colnames(aggPDS1)[1] <- "sample"
    
    # the rest of samples
    datt2 <- datt[ !(datt$sample %in%  c("146985", "146986", "143485" , "143486")), ]
    aggPDS2 <- aggregate( .~sample , FUN=mean,
                          data= datt2[ , c( "sample", "PDS" )] ) 
    # combine all samples
    aggPDS <- rbind(aggPDS2, aggPDS1)
    
    aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ 
      match( sub("SID","" ,aggPDS$sample) , annot_tab$CCG_Sample_ID)]
    aggPDS$group <- annot_tab$group[ 
      match( sub("SID","" ,aggPDS$sample ), annot_tab$CCG_Sample_ID)]
    aggPDS$gtype <- as.factor(annot_tab$Genotype[ 
      match( sub("SID","" ,aggPDS$sample ) , annot_tab$CCG_Sample_ID)])
    aggPDS$gtypeDE <- ifelse( aggPDS$gtype=="wt" , "control", "experimental")
    
    # PDSvec <- grep("PDS",names(datt1), value = T)
    aggPDS <- aggPDS[ !is.na(aggPDS$AlbCrRatio),]
    aggPDS$dataSet <- datt$dataSet[ match( aggPDS$sample, datt$sample )]
    

    # plot lm and correlation
    gg1 <- ggplot2::ggplot( data = aggPDS, aes( 
      x = aggPDS[,PDSvec[ii]], y = log10(AlbCrRatio) ) ) +
      geom_point(  size=6, aes(col=dataSet, shape=gtypeDE)) +
      theme_bw() +  theme( text = element_text(size = 22)) + 
      geom_smooth( method='lm', color="black", se = FALSE) + 
      ggtitle( "snRNAseq KFO data" ) + xlab("aggregated PDS")+
      stat_cor( size=7, method = "spearman" ) 
    # geom_text(aes(label = sample  ), size=6, position = "dodge")
    
    
    pdf(height = 4, width = 7, file = "PDS42vsAlbCr_KFO.scatter.pdf")
      gg1
    dev.off()
    png(height = 300, width = 600, file = "PDS42vsAlbCr_KFO.scatter.png")
      gg1
    dev.off() 
  }

### Suppl fig. pseudo-time DS vs Proteinura in KFO snRNAseq datasets
  {
  # read ptime signature
    pseudo <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1/Supl.Fig1/listSCSN.1K_sling.ptimeDE.rda")
    pseudoDS <- pseudo[pseudo$Gene.Symbol %in% allPodoGenes, ]
    pseudoDS <- data.frame( gene_symbol=pseudoDS$Gene.Symbol,
                            mean_rank = seq(pseudoDS$Gene.Symbol), 
                            direction_foldchange = ifelse(
                              pseudoDS$log2FoldChange<0 , -1, 1) )
    # convert gene names to human orthologs
    pseudoDS.HOMO <-  pseudoDS[1:300,]
    pseudoDS.HOMO$HOMO <- unlist(fun_homoTO.FROMmouse( 
      gns = pseudoDS.HOMO$gene_symbol, TO = F)$HOMO)                       
    
  listSCSN.1K.ptime <- lapply( seq(listSCSN.1K.sampl), function(ii)
    {
      print(ii)
      newSeu <- listSCSN.1K.sampl[[ii]]
      
      
      ## exclude genes non.expressed in more than certain percentage of cells
      expr <- newSeu@assays$RNA@counts
      expr <- expr[ rowSums( round(expr) > 0 ) > 0 , ]
      
      # calculate damage signatures
      newSeu@meta.data$PDS <- DS_calc.func( exprMatrices = expr , 
                                            DSignature = pseudoDS , 
                                            ntop = 42 , wghtd = T,
                                            ceilThrsh =  0.05 )
      
      # adjust ceilThrsh based on a total number of genes in the matrix!
      # the top should include ~ 1K genes
      
      return(newSeu)
    })
  names(listSCSN.1K.ptime) <- names(listSCSN.1K.sampl)
  

  # combine metadata from 3 experiments 
  datt <- Reduce( rbind , lapply( seq(listSCSN.1K.ptime), function(ii) {
    XX <- listSCSN.1K.ptime[[ii]]@meta.data 
    XX <- XX[,c( "group","sample","gtype", "PDS")]
    XX$dataSet <- names(listSCSN.1K.ptime)[ii]
    return(XX)
  }))
  
  # treat carefully 21 week Pdss2 samples since they have only per group measurements
  datt1 <- datt[ datt$sample %in% c( "146985", "146986", "143485" , "143486") ,]
  # aggregate  Pdss2 samples
  aggPDS1 <- aggregate( .~group , FUN = mean , 
                        data= datt1[ , c("group", "PDS")] )
  aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
  aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
  colnames(aggPDS1)[1] <- "sample"
  
  # the rest of samples
  datt2 <- datt[ !(datt$sample %in%  c("146985", "146986", "143485" , "143486")), ]
  aggPDS2 <- aggregate( .~sample , FUN=mean,
                        data= datt2[ , c( "sample", "PDS" )] ) 
  # combine all samples
  aggPDS <- rbind(aggPDS2, aggPDS1)
  
  aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ 
    match( sub("SID","" ,aggPDS$sample) , annot_tab$CCG_Sample_ID)]
  aggPDS$group <- annot_tab$group[ 
    match( sub("SID","" ,aggPDS$sample ), annot_tab$CCG_Sample_ID)]
  aggPDS$gtype <- as.factor(annot_tab$Genotype[ 
    match( sub("SID","" ,aggPDS$sample ) , annot_tab$CCG_Sample_ID)])
  aggPDS$gtypeDE <- ifelse( aggPDS$gtype=="wt" , "control", "experimental")
  
  # PDSvec <- grep("PDS",names(datt1), value = T)
  aggPDS <- aggPDS[ !is.na(aggPDS$AlbCrRatio),]
  aggPDS$dataSet <- datt$dataSet[ match( aggPDS$sample, datt$sample )]
  
  
  # plot lm and correlation
  gg1 <- ggplot2::ggplot( data = aggPDS, aes( 
    x = aggPDS[,PDSvec[ii]], y = log10(AlbCrRatio) ) ) +
    geom_point(  size=6, aes(col=dataSet, shape=gtypeDE)) +
    theme_bw() +  theme( text = element_text(size = 22)) + 
    geom_smooth( method='lm', color="black", se = FALSE) + 
    ggtitle( "snRNAseq KFO data" ) + xlab("aggregated PDS")+
    stat_cor( size=7, method = "spearman" ) 
  # geom_text(aes(label = sample  ), size=6, position = "dodge")
  
  # density plots 
  gg2 <- ggplot2::ggplot( data = datt, aes( x=datt[, "PDS"], color=gtype)) +
    geom_density(size=1.5) + ggtitle( PDSvec[ii])+ 
    theme_bw() +  theme( text = element_text(size = 22)) 
  
  # dotplot for samples of Nphs2mut
  aggPDS.nphs2 <- datt[datt$sample%in% listSCSN.1K.ptime$Nphs2$sample,]
  aggPDS.nphs2 <- aggregate( .~sample , FUN=mean,
             data= aggPDS.nphs2[ , c( "sample", "PDS" )] ) 
  aggPDS.nphs2$group <- annot_tab$group[ 
    match( sub("SID","" ,aggPDS.nphs2$sample ), annot_tab$CCG_Sample_ID)]
  gg3 <- ggplot2::ggplot( data = aggPDS.nphs2, 
                          aes( y=aggPDS.nphs2[,"PDS"], 
                               x=group, 
                               color=group)) +
    geom_jitter(size=1.5) + ggtitle( PDSvec[ii])+ 
    theme_bw() +  theme( text = element_text(size = 22)) +
    geom_label(aes(label=sample))
  
  # density plots for nphs2
  datt.nphs2 <- datt[datt$sample%in% listSCSN.1K.ptime$Nphs2$sample,]
  gg4 <- ggplot2::ggplot( data = datt.nphs2, 
                          aes( x=datt.nphs2[, "PDS" ], 
                               color=group)) +
    geom_density(size=1.5) + ggtitle( PDSvec[ii])+ 
    theme_bw() +  theme( text = element_text(size = 22)) 
  
  ggl <- cowplot::plot_grid(plotlist = list( gg1, gg2, gg3, gg4 ) , ncol=2)
  ggl
  
  pdf(height = 10, width = 18, file = "Supl.Fig2/ptimeDSvsAlbCr_KFO.scatter.pdf")
  ggl
  dev.off()

}

### Suppl fig. PDS vs Proteinura in bulk public datasets
  {
  stud <- c( "GSE117571", "GSE108629", "GSE17709" ,
             "GSE117987", "GSE126217" ,"GSE154955",
             "GSE112116", "GSE131266", "GSE134327", 
             "GSE110092" ,"GSE77717" , "KFO.Wt1"   )
  
  # read stage data
  annot_bulkStage <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/annot_bulkStage.rda")
  # read proteinuria data
  mgmg.Bulk <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/mgmg.Bulk.rda")
  # read bulk expression data PDS
  explist <- readRDS(  file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.sc.rda")
  # select studies
  blkData <- explist[ names(explist) %in% stud]
  
  # calculate PDS
  PDSsize.bulk <-  lapply( seq(blkData) , function(ii , ceilThrsh = 0.05  )
  {
    DS_calc.func( exprMatrices = blkData[[ii]], 
                  ceilThrsh = ceilThrsh , 
                  wghtd = T, useOrder = "mean_rank",
                  DSignature= DS_all , ntop = 42 )
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
  pdf(height = 4, width = 8, file = "Supl.Fig2/PDS.42vsAlbCr_bulk.pdf")
  print(gg)
  dev.off()
  # save the plot
  png(height = 400, width = 800, file = "Supl.Fig2/PDS.42vsAlbCr_bulk.png")
  print(gg)
  dev.off()
  
}


#### 3. relate PDS and Clinical traits in podocytes #### 
  ### snRNAseq KPMP data, correlation heatmap for PDS vs selected clinical trait
    {
      # load precomputed data, filtered for CKD and AKI samples with at least 3 podocytes
      library(corrplot)
      toPlot.filt <-readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/HumanValidation/KPMP.CKD.AKI_scsnRNAseq_podo.PDS.aggregated.rda")
      
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
      corrplot::corrplot(corMat.test$r,order="hclust", 
                         # type = "lower",
                         col= rev(COL2('RdBu', 200)),
                         tl.col = "black",
                         tl.pos = "l", 
                         method = 'color',
                         p.mat = corMat.test$p ,
                         sig.level = 0.05,
                         insig = 'label_sig',
                         tl.cex = 1.5, cl.cex = 1.5, pch.cex = 2)
      
    }

  ### Suppl Fig. correlation heatmap for PDS vs all clinical traits
    {         
    
      corMat <- Reduce( cbind.data.frame, 
                        lapply( 1:ncol(toPlot.filt),  function(ii){
                          X <- toPlot.filt[,ii]
                          if( !is.numeric(X)) { 
                            X <-  as.numeric( as.factor(X))
                            return(X)} else return(X)
                        }))
      rownames(corMat) <- toPlot.filt$partIDsnsc
      colnames(corMat) <- colnames(toPlot.filt)
      corMat <- as.data.frame( corMat[,apply(corMat, 2, sd,na.rm=T)>0] )
      
      corMat.test <-  psych::corr.test( corMat[ ,!colnames(corMat)%in% c("partIDsnsc","partID",
                                                                         "Participant.ID")], 
                                        method = "spearman", adjust = "none",
                                        use = "pairwise.complete.obs" )
      corMat.test$r[is.na(corMat.test$r)]<- 0
      
      # plot
      corrplot::corrplot(corMat.test$r,mar = c(0,0,0,0),
                         order="hclust",
                         col=rev(COL2('RdBu', 200)),
                         method = 'color', 
                         tl.pos = "l", tl.col = "black",
                         p.mat =corMat.test$p,
                         sig.level = 0.05,
                         # title = "KPMP sn.sc.RNAseq",
                         insig = 'label_sig',
                         tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
      
      
    
    
  }

  ### Suppl. fig. GSE176465 urine samples, PDS vs UPCR
    {
    seu.filt <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_seurat.rda" )
    
    Seurat::FeaturePlot( seu.filt ,
                         order=T , min.cutoff = 0,
                         features =  "WT1" ,
                         label = T, cols=c("lightgray","red"))
    
    seu.podo <- subset( seu.filt , subset= seurat_clusters == 9 )
    expMat <- seu.podo@assays$RNA@layers$counts
    rownames(expMat)<- rownames(seu.podo@assays$RNA)
    colnames(expMat)<- colnames(seu.podo@assays$RNA)
    expMat <- expMat[rowSums(expMat>0)> 0 ,] # >0 for scRNAseq, >10 for snRNAseq
    
    ### calcualte PDS, ceilThrsh must be >0.05 for HUMAN scRNAseq and snRNAseq
    set.seed(42)
    seu.podo$PDS  <- DS_calc.func( exprMatrices = expMat ,
                                   ceilThrsh = 0.1 ,
                                   DSignature = DS_all.HOMO,
                                   wghtd = T,
                                   progStat = T,
                                   geneIDname = "HOMO",
                                   ntop = 42)
    
    # aggregate
    seu.podo_pbulkM <- aggregate( .~patient+ biopsy_diagnosis+GSM,
                                  data = seu.podo@meta.data[
                                    , c("PDS","patient","biopsy_diagnosis","GSM")],
                                  FUN = mean )
    
    gg4 <- ggplot2::ggplot(  seu.podo_pbulkM  , aes(
      y= PDS, x=biopsy_diagnosis, color=patient)) +
      scale_color_colorblind() +
      geom_jitter(width = 0.1, size=6) +
      ggtitle("GSE176465 pseudo-bulk")+
      theme_bw() + theme(text = element_text(size = 24) ,
                         legend.position = "none")
    # stat_summary(fun.y= median, fun.ymin=median, fun.ymax=median, geom = "crossbar", width = .2, color = "red")
    # stat_compare_means( size = 8)
    gg4
    
    # correlate with UPCR
    seu.podo_pbulk.sbj <- aggregate( .~patient + biopsy_diagnosis,
                                     data = seu.podo@meta.data[
                                       , c("PDS","patient","biopsy_diagnosis")],
                                     FUN = mean )
    upcr.tab <- read.table(header = T, sep = "\t","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/patient_data.csv")
    seu.podo_pbulk.sbj$UPCR <- as.numeric( upcr.tab$UPCR..at.urine.collection.[
      match( seu.podo_pbulk.sbj$patient , upcr.tab$Study.participants) ])
    
    seu.podo_pbulk.sbj$UPCR.highest <- as.numeric( upcr.tab$UPCR..highest.recorded.[
      match( seu.podo_pbulk.sbj$patient , upcr.tab$Study.participants ) ])
    seu.podo_pbulk.sbj$
      
      cor.test( as.numeric(seu.podo_pbulk.sbj$UPCR),
                as.numeric(seu.podo_pbulk.sbj$PDS), method="spearman")
    
    ggplot( seu.podo_pbulk.sbj , aes( x = log10( UPCR ),
                                      y = PDS ))+
      geom_point(size=6, aes( color=biopsy_diagnosis))+theme_bw()+ 
      theme( text = element_text(size = 24), legend.position = "bottom")+ggtitle("GSE176465 pseudo-bulk")+
      geom_smooth( method = "lm", color="black",se = F)+
      stat_cor( method = "spearman", size=8)
    
  }

#### 4. validate PDS in Human spatial transcriptomics #### 
### visualise PDS in the highQ histo-image from KPMP patient 29.10282 
  
  # load Seurat object with annotated gloms
  kpmp.sptl_29.10282.seur <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/kpmp.sptl_29.10282.seur.rda")
  
  ### Spatial plot, PDS visualised 
    {
  kpmp.sptl_29.10282.seur$glom.PDS <- ifelse( is.na(kpmp.sptl_29.10282.seur$glom.annot),
                                             NA, kpmp.sptl_29.10282.seur$PDS.42)
  
      
      ### download barcode with coordinates from loupe browser
      coord.df <- read.csv(header = T, row.names=1, sep = ",","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/lowres.KPMP/Spatial-Projection.csv")
      coord.df <- coord.df[colnames(kpmp.sptl_29.10282.seur),]
      colnames(coord.df) <- c("imagerow" ,"imagecol")
      kpmp.sptl_29.10282.seur@images$image =  new(
        Class = 'SlideSeq',
        assay = "Spatial",
        key = "image_",
        coordinates = coord.df
      )

  # 
  gg2 <-kpmp.sptl_29.10282.seur@images$image@coordinates %>%
    bind_cols(kpmp.sptl_29.10282.seur@meta.data) %>%
    # filter(tissue==1) %>%
    mutate(x=imagecol,y=imagerow) %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(col=glom.PDS)) +
    scale_color_distiller(palette = "Spectral")+
    scale_x_reverse()+ scale_y_reverse()+ coord_fixed()+ theme_void() +
    theme(legend.position = "bottom", text = element_text(size=20))
  

  gg3 <-kpmp.sptl_29.10282.seur@images$image@coordinates %>%
    bind_cols(kpmp.sptl_29.10282.seur@meta.data) %>%
    mutate(x=imagecol,y=imagerow) %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(col=glom.annot)) +
    # scale_color_distiller(palette = "Spectral")+
    scale_x_reverse()+ scale_y_reverse()+ coord_fixed()+ theme_void()+
    theme(legend.position = "bottom", text = element_text(size=20))
  
  
  # gg4<- SpatialFeaturePlot( kpmp.sptl_29.10282.seur, 
  #                           features =  "WT1NPHS12",
  #                           pt.size.factor = 10) + 
  #   scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
  #   theme(legend.position = "bottom", text = element_text(size=20))+
  #   scale_x_reverse()
  
  gg4 <-kpmp.sptl_29.10282.seur@images$image@coordinates %>%
    bind_cols(kpmp.sptl_29.10282.seur@meta.data) %>%
    mutate(x=imagecol,y=imagerow) %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(col=WT1NPHS12)) +
    scale_color_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
    scale_x_reverse()+ scale_y_reverse()+ coord_fixed()+ theme_void()+
    theme(legend.position = "bottom", text = element_text(size=20))
  
  # scale_fill_viridis_b(option = "A",direction = -1,n.breaks=7,)
  
  ggl <- cowplot::plot_grid(plotlist = list( gg3 , gg2 , gg4), ncol = 1)
  
  pdf( height = 10, width = 10, file="KPMP_29.10282_Spatial.PDS.pdf")
    ggl
  dev.off()
  }
  
  
  ### Suppl. fig  summary boxplot: PDS by glom type
  # boxplots
  gg0 <- ggplot(kpmp.sptl_29.10282.seur@meta.data, 
         aes(x=glomType, y=PDS.42, fill=glomType))+ 
    geom_boxplot()+theme_bw()+
    geom_point(size=3,alpha = 0.3)+
    theme( text=element_text(size = 20),
           axis.text.x = element_blank())+
    annotate(geom="text",label= "AUCell.thrsh=0.3", x=2, y=0)
  
  pdf( height = 10, width = 10, file="Supl.Fig2/KPMP_29.10282_gloms.PDS42_boxplots.pdf")
  gg0
  dev.off()
