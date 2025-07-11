###=====# Single-Cell Resolution of Cellular Damage Illuminates Disease Progression #=======###
###=====# supervised vs unsupervised method for disease signature #=======###

listSCSN.1K.sampl <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")

#### sling unsupervised sc/sn trajectory ####
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

### pseudo-time DS vs Proteinura in KFO snRNAseq datasets
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
  
  PDSvec <- grep("PDS",names(datt1), value = T)
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

