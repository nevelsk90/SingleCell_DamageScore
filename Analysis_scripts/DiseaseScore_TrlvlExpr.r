
#### viz.  expression of ind. genes in single-cell podocyte datasets ####
## load data 
  listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
  allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")
  
### one gene 
  Ggene <- "Dpp4"
  Ggene <- "Pkm"
    # Ggene <- c("Clock","Bmal1","Arntl","Bmal2","Arntl2","Npas2",
  #            "Nr1d1" , "Nr1d2" , "Rora", "Rorb", "Rorc" )
  GeneExprSCSN <- Reduce( rbind, lapply( seq(listSCSN.1K.sampl), function(ii,
                                                                          ggene=Ggene)
  {
    print( ii )
    
    datt <-  listSCSN.1K.sampl[[ii]]
    
    datt1 <- data.frame( t(datt@assays$RNA@data[ 
      rownames(datt@assays$RNA@data) %in% ggene ,,drop=F] ) )
    
    
    datt1$gtypeDE <- datt@meta.data$gtypeDE
    datt1$dataSet <- names(listSCSN.1K.sampl)[ii]
    return(datt1)
  }) )
  
  toPlot <- GeneExprSCSN[!GeneExprSCSN$dataSet%in% c("doxo","nephr.D5","Lmx1b"),]
  # toPlot <- GeneExprSCSN[GeneExprSCSN$dataSet%in% c("Nphs2","Wt1","Pdss2"),]
  # toPlot <- GeneExprSCSN[!GeneExprSCSN$dataSet%in% c("Nphs2","Wt1","Pdss2","doxo","nephr.D5","Lmx1b"),]
  # 
  toPlot.m <- reshape2::melt( toPlot)
  toPlot.m$variable <- factor(toPlot.m$variable , levels = Ggene)
  
  
  gg0 <- ggplot(toPlot.m, aes(y= value,
                              x=dataSet, fill=gtypeDE))+
    geom_violin()+theme_bw()+
    theme(text = element_text(size = 20))+
    facet_wrap(vars(variable))+
    stat_compare_means(label = "p.format")+
    stat_summary(fun = mean, geom = "crossbar",
                 colour = "red",
                 position = position_dodge(width = 0.9))+
    scale_fill_colorblind()+ ylab("Dpp4 expr")
  
  pdf( width = 12, height = 6, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Revision1/Pkm.exprSCSN_trnscrpt.lvl.pdf")
  print(gg0)
  dev.off()
  
  

