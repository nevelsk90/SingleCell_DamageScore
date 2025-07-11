
options( connectionObserver = NULL )
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", 
            "/media/tim_nevelsk/WD_tim/SOFT/R"))

mallinfo::malloc.trim()
gc()


library( GSEABase )
library( biomaRt )
library( ggplot2 )
library( ggthemes)
library( ggpubr)
library( DEP )

source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")

### load kegg paths
keggPath_pathlist <- readRDS( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/kegg_pathsList_gName.07.09.23.rda")
# podocyte genes
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")


#### process proteomics data ####
filt.thrsh <- 0.49
###   read in Biognosys Spectronaut data, PXD04741
{
  library(tidyverse)
  
  BSproteome_data <- readr::read_tsv("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/PROTEOME/PXD047417_Nephrin-Ab-Proteome_Report_Protein_Quant.tsv")
  
  # Create a new column that uniquely identifies each sample. 
  # This could be built from R.FileName or a combination of R.Condition and R.Replicate.
  proteome_data <- BSproteome_data %>%
    mutate(sample = paste(R.Condition, R.Replicate, sep = "_"))
  
  # Pivot the data so each sample becomes a column.
  # We assume that each row corresponds to one protein measurement.
  wide_data <- proteome_data %>%
    select(PG.ProteinAccessions, PG.Genes, sample, PG.Quantity) %>%
    pivot_wider(names_from = sample, values_from = PG.Quantity, values_fn=min)
  wide_data <- wide_data %>%
    mutate(across(where(is.numeric), ~ ifelse(!is.finite(.), NA, .)))
  # View the wide-formatted data
   
  # Assuming that sample names follow the pattern "Condition_Replicate"
  # First, extract unique sample names from the long-format data.
  sample_names <- unique(proteome_data$sample)
  
  # Create a metadata data frame
  metadata <- data.frame(
    label = sample_names,
    condition = sapply(strsplit(sample_names, "_"), `[`, 1),
    replicate = sapply(strsplit(sample_names, "_"), `[`, 2),
    stringsAsFactors = FALSE
  )
  
  print(metadata) 
  
  # ## split rows 
  # proteome_unique <- wide_data %>%
  #   separate_rows( PG.ProteinAccessions , sep = ";")
  ## make unique
  proteome_unique <- DEP::make_unique( wide_data , names="PG.Genes", 
                                       ids="PG.ProteinAccessions", delim = ";" )

  
  proteome_unique$name <- toupper( proteome_unique$name )
  
  ## 
  MRpodo_proteome_se <- DEP::make_se( proteome_unique , 3:14, metadata )


    # Filter for proteins that are identified in all replicates of at least one condition
  MRpodo_proteome_se <- filter_missval(MRpodo_proteome_se, thr = round(
    min(summary( as.factor(metadata$condition)))*filt.thrsh))
    
    # Plot a barplot of the protein identification overlap between samples
    print( plot_coverage(MRpodo_proteome_se) )
    
    # Normalize the data
    MRpodo_proteome_se <- normalize_vsn(MRpodo_proteome_se)
    meanSdPlot( MRpodo_proteome_se )
    
    # # Visualize normalization by boxplots for all samples before and after normalization
    # plot_normalization(Nphs2_proteome_filt, Nphs2_proteome_Norm)
    
    # imputation
    plot_detect( MRpodo_proteome_se )
    MRpodo_proteome_se  <- impute( MRpodo_proteome_se, fun = "MLE")
    
  
  
}

### read human DIA MS PXD054062
{
  metadata2 <- readxl::read_excel(c, sheet = 1)
  metadata2 <- data.frame(
    label = metadata2$`LC-MS run`,
    condition = metadata2$`sample type`,
    replicate = paste( sapply(strsplit(metadata2$`LC-MS run`, "_"), `[`, 2),
                       metadata2$`Experiment #` , sep = "."),
    batch= metadata2$`Experiment #`,
    sample = metadata2$`LC-MS run`,
    stringsAsFactors = FALSE
  )
  metadata2$condition <- ifelse( metadata2$condition== "Healthy control", "control",
                                 "PLA2R.MN" )

  ## read the data
  proteome_data1 <- readxl::read_excel("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/PROTEOME/PXD054062_Report_unique_genes_NEFR1028_DIA.xlsx", sheet = 1)
  colnames(proteome_data1) <- gsub( "D:|\\\\NEFR\\\\|.raw|\\\\", "", colnames(proteome_data1))
  
  proteome_data2 <- readxl::read_excel("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/PROTEOME/PXD054062_Report_unique_genes_NEFR849_DIA.xlsx", sheet = 1)
  colnames(proteome_data2) <- sub( "_2022.*", "", 
                                   gsub( "X:\\\\Didier\\\\Data\\\\NEFR\\\\|.raw|_2022.*|_DIA", "", 
                                         colnames(proteome_data2)) )
 
  ## make unique
  proteome_unique1 <- DEP::make_unique( proteome_data1 , names="Genes", 
                                       ids="Genes", delim = ";")
  proteome_unique2 <- DEP::make_unique( proteome_data2 , names="Genes", 
                                        ids="Genes", delim = ";")
  ## 
  MRpodo_proteome_se1 <- DEP::make_se( proteome_unique1 ,
                                       columns = 2:32, 
                                       expdesign =  metadata2[metadata2$batch==5,] )
  MRpodo_proteome_se2 <- DEP::make_se( proteome_unique2 ,
                                       columns = 2:29, 
                                       expdesign =  metadata2[metadata2$batch==4,] )
  
  # Filter for proteins that are identified in all replicates of at least one condition
  MRpodo_proteome_se1 <- filter_missval(MRpodo_proteome_se1, 
                                        thr = round(min(summary(
                                          as.factor(metadata2[metadata2$batch==5,
                                                              "condition"])))*filt.thrsh))
  MRpodo_proteome_se2 <- filter_missval(MRpodo_proteome_se2, 
                                        thr = round(min(summary(
                                          as.factor(metadata2[metadata2$batch==4,
                                                              "condition"])))*filt.thrsh))
  # # Plot a barplot of the protein identification overlap between samples
  # print( plot_coverage(MRpodo_proteome_se1) )
  # print( plot_coverage(MRpodo_proteome_se2) )
  
  # Normalize the data
  MRpodo_proteome_se1 <- normalize_vsn(MRpodo_proteome_se1)
  MRpodo_proteome_se2 <- normalize_vsn(MRpodo_proteome_se2)
  
  # meanSdPlot( MRpodo_proteome_se1 )
  # meanSdPlot( MRpodo_proteome_se2 )
  # 
  # # Visualize normalization by boxplots for all samples before and after normalization
  # plot_normalization(Nphs2_proteome_filt, Nphs2_proteome_Norm)
  
  # imputation
  # plot_detect( MRpodo_proteome_se1 )
  #   plot_detect( MRpodo_proteome_se2 )
  # 
    MRpodo_proteome_se1  <- impute( MRpodo_proteome_se1, fun = "MLE")
  MRpodo_proteome_se2  <- impute( MRpodo_proteome_se2, fun = "MLE")
  
}

### read human DDA MS PXD054062
{
  
  ## read the data
  library(mzID)
  
  # Read and parse the mzIdentML file (automatically handles gzip files)
  mzidData <- mzID("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/PROTEOME/NEFR738_all_LFQ_Paper.mzid.gz")
  psmData <- flatten(mzidData)
  
  # Inspect the PSM-level data (peptide spectrum matches)
  head( psmData )
 
  
    colnames( psmData ) <- sub( "_2022.*", "", 
                                   gsub( "X:\\\\Didier\\\\Data\\\\NEFR\\\\|.raw|_2022.*|_DIA", "", 
                                         colnames(proteome_data2)) )
  
  ## make unique
  proteome_unique1 <- DEP::make_unique( proteome_data1 , names="Genes", 
                                        ids="Genes", delim = ";")

  ## 
  MRpodo_proteome_se1 <- DEP::make_se( proteome_unique1 ,
                                       columns = 2:32, 
                                       expdesign =  metadata2[metadata2$batch==5,] )
 
  # Filter for proteins that are identified in all replicates of at least one condition
  MRpodo_proteome_se1 <- filter_missval(MRpodo_proteome_se1, 
                                        thr = round(min(summary(
                                          as.factor(metadata2[metadata2$batch==5,
                                                              "condition"])))*filt.thrsh))

  # # Plot a barplot of the protein identification overlap between samples
  # print( plot_coverage(MRpodo_proteome_se1) )

  # Normalize the data
  MRpodo_proteome_se1 <- normalize_vsn(MRpodo_proteome_se1)

  # meanSdPlot( MRpodo_proteome_se1 )
  # 
  # # Visualize normalization by boxplots for all samples before and after normalization
  # plot_normalization(Nphs2_proteome_filt, Nphs2_proteome_Norm)
  
  # imputation
  # plot_detect( MRpodo_proteome_se1 )
  # 
  MRpodo_proteome_se1  <- impute( MRpodo_proteome_se1, fun = "MLE")

}

### read MaxQuant, KFO
{
  ll <- list.files( pattern = "proteinGroups", full.names = T , 
                    path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/data/")
  
  podoGlom_proteome.raw <- lapply( ll , read.csv, header = TRUE,
                                   stringsAsFactors = FALSE, sep = "\t")
  # add gene names to first dataset
  X <- (strsplit(podoGlom_proteome.raw[[1]]$Protein.IDs, split = "\\||;"))
  X1 <- unlist( lapply(X, function(x) gsub( "_MOUSE", "",paste( unique(x[grep("MOUSE",x)]), collapse = ";") )))
  podoGlom_proteome.raw[[1]]$Gene.names <- X1
  
  # read the data, make unique protein names and filter out  contaminants
  podoGlom_proteome <- lapply( podoGlom_proteome.raw , function(proteome) {
    # remove decoy matches and matches to contaminant
    proteome = proteome[!proteome$Reverse=="+",]
    proteome = proteome[!proteome$Potential.contaminant=="+",]
    
    # Make unique names using the annotation in the "Gene.names" column as primary names 
    # and the annotation in "Protein.IDs" as name for those that do not have an gene name.
    proteome_unique <- DEP::make_unique( proteome , names="Gene.names", 
                                         ids="Protein.IDs", delim = ";")
    proteome_unique$name <- toupper( proteome_unique$name )
    return(proteome_unique)
    
  })
  
  # Generate a SummarizedExperiment object using an experimental design
  ll.meta <- list.files( pattern = "experimental_design", full.names = T , 
                         path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/data/")
  
  # creating se objects
  podoGlom_proteome_se <- lapply(seq(podoGlom_proteome), function( ii )
  {
    proteome <- podoGlom_proteome[[ii]]
    print(ii)
    MRpodo_LFQ_columns <- grep("LFQ.", colnames( proteome )) # get LFQ column numbers
    experimental_design <- read.table(header=T, sep = "\t",ll.meta[[ii]])
    MRpodo_proteome_se <- make_se( proteome , MRpodo_LFQ_columns, experimental_design )
    return(MRpodo_proteome_se)
  })
  
  # filtering
  podoGlom_proteome_filt <- lapply(seq(podoGlom_proteome_se), function( ii )
  {
    print(ii)
    proteome <- podoGlom_proteome_se[[ii]]
    experimental_design <- read.table(header=T, sep = "\t",ll.meta[[ii]])
    
    # Filter for proteins that are identified in all replicates of at least one condition
    proteome <- filter_missval(proteome, thr = round(
      min(summary( as.factor(experimental_design$condition)))*filt.thrsh))
    
    # Plot a barplot of the protein identification overlap between samples
    # print( plot_coverage(proteome) )
    
    # Normalize the data
    proteome <- normalize_vsn(proteome)
    # meanSdPlot( proteome )
    
    # # Visualize normalization by boxplots for all samples before and after normalization
    # plot_normalization(Nphs2_proteome_filt, Nphs2_proteome_Norm)
    
    # imputation
    # plot_detect( proteome )
    proteome  <- impute( proteome, fun = "MLE")
    
  })
  
  names(podoGlom_proteome_filt)<- sub("_.*","",basename(ll))
  
}


#### check expression of Oxphos and Glycolisis genes ####
podoGlom_proteome_filt <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/PROTEOME/Prot.Glom.datt_.rda")

### select gene and protein names 
  {
  ggenes1 <- keggPath_pathlist[["Glycolysis / Gluconeogenesis"]][ keggPath_pathlist[["Glycolysis / Gluconeogenesis"]] %in%  allPodoGenes ]
  ggenes2 <- keggPath_pathlist[["Oxidative phosphorylation"]][ keggPath_pathlist[["Oxidative phosphorylation"]] %in%  allPodoGenes ]
  ggenes3 <- c("Pklr", "Pkm","Cycs","Atp5f1a","Tomm20",
               "ATP5A1", "TOM20","CYC","Atp5a1", "Tom20","Cyc")
  
  # convert gene to protID
  mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl" )
  gene2prot <- biomaRt::getBM( attributes=c( 'external_gene_name',"uniprotsptrembl","uniprotswissprot"),  
                               mart = mart_mouse)
  gene1.prot <-  c( toupper( ggenes1 ), 
                    gene2prot$uniprotsptrembl[ ( gene2prot$external_gene_name ) %in%  ggenes1 ] , 
                    gene2prot$uniprotswissprot[ ( gene2prot$external_gene_name ) %in%  ggenes1 ] )
  gene1.prot <- gene1.prot[gene1.prot!=""]
  gene2.prot <-  c( toupper( ggenes2 ), 
                    gene2prot$uniprotsptrembl[ ( gene2prot$external_gene_name ) %in%  ggenes2 ] , 
                    gene2prot$uniprotswissprot[ ( gene2prot$external_gene_name ) %in%  ggenes2 ] )
  gene2.prot <- gene2.prot[gene2.prot!=""]
  
  gene3.prot <-  c( toupper( ggenes3 ), 
                    gene2prot$uniprotsptrembl[ ( gene2prot$external_gene_name ) %in%  ggenes3 ] , 
                    gene2prot$uniprotswissprot[ ( gene2prot$external_gene_name ) %in%  ggenes3 ] )
  gene3.prot <- c( gene3.prot[gene3.prot!=""] ,"ATP5A1", "TOM20","CYC")   
  
  oxphos.marks <- c(gene1.prot, gene2.prot,gene3.prot)
  oxphos.marks.genes <- Reduce(union, list(ggenes1, ggenes2,ggenes3))
  
  fun_homoTO.FROMmouse(oxphos.marks.genes, TO = F)
  oxphos.marks.genes.H <- union( oxphos.marks.genes , 
                               unlist( fun_homoTO.FROMmouse(
                                 oxphos.marks.genes, TO = F)$HOMO ))
}
  
  ### transcript lvls 
  {
    listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
    
    GlyOxph.pathMean <- lapply( seq(listSCSN.1K.sampl), function(ii)
    {
      print( ii )
      
      datt <-  Seurat::ScaleData( listSCSN.1K.sampl[[ii]], features = c( ggenes1, ggenes2))
      
      datt1 <- cbind( rowMeans( datt@assays$RNA@data[ ggenes1 ,  datt$gtypeDE=="control"] ), 
                      rowMeans(  datt@assays$RNA@data[ ggenes1 ,  datt$gtypeDE!="control"] ) )
      
      datt2 <- cbind( rowMeans( datt@assays$RNA@data[ ggenes2 ,  datt$gtypeDE=="control"] ),
                      rowMeans(  datt@assays$RNA@data[ ggenes2 ,  datt$gtypeDE!="control"] ) )
      
      ll <- list( datt1 , datt2)
      names(ll) <- c("glyc","oxphos")
      return(ll)
    })
    
    Gly.pathMean <- apply(simplify2array(lapply( GlyOxph.pathMean, "[[",1) ), 
                          c(1, 2), FUN = mean, na.rm = TRUE)
    
    colnames(Gly.pathMean) <-  paste( c("ctrl","xprmnt"),"Glyc" ,  sep = "_")
    
    Oxph.pathMean <-apply(simplify2array(lapply( GlyOxph.pathMean, "[[",2) ), 
                          c(1, 2), FUN = mean, na.rm = TRUE)
    colnames(Oxph.pathMean) <-  paste(c("ctrl","xprmnt"),  "Oxphos" ,  sep = "_")
    
    
    gg1 <- ComplexHeatmap::Heatmap(( Gly.pathMean), col=viridis(10), 
                                   cluster_columns  = F, name = "Glycolisis KEGG")
    gg2 <- ComplexHeatmap::Heatmap(( Oxph.pathMean), col=viridis(10), 
                                   cluster_columns = F, name="Oxphos KEGG")
    
    pdf( width = 4, height = 24, file = "XX.pdf")
    plot_grid(grid.grabExpr(draw(gg1)), grid.grabExpr(draw(gg2)), nrow = 2, rel_heights = c(1,3))
    dev.off()
    
    ### one gene 
    Ggene <- "Dpp4"
    Ggene <- c("Clock","Bmal1","Arntl","Bmal2","Arntl2","Npas2",
               "Nr1d1" , "Nr1d2" , "Rora", "Rorb", "Rorc" )
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
    
 pdf( width = 18, height = 12, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Revision1/Circ.TFs.exprSCSN_trnscrpt.lvl.pdf")
 print(gg0)
 dev.off()
    
    
  }
  
  ### plot protein expression of markers
  {
    
    
    ### custom jitterplot
    toPlot <- Reduce( rbind, lapply( seq(podoGlom_proteome_filt) , 
                                     function( ii, # ggenes= c( gene1.prot, gene2.prot)
                                               ggenes= toupper(Ggene ))
    {
      print(ii)
      
      ### prepare the data
      
      XX <- podoGlom_proteome_filt[[ii]]@assays@data[[1]]
      # print(table(rownames(XX) %in% ggenes ))
      
      if( dim(table(rownames(XX) %in% ggenes ))>1 ){
        
        XX <- XX[  rownames(XX) %in% ggenes , , drop=F]
        # remove genes with less than 3 measurments 
        XX <- as.data.frame(  XX[rowSums( !is.na(XX) )>=3, , drop=F] )
        
        XX$name <- rownames(XX)
        
        XX.m <- reshape2::melt( XX )
        ## add metadata
        XX.m$condition <- sub( "_.*", "", XX.m$variable)
        XX.m$study <- names(podoGlom_proteome_filt)[ii]
        # print( XX.m )
      } else  XX.m <- as.data.frame( matrix(NA, nrow = 0, ncol = 5, 
                                            dimnames = list(NULL, c("name","variable",
                                                                    "value","condition", "study"))))
      
      return(XX.m)
      
    })     )
    
    # toPlot$pathN <- ifelse(toPlot$name %in%  gene1.prot , "Glycolisis",ifelse(toPlot$name %in%  gene2.prot ,"Oxphos", "other" ))
    
    toPlot$condition <- ifelse(  toPlot$condition%in% c("c", "Con","control","Ctrl") , "control",  toPlot$condition)
    toPlot$condition <- ifelse(  toPlot$condition%in% c("X5day","X9day") , "5to9_day", 
                                 toPlot$condition)
    toPlot$condition[ toPlot$condition=="d"] <- "doxo"
    toPlot$condition[ toPlot$condition=="Podocin"] <- "Nphs2mut"
    
    toPlot$condition <- relevel( as.factor( toPlot$condition ) , ref = "control")
    toPlot$conditionBin <- ifelse(    toPlot$condition%in% c("control","Podo"), "control", "experimental")
    
    # plot for multiple genes
    ggl <- ggplot2::ggplot( data = toPlot ,  
                            aes( x = value , 
                                 # y = reorder(name, value, FUN = mean, na.rm=T) ) 
                                 y =condition)
                            ) + 
      geom_jitter(shape=16, alpha=0.5, position=position_jitter(width =0.1, height = 0.2),
                  aes( color=conditionBin),
                 size = 3) + 
      facet_grid(
         rows = vars(name) , 
        # rows = vars(condition) , 
        
        cols = vars(study),
        scales = "free_y",
        space = "free") + 
      scale_colour_colorblind(name = "proteomics\ndataset") +
      theme_bw() + 
      # labs(title = sub("_.*","",basename(ll)[ii]) ) +
      labs(y = "Gene symbol", x = "LFQ itensity") +
      theme(text=element_text(size=20 ), axis.text.y = element_text(size=14 )) +
      guides(color = guide_legend(override.aes = list(size=5)))
    
    pdf( height = 4, width = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/figures/NewGlyPhos_protein.lvl.pdf")
    print(ggl)
    dev.off()
    
    
    ### one gene plot summary
    # plot for one genes
    ggl2 <- ggplot2::ggplot( data = toPlot ,  
                            aes( x =conditionBin   , 
                                 # y = reorder(name, value, FUN = mean, na.rm=T) ) 
                                 y = value )  ) + 
      facet_wrap( vars(study),
        scales = "free", nrow = 2,
        # space = "free_x"
        ) + 
      geom_jitter(shape=16, alpha=0.5, position=position_jitter(width =0.1, height = 0.2),
                  aes( color=conditionBin),
                  size = 3) + 
        scale_colour_colorblind(name = "proteomics\ndataset") +
      theme_bw() +       stat_compare_means(label = "p.format")+

      # labs(title = sub("_.*","",basename(ll)[ii]) ) +
      labs( y = "DPP4 LFQ itensity", x = "studies") +
      theme(text=element_text(size=14 ), axis.text.y = element_text(size=14 ),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      guides(color = guide_legend(override.aes = list(size=5)))
    
    pdf( height = 8, width = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/figures/DPP4.Expr_protein.lvl.pdf")
    print(ggl2)
    dev.off()
    
    ###
    toPlot.sel <- toPlot[toPlot$study=="PXD047417",]
    ggl3 <- ggplot2::ggplot( data = toPlot.sel ,  
                             aes( x =condition   , 
                                  # y = reorder(name, value, FUN = mean, na.rm=T) ) 
                                  y = value )  ) + 
      geom_jitter(shape=16, alpha=0.5, position=position_jitter(width =0.1, height = 0.2),
                  aes( color=condition),
                  size = 3) + 
      scale_colour_colorblind(name = "proteomics\ndataset") +
      theme_bw() +       
      stat_compare_means(label = "p.format", method = "wilcox.test")+
      # labs(title = sub("_.*","",basename(ll)[ii]) ) +
      labs( y = "DPP4 LFQ itensity", x = "studies") +
      theme(text=element_text(size=14 ), axis.text.y = element_text(size=14 ),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap( vars(name),
                  scales = "free", nrow = 2,
                  # space = "free_x"
      ) + 
      guides(color = guide_legend(override.aes = list(size=5)))
    
    pdf( height = 8, width = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/figures/Circ.TF.Expr_protein.lvl_PXD047417.pdf")
    print(ggl3)
    dev.off()
    
  }
  
  ### test DE, show GlyPhos genes
  {
    
    ttest <- c("X9day_vs_control", "d_vs_c", "LPS_vs_Con", "Endo_vs_Podo", "Podocin_mut_vs_control",
               "Nephrin_vs_Ctrl", "PLA2R.MN_vs_control", "PLA2R.MN_vs_control")
    
    datt <- c( podoGlom_proteome_filt , MRpodo_proteome_se , MRpodo_proteome_se1, MRpodo_proteome_se2)
    names(datt ) <- c(sub("_.*","",basename(ll)),
                      "PXD047417", "PXD054062_NEFR849","PXD054062_NEFR1028" )
    
    protDE.list  <- lapply( seq(datt), function(ii){
      datt <- datt[[ii]]
      data_diff <- DEP::test_diff( datt ,
                                   type = "manual",
                                   test = ttest[ii])
      dep <- add_rejections(data_diff, alpha = 1, lfc = 0)
      
      return( dep )
    }) 
    
   
    
    plotList <- Reduce( rbind, lapply( seq(protDE.list), function(ii)
      {
      
      dep <- protDE.list[[ii]]
      
      # XX <- DEP::get_results( dep  )
      # 
      # if( ii == 1){
      #   X <- (strsplit(XX$ID, split = "\\||;"))
      #   X1 <- unlist( lapply(X, function(x) gsub( "_MOUSE", "",paste( unique(x[grep("MOUSE",x)]), collapse = ";") )))
      #   X2 <- unlist( lapply(X, function(x) gsub( ".*CON__", "",paste( unique(x[grep("P\\d.*|Q\\d.*",x)]), collapse = ";") )))
      #   XX$ID <- ifelse( X1=="", X2 , X1)  
      # }
      # oxphos.DE <-XX[ XX$ID %in% c(gene1.prot, gene2.prot,gene3.prot), ] 
      
  
      toPlot <- plot_volcano( dep, contrast = ttest[ii] , plot = F)
  
  toPlot$OxPhos <- ifelse(toPlot$protein   %in%  toupper( oxphos.marks.genes.H ),
                          "oxphos","")
  
  toPlot$study <- names(datt )[ii]
  toPlot$label2 <- ifelse(toPlot$protein   %in%  oxphos.marks.genes.H ,
                          toPlot$protein ,"")
  return(toPlot)
} ) )
    
    plotList <- plotList[ order(plotList$OxPhos),]
    
    gg1 <- ggplot(plotList, aes(x=log2_fold_change, y=plotList$'p_value_-log10',
                                label=label2, 
                                color= OxPhos )) + 
      geom_vline(xintercept=0, color="red")+
      geom_point( alpha=0.5)+ theme_bw()+ 
      # scale_color_colorblind()+
      ggrepel::geom_text_repel()+ 
      facet_wrap(vars(study), scales = "free")+
      theme( text = element_text( size=20))+ ggtitle("protein DE in Nphs2mut.")
    
    
    toPlot <- plotList[ plotList$OxPhos != "" ,]
    toPlot$newMark <- ifelse(toPlot$protein   %in% (toupper( ggenes3))  ,
                               "newMark" ,"")
    toPlot <- toPlot[ order(toPlot$newMark),]
    toPlot$label2 <- ifelse(toPlot$protein   %in%  (toupper( ggenes3)) ,
                            toPlot$protein ,"")
    
    gg2 <- ggplot(toPlot, aes(x=log2_fold_change, y=toPlot$'p_value_-log10',
                                label=label2, 
                                color= newMark )) + 
      geom_vline(xintercept=0, color="red")+
      geom_point( alpha=0.5)+ 
      theme_bw()+ 
      scale_color_colorblind()+
      ggrepel::geom_text_repel()+ 
      facet_wrap(vars(study), scales = "free")+
      theme( text = element_text( size=20))+ ggtitle("protein DE in Nphs2mut.")
    
    pdf( height = 24, width = 24, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/figures/GlyPhos.newMarks_protein.DE_M.H.pdf")
    print(gg2)
    dev.off()
    
    
    ###
    protDE.list_tabs <-  lapply( seq(protDE.list), function(ii)
       DEP::get_results( protDE.list[[ii]]  )  ) 
    names(protDE.list_tabs)<- names(datt)

  }

  
 
  


#### apply PDS to proteomics ####
MtoH <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/DSpodo_musTOhomo_oneTOone.tsv")

DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
DS_all$gene_symbol <- MtoH$Mouse_Symbol
DS_all.PROT <- DS_all
DS_all.PROT$gene_symbol <- toupper(DS_all.PROT$gene_symbol)
DS_all.PROT$gene_symbol.HS <- MtoH$Human_Symbol

###
# bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
datt <- c( podoGlom_proteome_filt , MRpodo_proteome_se , MRpodo_proteome_se1, MRpodo_proteome_se2)
names(datt ) <- c(sub("_.*","",basename(ll)),
                  "PXD047417", "PXD054062_NEFR849","PXD054062_NEFR1028" )

PDS_prot.M <- Reduce( rbind, lapply( c(1,2,4:6), function(ii){
  print(ii)
  
  XX <- datt[[ii]]@assays@data[[1]]
  # print( table(  DS_all.PROT$gene_symbol[1:42]%in%  rownames(XX)))
  # print( table(  DS_all.PROT$gene_symbol[1:60]%in%  rownames(XX)))
  
  dat <-tryCatch( { 
    data.frame( PDS=DS_calc.func(exprMatrices = XX , DSignature = DS_all.PROT,
               ntop = 60, ceilThrsh = 0.1) ,
               sample= colnames(XX), 
               study = names(datt )[ii])
    }, error = function(e) NA )
}) )

PDS_prot.H <- Reduce( rbind, lapply( c(7:8), function(ii){
  print(ii)
  
 
  
  XX <- datt[[ii]]@assays@data[[1]]
  # print( table(  DS_all.PROT$gene_symbol.HS[1:42]%in%  rownames(XX)))
  # print( table(  DS_all.PROT$gene_symbol.HS[1:60]%in%  rownames(XX)))
  # print( table(  DS_all.PROT$gene_symbol.HS[1:100]%in%  rownames(XX)))  
  
  dat <-tryCatch( { 
    data.frame( PDS= DS_calc.func(exprMatrices = XX , 
                                 geneIDname = "gene_symbol.HS",
                                 DSignature = DS_all.PROT,
                                 ntop = 60, ceilThrsh = 0.1) ,
                sample= colnames(XX), 
                study = names(datt )[ii])
  }, error = function(e) NA )
}) )

## combine mouse and human
PDS_prot <- rbind(PDS_prot.M , PDS_prot.H)
PDS_prot$condition <- sub( "_.*", "", PDS_prot$sample)
# PDS_prot$condition <- ifelse( grepl( "_.*\\.", PDS_prot$sample), 
#                               sub("_.*\\.", "_", PDS_prot$sample), 
#                               PDS_prot$condition )

ggl <- lapply( seq( unique(PDS_prot$study)), function(ii){
  toPlot <- PDS_prot[ PDS_prot$study == unique(PDS_prot$study)[ii] , ]
  
  if(ii==3){ toPlot$condition <- relevel(as.factor(toPlot$condition ), ref = "Podo") }
  ggplot(toPlot, aes(y=PDS, x=condition,
                       color= condition ), colo="darkgrey") + 
    geom_boxplot()+
    geom_point( alpha=0.5, size = 3)+ 
    stat_compare_means(method = "t.test" )+
    theme_bw()+ 
    scale_color_colorblind()+
    # facet_grid(. ~ study, scales = "free_x", space = "free") +
    theme( text = element_text( size=20), legend.position = "bottom")+
    ggtitle(unique(PDS_prot$study)[ii] )
})


 cowplot::plot_grid(plotlist = ggl)
# dev.off()
# dev.off()

pdf( height = 15, width = 12, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/figures/PDS.protein.test.Updt.pdf")
print(cowplot::plot_grid(plotlist = ggl[c(3,3,1,2,4,5,6,7)], ncol = 2) )
dev.off()
