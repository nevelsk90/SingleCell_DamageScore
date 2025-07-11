###=====# Single-Cell Resolution of Cellular Damage Illuminates Disease Progression #=======###
###=====# supplementary figure 2 PDS code #=======###
mallinfo::mallinfo()
mallinfo::malloc.trim()

options(connectionObserver = NULL)
library(DEP)
library(enrichplot)


# load necessary code
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")  

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")


### top 42 markers
DS_all.42 <- DS_all[ 1:42, ]
# all genes UP
upp <- DS_all.42$gene_symbol[DS_all.42$direction_foldchange==1 ]
upp  <-  c( toupper( upp ), 
            tx2gene$uniprotsptrembl[ ( tx2gene$external_gene_name ) %in%  upp ] , 
            tx2gene$uniprotswissprot[ ( tx2gene$external_gene_name ) %in%  upp ] )
# all genes DOWN
downn <- DS_all.42$gene_symbol[DS_all.42$direction_foldchange==-1 ]
downn <- c( toupper( downn ) ,
            tx2gene$uniprotsptrembl[ ( tx2gene$external_gene_name ) %in%  downn ] , 
            tx2gene$uniprotswissprot[ ( tx2gene$external_gene_name ) %in%  downn ] )


#### Suppl Fig. 2A protein lvl expression of PDS markers ####
  

   # read PXD016238 data
    proteome <- read.csv( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/data/PXD016238_proteinGroups.txt", header = TRUE,
                                   stringsAsFactors = FALSE, sep = "\t")
 
  
  # read the data, make unique protein names and filter out  contaminants
  # remove decoy matches and matches to contaminant
    proteome = proteome[!proteome$Reverse=="+",]
    proteome = proteome[!proteome$Potential.contaminant=="+",]
    
    # Make unique names using the annotation in the "Gene.names" column as primary names 
    # and the annotation in "Protein.IDs" as name for those that do not have an gene name.
    proteome_unique <- DEP::make_unique( proteome , names="Gene.names", 
                                         ids="Protein.IDs", delim = ";")
    proteome_unique$name <- toupper( proteome_unique$name )

 
  # Generate a SummarizedExperiment object using an experimental design

    MRpodo_LFQ_columns <- grep("LFQ.", colnames( proteome_unique )) # get LFQ column numbers
    experimental_design <- read.table(header=T, sep = "\t","/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/data/PXD016238_experimental_design.txt")
    MRpodo_proteome_se <- make_se( proteome_unique , MRpodo_LFQ_columns, experimental_design )

  # Filter for proteins that are identified in all replicates of at least one condition
    proteome <- filter_missval(MRpodo_proteome_se, thr = round(min(summary( as.factor(experimental_design$condition)))*0.49))
    
    # Plot a barplot of the protein identification overlap between samples
    print( plot_coverage(proteome) )
    
    # Normalize the data
    proteome <- normalize_vsn(proteome)
    meanSdPlot( proteome )
    
    # # Visualize normalization by boxplots for all samples before and after normalization
    # plot_normalization(Nphs2_proteome_filt, Nphs2_proteome_Norm)
    
    # imputation
    plot_detect( proteome )
    proteome  <- impute( proteome, fun = "MLE")
    

  # DE analysis
    podoGlom_proteome_DE <- test_diff( proteome, type = "all")
    podoGlom_proteome_DE <- add_rejections(podoGlom_proteome_DE, alpha = 0.1, lfc=0)
  
  
 
  ### boxplots of proteins of interest
  lapply( seq(podoGlom_proteome_DE) , function(ii){
    print(ii)
    
    # check if any proteins (up) of interest are detected and if yes - plot
      plot_single( podoGlom_proteome_DE , type= "contrast",
                          proteins = c("WT1","NPHS2") )
    

      pp2 <- plot_single( podoGlom_proteome_DE,  type= "contrast",
                          proteins =  downn )
 
      
    # plot 
    if( exists("pp1") && exists("pp2")) {
      gg <- cowplot::plot_grid( pp1, pp2 , labels = c("UP","DOWN") )
    } else if( exists("pp1")){
      gg <- cowplot::plot_grid( pp1, labels = "UP" )
    } else if( exists("pp2")){
      gg <- cowplot::plot_grid( pp2, labels = "DOWN" )
    } else gg <- NULL
    
    
    # save plots
    pdf(width = 10, height = 10, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                               sub("_.*","",basename(ll)[ii]),"_PDS.50_contrasts.pdf"))
    print( gg )
    dev.off()
    png(width = 800, height = 800, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                 sub("_.*","",basename(ll)[ii]),"_PDS.50_contrasts.png"))
    print( gg )
    dev.off()
  } )
  
  ### custom jitterplot
  lapply( seq(podoGlom_proteome_se) , function( ii )
  {
    
    ### prepare the data
    XX <- podoGlom_proteome_se[[ii]]@assays@data[[1]]
    XX <- XX[  rownames(XX) %in% c(upp, downn), ]
    # remove genes with less than 3 measurments 
    XX <- as.data.frame(  XX[rowSums( !is.na(XX) )>=3, ] )
    
    XX$name <- rownames(XX)
    
    XX.m <- reshape2::melt( XX )
    ## add metadata
    XX.m$condition <- sub( "_.*", "", XX.m$variable)
    XX.m$direction <- ifelse( XX.m$name %in%  downn ,
                              "DOWN", "UP")
    
    # plot 
    gg <- ggplot2::ggplot( data = XX.m ,  
                           aes( x = value , 
                                y = reorder(name, value, FUN = mean, na.rm=T) ) ) + 
      geom_jitter(shape=16, position=position_jitter(width =0.1, height = 0.1),
                  size = 5, aes( color=condition)) + 
      facet_grid(rows = vars(direction) ,scales = "free_y", space = "free_y") + 
      scale_colour_colorblind(name = "proteomics\ndataset") +
      theme_bw() + labs(title = sub("_.*","",basename(ll)[ii]) ) +
      labs(y = "Gene symbol", x = "LFQ itensity") +
      theme(text=element_text(size=20 )) +
      guides(color = guide_legend(override.aes = list(size=5)))
    
    # save plot
    pdf(width = 10, height =  nrow( XX)*0.4, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                           sub("_.*","",basename(ll)[ii]),"_PDS.50_jitterplotNOTscaled.pdf", sep = ""))
    print(gg)
    dev.off()
    png(width = 800, height =  nrow( XX)*40, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                           sub("_.*","",basename(ll)[ii]),"_PDS.50_jitterplotNOTscaled.png", sep = ""))
    print(gg)
    dev.off()
  })
  
  ### average expression in podocytes
  controlIDs <- c("control", "c_", "Con_","Podo","control")
  podoGlom_exprLVL <- Reduce( function( x , y ) merge( x, y, by="Row.names", all=T) , 
                              lapply( seq(podoGlom_proteome_se), function(ii){
                                XX <- podoGlom_proteome_se[[ii]]@assays@data[[1]]
                                XX <- scale(XX ,center = F )
                                XX <- XX[  rownames(XX) %in% c(upp, downn) , ]
                                XX <- XX[ rowSums( is.na(XX)) < ncol(XX)-1,]
                                print( dim(XX))
                                XX <- data.frame( value = rowMedians( XX, na.rm = T) ,
                                                  row.names = rownames(XX))
                                XX$Row.names <- rownames(XX) 
                                return(XX)
                              }) )
  rownames(podoGlom_exprLVL) <- podoGlom_exprLVL$Row.names
  podoGlom_exprLVL <- podoGlom_exprLVL[,-1]
  
  X <- rowMeans( podoGlom_exprLVL, na.rm = T)
  X <- X[order(X)]
  par(mar=c(4,6,2,2))
  barplot( X[order(X)], las=2, horiz = T, xpd = FALSE ,xlim = c(0.8, max(X)),
           main="average protein expression (scaled)\nof PDS markers")
  
  
  

#### Suppl Fig. 2F functional annotation of PDS genes ####
  
  gene.bckgrnd_eID <- tx2gene$entrezgene_id[ match( allPodoGenes, tx2gene$external_gene_name)]
  gene.bckgrnd_eID <- as.character( gene.bckgrnd_eID[!is.na(gene.bckgrnd_eID)] )
  tx2gene <- getBM(attributes=c( "entrezgene_id", "external_gene_name"),  mart = mart_mouse)
  DS.42.eID <- as.character( tx2gene$entrezgene_id[ match( DS_all$gene_symbol[1:42],
                                                           tx2gene$external_gene_name)] )
  
  # test with cluster profiler
  DS.42_cProfiler <- cProfiler.GKR( ggenes.eID=DS.42.eID , 
                                    gene.bckgrnd_eID = gene.bckgrnd_eID )
  
  
  ### cluster terms
  gg<- enrichplot::emapplot_cluster( enrichplot::pairwise_termsim(
    DS.42_cProfiler$GO.enrich), color = "pvalue", cex_label_group=1.2,
    showCategory = 15)
  # p3 <-  enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
  #   datt$KEGG.enrich), color = "pvalue")
  # p2  <- enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
  #   DS.42_cProfiler$REACT.enrich), color = "pvalue", cex_label_group=1.2)
  # cowplot::plot_grid( plotlist =  list(p1,p2), nrow = 1, labels = c("GO","REACT"))
  pdf(height = 4, width = 6, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/DS.42_cProfiler_GOclustered")
  print(gg)
  dev.off()
  

