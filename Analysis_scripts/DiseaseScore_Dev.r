# ############################## #
### disease score development ###
# ############################## #

library(ggplot2 )
library(ggthemes)
library(plyr	)
library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
tx2gene <- getBM( attributes=c('ensembl_gene_id', 'external_gene_name',"refseq_mrna"),  mart = ensembl)

#### analyse public transcriptomic data from podocyte damage models #### 
  {
  
    # ###  GSE969
    # library("GEOquery")
    # gse=getGEO(filename="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA/GSE969_series_matrix.txt")
    # GSE969_expr <- gse@assayData$exprs[,-c(1,6,12)]
    # GSE969_expr <- data.frame ( GB_ACC = gse@featureData@data$GB_ACC , lfc = rowMeans( GSE969_expr ) )
    # GSE969_expr <- GSE969_expr[GSE969_expr$GB_ACC!="",]
    # GSE969_expr$Gene.symbol <- read.table( sep = "\t",header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA/DE/bioDBnet_db2db_GSE969.txt")[,2]

      
  #### MA studies analysed with GEO2r
  {
    ### read the data
    ll <- list.files( pattern = ".tsv",path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA/DE", full.names = T)
    FSGS_MA_DE <- lapply( ll, read.delim, header=T , sep="\t", fill=T)
    # add gene names if not available
    # FSGS_MA_DE[[2]]$Gene.symbol <- tx2gene$external_gene_name[ match( FSGS_MA_DE[[2]]$GB_ACC, tx2gene$refseq_mrna)]
    FSGS_MA_DE[[3]] <-  FSGS_MA_DE[[3]][ FSGS_MA_DE[[3]]$GB_ACC!="",]
    FSGS_MA_DE[[3]]$Gene.symbol <- read.table( sep = "\t",header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA/DE/bioDBnet_db2db_ob1.txt")[,2]

    ### unify col names 
    FSGS_MA_DE <- lapply(FSGS_MA_DE, function(x) {
      colnames(x) <- sub("GENE_SYMBOL","Gene.symbol",colnames(x)) 
      return(x)
    })
    
    ### split gene names for probes with multigenes
    # define a function to crreate duplicate rows with splitted genes
    # https://stackoverflow.com/questions/38499032/repeat-the-rows-in-a-data-frame-based-on-values-in-a-specific-column
    duplicate_rows <- function( P.Value , logFC , gNames ) {
      expanded_samples <- unlist( strsplit(gNames , split = "///") )
      repeated_rows <- data.frame( "P.Value"= P.Value, "logFC"=logFC ,"Gene.symbol" = expanded_samples )
      repeated_rows
    }
    
    # apply the function
    FSGS_MA_DE <- lapply(FSGS_MA_DE, function(x) {
      x <- x[x$Gene.symbol!="",]
      expanded_rows <- Map(f = duplicate_rows, x$P.Value, x$logFC , x$Gene.symbol)
      })
    
    # and assemble dataframes back 
    FSGS_MA_DEdf <-lapply(FSGS_MA_DE, function(x) do.call( rbind , x) )
    
    
    ### aggregate over genes
    # using geometric or arythmetic mean for combining p-values is OK for dependent p-s 
    # https://arxiv.org/pdf/1212.4966.pdf
    FSGS_MA_DEdf <- lapply( FSGS_MA_DEdf , function(x){
      # print( dim(x) )
      aggregate( . ~ Gene.symbol , data=x[,c("Gene.symbol", "logFC","P.Value")] ,  mean )
    } )
    ### make a df
    # union of all gNames
    gg <- Reduce( union, lapply(FSGS_MA_DEdf, function(x) x[,grepl("Gene.symbol",colnames(x))]))
    ### combine datasets
    # combine LFCs
    FSGS_MA_lfcDE <- Reduce( cbind, lapply(FSGS_MA_DEdf, function(x) x[ match(gg, x[,grepl("Gene.symbol",colnames(x))] ),"logFC"]))
    # combine p-values
    FSGS_MA_pvalDE <- Reduce( cbind, lapply(FSGS_MA_DEdf, function(x)  x[ match(gg, x[,grepl("Gene.symbol",colnames(x))] ),"P.Value"]))
    # combine ranks
    FSGS_MA_ranksDE <- Reduce( cbind, lapply(FSGS_MA_DEdf, function(x) rank ( na.last="keep", x[ match(gg, x[,grepl("Gene.symbol",colnames(x))] ),"P.Value"])) )
    # assign col and row names
    rownames(FSGS_MA_pvalDE) <- rownames(FSGS_MA_lfcDE) <- rownames( FSGS_MA_ranksDE) <-  gg
    colnames(FSGS_MA_pvalDE) <-  colnames(FSGS_MA_lfcDE) <- colnames(FSGS_MA_ranksDE) <- sub(".tsv","_MA", basename( ll ))
    
    saveRDS( list( FSGS_MA_ranksDE , FSGS_MA_lfcDE), file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA_lfcRank15.12.20.rda")
  
    ### make correlation plot
    XXX <- FSGS_MA_lfcDE
    # make all entries with p-value < 0.1 (==noise)  equal 0
    XXX <- XXX * !( FSGS_MA_pvalDE > 0.1 | is.na(FSGS_MA_pvalDE) ) 
    ## add antiGBM kidney data
    # XXX <- cbind( XXX , GSE969_antiGBM_kidney=GSE969_expr$lfc[ match( gg , GSE969_expr$Gene.symbol )])
    XXX[XXX ==0] <- NA # convert 0 entries to NAs so that they don't influence correlations
    XXX <- as.matrix( XXX )
    XXX[XXX <= -Inf] <- min(XXX[is.finite(XXX)], na.rm = T)
    XXX[XXX >= Inf] <- max(XXX[is.finite(XXX)], na.rm = T)
    
    #
    toPlot <- cor( XXX , use = "pairwise.complete.obs", method = "spearman" )
    res1 <- matrix(p.adjust (cor.mtest( XXX , method = "spearman" )$p),ncol=15)
    
    library("corrplot")
    corrplot::corrplot( toPlot , hclust.method = "ward.D",p.mat = res1,
                        sig.level = 0.2 ,pch.cex = 1,
                        order="hclust", tl.cex=0.8, cl.cex=0.8, tl.col="black"  )
    
    }
  
  #### read bulk RNAseq studues
  {
    
    # GSE126217 , cosmc
    Cosmc <- readxl::read_xlsx( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE126217_Cosmc_Glom_RNAseq_processed_readcounts.xlsx")
    
    ### DE analysis using characteristic directions   library(GeoDE)
    {
      ### GSE154955 , adriomyc
      XX <- read.table( header=T , "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE154955_gene.FPKM.table.control.D9.D14.txt")
      Adriom <- aggregate( . ~ Gene.symbol , data=XX ,  mean ) # collapse duplicate genes
      Adriom <- Adriom[ rowSums(Adriom[,2:ncol(Adriom)]) > 0.1 , ] # remove noise
      # perform DE analysis
      Adriom_DE <- chdirAnalysis ( Adriom , sampleclass = as.factor(c(1,1,1,rep(2,6))) , nnull=20 , CalculateSig=F )
      Adriom_DE$lfc <- log2( rowMeans(Adriom[,5:ncol(Adriom)])/rowMeans(Adriom[,2:4]))
      
      ### GSE123179 , Cd2apFyn 
      XX <- readxl::read_xlsx( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE123179_MAFB_NRPKM_2_15_17.xlsx")
      Cd2apFyn <- as.data.frame( XX ) [,-c(1,9)]
      Cd2apFyn <- aggregate( . ~ Gene.symbol , data=Cd2apFyn ,  mean ) # collapse duplicate genes
      Cd2apFyn <- Cd2apFyn[ rowSums(Cd2apFyn[,2:ncol(Cd2apFyn)]) > 0.1 , ]
      Cd2apFyn_DE <- chdirAnalysis ( Cd2apFyn , sampleclass = as.factor(c(2,1,1,2,2,1)) , nnull=20 , CalculateSig=F )
      Cd2apFyn_DE$lfc <- log2(rowMeans(Cd2apFyn[,c(2,5,6)])/rowMeans(Cd2apFyn[,c(3,4,7)]))
      
      ### GSE117987 , HIV(Tg26) adriomycin
      XX <- read.table( header=T ,sep = ",", "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE117987_repRpkmMatrix_featureCounts.csv")
      HIVtg26.Adriom <- XX[ , c(1,2,4,6,11,12,13)]
      HIVtg26.Adriom <- aggregate( . ~ Gene.symbol , data=HIVtg26.Adriom ,  mean ) # collapse duplicate genes
      HIVtg26.Adriom <- HIVtg26.Adriom[ rowSums(HIVtg26.Adriom[,2:ncol(HIVtg26.Adriom)]) > 0.1 , ]
      HIVtg26.Adriom_DE <- chdirAnalysis ( HIVtg26.Adriom , sampleclass = as.factor(c(1,2,2,2,1,1)) , nnull=20 , CalculateSig=F )
      HIVtg26.Adriom_DE$lfc <- log2(rowMeans(HIVtg26.Adriom[,c(3,4,5)])/rowMeans(HIVtg26.Adriom[,c(2,6,7)]))
      
      
     
      ### GSE79291 streptoz 
      # !! not all podocyte samples are available on GEO/SRA, use gloms for now
      XX <- readr::read_delim(delim = "\t","/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE79291_Final_Normalized_Data.txt")
      XX <- cbind( Gene.symbol = XX$Symbol, XX[,c(29:36)])
      streptoz <- aggregate( . ~ Gene.symbol , data=XX ,  mean ) # collapse duplicate genes
      streptoz <- streptoz[ rowMeans(streptoz[,2:ncol(streptoz)]) > 0.1 , ]
      streptoz[,-1] <- streptoz[,-1]^2
      streptoz_DE <- chdirAnalysis ( streptoz , sampleclass = as.factor(c(1,1,1,1,2,2,2,2)) , nnull=20 , CalculateSig=F )
      streptoz_DE$lfc <- log2(rowMeans(streptoz[,2:5], na.rm = T)/rowMeans(streptoz[,6:ncol(streptoz)], na.rm = T))
      
      ### GSE134327 streptozotocin
      ll <- list.files( full.names = T , path= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE134327_RAW/")
      stz2 <- Reduce( cbind.data.frame , lapply(ll, function(x) read.delim(gzfile(x), header = T )[,c(1,20)] ) )
      stz2 <- cbind( Gene.symbol = stz2[,1], stz2[,which(colnames(stz2)!="Symbol")])
      colnames(stz2)[-1] <- sub( "_.*", "", basename(ll) )
      stz2 <- aggregate( . ~ Gene.symbol , data=stz2 ,  mean ) # collapse duplicate genes
      stz2 <- stz2[, colnames(stz2) %in% c("Gene.symbol" ,"GSM3942312" , "GSM3942313" , "GSM3942314" , "GSM3942315" ,
        "GSM3942316", "GSM3942320" ,"GSM3942321" ) ]
      stz2[,-1] <- stz2[,-1]^2
      stz2_DE <- chdirAnalysis ( stz2 , sampleclass = as.factor(c(1,1,2,2,2,1,1)) , nnull=20 , CalculateSig=F )
      stz2_DE$lfc <- log2(rowMeans(stz2[,4:6], na.rm = T)/rowMeans(stz2[,c(2,3,7,8)], na.rm = T) )
      
      ### GSE77717  db/db
      ll <- list.files( full.names = T , path= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE77717_RAW/")
      dbdb <- Reduce( cbind.data.frame , lapply(ll, function(x) read.delim(gzfile(x), header = T )[,c(2,11)] ) )
      dbdb <- cbind( Gene.symbol = dbdb[,1], dbdb[,which(colnames(dbdb)!="gene_Name")])
      dbdb <- aggregate( . ~ Gene.symbol , data=dbdb ,  mean ) # collapse duplicate genes
      dbdb <- dbdb[ rowMeans(dbdb[,2:ncol(dbdb)]) > 0.1 , ]
      dbdb_DE <- chdirAnalysis ( dbdb , sampleclass = as.factor(c(2,2,2,1,1,1)) , nnull=20 , CalculateSig=F )
      dbdb_DE$lfc <- log2(rowMeans(dbdb[,2:4], na.rm = T)/rowMeans(dbdb[,5:7], na.rm = T) )
      
      
      ### GSE123853 db/db
      ll  <- list.files( full.names = T , path= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE123853_RAW/")
      dbdb2 <- Reduce( cbind.data.frame , lapply(ll[c(48:53,60:65)], function(x) read.delim(gzfile(x), header = F ) ) )
      dbdb2 <- cbind( Gene.symbol = dbdb2[,1], dbdb2[,which(colnames(dbdb2)!="V1")])
      dbdb2 <- aggregate( . ~ Gene.symbol , data=dbdb2 ,  mean ) # collapse duplicate genes
      dbdb2 <- dbdb2[ rowMeans(dbdb2[,2:ncol(dbdb2)]) > 0.1 , ]
      dbdb2_DE <- chdirAnalysis ( dbdb2 , sampleclass = as.factor(c(rep(1,6), rep(2,6))) , nnull=20 , CalculateSig=F )
      dbdb2_DE$lfc <- log2(rowMeans(dbdb2[,2:4], na.rm = T)/rowMeans(dbdb2[,5:7], na.rm = T) )
      
      
      ### GSE110092  shroom3
      ll <- list.files( full.names = T , path= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE110092_RAW/")
      shroom3 <- Reduce( cbind.data.frame , lapply(ll , function(x) read.delim(gzfile(x), header = T ) ) )
      shroom3 <- cbind( Gene.symbol = shroom3[,1], shroom3[,which(colnames(shroom3)!="Symbol")])
      shroom3 <- aggregate( . ~ Gene.symbol , data=shroom3 ,  mean ) # collapse duplicate genes
      shroom3[,-1] <- shroom3[,-1]^2
      shroom3_DE <- chdirAnalysis ( shroom3 , sampleclass = as.factor(c(2,2,2,1,1,1,1)) , nnull=20 , CalculateSig=F )
      shroom3_DE$lfc <- log2(rowMeans(shroom3[,2:4], na.rm = T)/rowMeans(shroom3[,5:8], na.rm = T) )
      
      bulkCHRDIR <- list( Adriom_DE , Cd2apFyn_DE , HIVtg26.Adriom_DE, streptoz_DE , stz2_DE, dbdb_DE, dbdb2_DE , shroom3_DE)
      
      lfc_bulkCHRDIR <- lapply( bulkCHRDIR , function(x) {
        print(length(x$lfc ))
        names( x$lfc ) <- rownames( x$chdirprops$chdir[[1]])
        return( x$lfc) })
      
      rank_bulkCHRDIR <- lapply( bulkCHRDIR , function(x) {
        xx <- rank(-abs(x$chdirprops$chdir[[1]])[,1] , na.last="keep") } )
      
      gg <- Reduce(union, lapply(rank_bulkCHRDIR, names))
      
      lfc_bulk.chrdirDF <- Reduce( cbind , lapply(lfc_bulkCHRDIR, function(x) x[match(gg, names(x))]))
      rank_bulk.chrdirDF <- Reduce( cbind , lapply(rank_bulkCHRDIR, function(x) x[match(gg, names(x))]))
      
      rownames(lfc_bulk.chrdirDF) <- rownames(rank_bulk.chrdirDF) <-  gg
      colnames(lfc_bulk.chrdirDF) <- colnames(rank_bulk.chrdirDF) <-  c( "GSE154955_adriomyc" , "GSE123179_Cd2apFyn" , "GSE117987_HIV(Tg26)adriomycin", "GSE79291_streptoz", "GSE134327_streptoz", "GSE77717_dbdb" ,"GSE123853_dbdb", "GSE110092_shroom3")
      
      }

    ### analysis of raw reads with DEseq2
    {
      library(DESeq2)
      ### GSE138774 bpix 
      ll <- list.files( full.names = T , path= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE138774_RAW")
      bpix <- Reduce( cbind.data.frame , lapply(ll, function(x) read.table(gzfile(x))) )
      bpix <- bpix[- c((nrow(bpix)-4):nrow(bpix)) , ]
      rownames(bpix) <- bpix$V1
      colnames(bpix) <- sub( "_counts_file.txt.gz" ,"", basename(ll))
      bpix <- bpix[,which(colnames(bpix) !="V1")]
      bpix <- bpix[ rowMeans(bpix) > 1 , ]
      coldata <- cbind.data.frame(samples=colnames(bpix) , 
                                  gt=as.factor( c("KO","KO","KO","CTRL","CTRL","CTRL")) )
      # use DEseq2
      dds <- DESeqDataSetFromMatrix(countData = bpix,
                                    colData = coldata,
                                    design= ~ gt)
      dds <- DESeq(dds)
      bpix_DE <- results(dds)
      bpix_DE$Gene.symbol <-  tx2gene$external_gene_name [ match( rownames(bpix_DE), tx2gene$ensembl_gene_id)]
      # aggregate duplicated genes
      bpix_DE <- aggregate( . ~ Gene.symbol , data=bpix_DE[,c("Gene.symbol", "log2FoldChange","pvalue")] ,  mean )
      
      ### GSE119049 LCDD
      XX <- read.table( header=T ,sep = "\t", "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE119049_table_featurecounts_GLO.txt")
      LCDD <- cbind( row.names= make.unique( XX$Geneid ), XX[, 7:ncol(XX)])
      LCDD <- LCDD[ rowMeans(LCDD) >  1, ]
      coldata <- cbind.data.frame(samples=colnames(LCDD) , gt=as.factor(c(2,1,1,1,2,1,2,1,1,2,1,1,1)))
      # use DEseq2
      dds <- DESeqDataSetFromMatrix(countData = LCDD,
                                    colData = coldata,
                                    design= ~ gt)
      dds <- DESeq(dds)
      LCDD_DE <- results(dds)

      LCDD_DE <- cbind.data.frame( Gene.symbol = rownames (LCDD_DE), 
                                        log2FoldChange= LCDD_DE$log2FoldChange , 
                                        pvalue=LCDD_DE$pvalue)
      
      ### GSE136138 norm aging
      # annotation
      # gse= GEOquery::getGEO(filename="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE136138_series_matrix.txt.gz")
      sInfo <- read.table( header = T, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE136138_design.csv")
      sInfo$ageCode <- ifelse( sInfo$Sample_age=="Young" , 1,2)
      # DE
      XX <- read.table( header=T , row.names=1, sep = ",", "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk/GSE136138_Counts.csv")
      rownames(XX) <- make.unique(XX$Gene)
      agingBulk <- XX[,-1]
      agingBulk <- agingBulk[ rowMeans(agingBulk) > 1 , ]
      coldata <- cbind.data.frame(samples= colnames(agingBulk), gt= as.factor( sInfo$ageCode[ match( colnames(agingBulk), gsub("-", "\\.",sInfo$Sample_title))] ) )
      
      dds <- DESeqDataSetFromMatrix(countData = agingBulk,
                                    colData = coldata,
                                    design= ~ gt)
      dds <- DESeq(dds)
      agingBulk_DE <- results(dds)
      agingBulk_DE <- cbind.data.frame( Gene.symbol = rownames (agingBulk_DE), 
                                        log2FoldChange= agingBulk_DE$log2FoldChange , 
                                        pvalue=agingBulk_DE$pvalue)
      
      gg <- Reduce(union, lapply( list(agingBulk_DE), rownames))
      
      
    ### combine experiments
      gg <- Reduce(union, lapply( list(bpix_DE , LCDD_DE, agingBulk_DE), function(x) x$Gene.symbol))
      
      lfc_bulk.deseqDF <- Reduce( cbind , lapply(list(bpix_DE , LCDD_DE, agingBulk_DE), function(x) x$log2FoldChange[match(gg, (x$Gene.symbol))]))
      rank_bulk.deseqDF <- Reduce( cbind , lapply(list(bpix_DE , LCDD_DE, agingBulk_DE), function(x) x$pvalue[match(gg, (x$Gene.symbol))]))
      rank_bulk.deseqDF <- apply(rank_bulk.deseqDF, 2, rank)
      rownames(lfc_bulk.deseqDF) <- rownames(rank_bulk.deseqDF) <-  gg
      colnames(lfc_bulk.deseqDF) <- colnames(rank_bulk.deseqDF) <- c( "GSE138774_bpix", "GSE119049_LCDD", "GSE136138_aging" )
    
      }
    
    # # old chrdir results
    # bulk_DE <- readRDS(file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk_DE.rda")
    
    ### combine ChrDir and DEseq2
    {
      gg <- union( Cosmc$ID , union(rownames(lfc_bulk.chrdirDF), rownames(lfc_bulk.deseqDF)))
      lfc_bulk <- cbind.data.frame( lfc_bulk.chrdirDF[match( gg, rownames(lfc_bulk.chrdirDF)),],
                                    lfc_bulk.deseqDF[match( gg, rownames(lfc_bulk.deseqDF)),],
                                    Cosmc$log2FoldChange[match( gg,Cosmc$ID )])
      rank_bulk <- cbind.data.frame( rank_bulk.chrdirDF[match( gg, rownames(rank_bulk.chrdirDF)),],
                                     rank_bulk.deseqDF[match( gg, rownames(rank_bulk.deseqDF)),],
                                     rank(Cosmc$pvalue[match( gg,Cosmc$ID )]) )
      rownames(lfc_bulk) <- rownames(rank_bulk) <-  gg
      colnames(lfc_bulk) <- colnames(rank_bulk) <- paste(c( colnames(rank_bulk.chrdirDF),
                                                      colnames(rank_bulk.deseqDF),
                                                      "GSE126217_cosmc"), "bulk",sep = "_")
      
      saveRDS( list(rank_bulk, lfc_bulk) , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk_DE15.12.20.rda")
      
    }
   
    # library(tidyr)
    # library(dplyr)
    # datt <- rank_bulk
    # rank_bulk <- lapply( datt , function(x) {
    #   xx <- data.frame( Gene.symbol= names(x), lfc=x )
    #   xx <- separate_rows(xx, Gene.symbol, sep = ",")
    #   
    #   yy <- as.numeric(xx$lfc) 
    #   names(yy) <- (xx$Gene.symbol)
    #   return( yy )
    # })
    
    
   }
  
  #### eNOS/strept
  {
    library(Seurat)
    expr1 <- read.table("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/GSE127235/GSE127235_exp1_counts.txt", header = T, row.names = 1)
    expr2 <- read.table("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/GSE127235/GSE127235_exp2_counts.txt", header = T, row.names = 1)
    
    ddat <- cbind( expr1 , expr2 ) 
    rownames(ddat) <- sub("\\..*","",rownames(ddat) )
    
    cellqual <- read.delim(row.names = 1, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/sc/GSE127235_cellqual.csv")
    
    ddat <- ddat[,which(cellqual==1)]
    
    # ddat$Gene.symbol <- tx2gene$external_gene_name[ match(  rownames(ddat) , tx2gene$ensembl_gene_id)]
    # ddat <- aggregate( .~Gene.symbol , data=ddat, mean)
    # rownames(ddat) <- ddat$Gene.symbol
    # ddat <- ddat[,-1]
    eNOS.strepto <- Seurat::CreateSeuratObject( ddat , project = "eNOS.strepto" , min.cells = 100, min.genes = 500 )
    eNOS.strepto <- NormalizeData(eNOS.strepto)
    eNOS.strepto <- FindVariableFeatures(eNOS.strepto, verbose = T )
    eNOS.strepto <- ScaleData( eNOS.strepto ) 
    
    eNOS.strepto <- RunPCA( eNOS.strepto, verbose = T )
    ElbowPlot(eNOS.strepto, ndims = 50)
    eNOS.strepto <- FindNeighbors( eNOS.strepto, dims = 1:20)
    eNOS.strepto <- FindClusters( eNOS.strepto, resolution = 0.5)
    eNOS.strepto <- RunUMAP( eNOS.strepto, dims = 1:20 )
    
    ## add  group
    group <- as.factor( c(rep("strptzt",400), rep("wt",400), rep("strptzt",400), rep("wt",400) ) )
    group <- group[which(cellqual==1)]
    names( group ) <- colnames( eNOS.strepto )
    eNOS.strepto <- AddMetaData(object = eNOS.strepto, metadata = group, col.name = "group")
    # plot 
    p1 <- DimPlot( eNOS.strepto, reduction = "umap", label = TRUE)
    p2 <- DimPlot( eNOS.strepto, reduction = "umap", group.by ="group")
    
    p3 <- Seurat::FeaturePlot( eNOS.strepto , order=T , min.cutoff = 1, features = "WT1" , 
                               label = T, cols=c("lightgray","red")) 
    cowplot::plot_grid(p1, p2, p3, ncol = 3)
    
    # analyze podocytes
    eNOS.strepto_podo <- subset(eNOS.strepto, idents = 2)
    DefaultAssay( eNOS.strepto_podo ) <- "RNA"
    eNOS.strepto_podo <- NormalizeData(eNOS.strepto_podo)
    eNOS.strepto_podo <- FindVariableFeatures(eNOS.strepto_podo, verbose = T )
    eNOS.strepto_podo <- ScaleData( eNOS.strepto_podo ) 
    
    eNOS.strepto_podo <- RunPCA( eNOS.strepto_podo, verbose = T )
    ElbowPlot(eNOS.strepto_podo, ndims = 50)
    eNOS.strepto_podo <- FindNeighbors( eNOS.strepto_podo, dims = 1:30)
    eNOS.strepto_podo <- FindClusters( eNOS.strepto_podo, resolution = 0.5)
    eNOS.strepto_podo <- RunUMAP( eNOS.strepto_podo, dims = 1:30 )
    DimPlot( eNOS.strepto_podo, reduction = "umap", label = TRUE)
    Seurat::FeaturePlot( eNOS.strepto_podo , order=T , min.cutoff = 0, features = "Wt1" , 
                         label = T, cols=c("lightgray","red")) 
    Seurat::FeaturePlot( eNOS.strepto_podo , order=T , min.cutoff = 0, 
                         features = "Nos3" , 
                         label = T, cols=c("lightgray","red")) 
    
  }
  }

#### find sets of Universal FSGS marker genes #### 
  {
  #### Combine MA, bulk, sc RNAseq to get universal FSGS markers
  {
    ### make correlation plot
    {
      older <- FSGSexpr11[[1]][,c(1,2,6,9,10,11 )]
      ### load sc, analysis in GSE146912_scAnalysis.r script
      scDR<- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/scRNAseq_GSE146912/scRNAseq_GSE146912_DEtab.rda")
      sc_deDF <- scDR[[1]][,colnames(scDR[[1]])[c(F, T)]]
      ### load MA
      FSGS_MA_DE <-  readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA_lfcRank.rda")
      FSGS_MA_DEdf <- FSGS_MA_DE[[2]]
      FSGS_MA_DEdf <- FSGS_MA_DEdf[,!colnames(FSGS_MA_DEdf)%in%
                                     c("GSE117571__Foxc1.2_primPodo_MA","GSE106828__Tsc2_primPodo_MA" )]
      ### load bulk
      bulkDE <-  readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk_RankLfc.rda")
      lfc_bulkDF <- bulkDE[[2]]
      lfc_bulkDF <- lfc_bulkDF[,!colnames(lfc_bulkDF)%in%
                                 c("GSE117331__NPHS2cre_gloms_bulk" )]
      # colnames(older) <-  c( "KFO.Nphs2.hetdel_bulk", "KFO.adriomycin_bulk" , "GSE127235_eNOS.strep._podo_scRNAseq", "KFO.Wt1.hetdel_4w_bulk" ,
      #                        "KFO.Wt1.hetdel_12w_bulk" , "KFO.Wt1.hetdel_scRNAseq")
      
      gg2 <- Reduce(union, lapply(list(rownames(FSGS_MA_DEdf),
                                       rownames(older),
                                       rownames(lfc_bulkDF),
                                       rownames(sc_deDF) ), toupper ) )
      FSGS_XX <- cbind( older[match(gg2, toupper(rownames(older)) ),],
                        FSGS_MA_DEdf[match(gg2, toupper(rownames(FSGS_MA_DEdf))),] ,
                        lfc_bulkDF[match(gg2, toupper(rownames(lfc_bulkDF))),] ,
                        sc_deDF[match(gg2, toupper(rownames(sc_deDF)) ), ] )
      rownames(FSGS_XX) <- gg2
      
      ## plot corr
      colnames(FSGS_XX) <- sub( "podoDE","podo_sc",colnames(FSGS_XX) )
      colnames(FSGS_XX)[3] <- "GSE127235__eNOS.streptozt_podo_sc"
      XXX <- as.matrix( FSGS_XX )
      XXX[XXX <= -Inf] <- min(XXX[is.finite(XXX)], na.rm = T)
      XXX[XXX >= Inf] <- max(XXX[is.finite(XXX)], na.rm = T)
      #
      XXX <- cor( XXX , use = "complete", method = "spearman" )
      
      library("corrplot")
      corrplot::corrplot( XXX , hclust.method = "ward.D", order="hclust", tl.cex=0.8, cl.cex=0.8, tl.col="black"  )
      
    }
    
    ### combine ranks
    {
      ### 11 studies used previously 
      FSGSexpr11 <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/FSGSexpr11.rda")
      older <- FSGSexpr11[[1]][,c(1,2,5,6,9:11 )]
      older_rank <- FSGSexpr11[[2]][,c(1,2,5,6,9:11 )]
      older_rank <- apply(older_rank, 2, rank, na.last="keep")
      
      ### load sc, analysis in GSE146912_scAnalysis.r script
      scDR<- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/scRNAseq_GSE146912/scRNAseq_GSE146912_DEtab.rda")
      sc_deDF <- scDR[[1]][,colnames(scDR[[1]])[c(F, T)]]
      sc_rankDF <- scDR[[2]][,colnames(scDR[[2]])[c(F, T)]]
      sc_rankDF <- apply(sc_rankDF, 2, rank, na.last="keep") 
      
      ### load MA
      FSGS_MA_DE <-  readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA_lfcRank15.12.20.rda")
      FSGS_MA_DEdf <- FSGS_MA_DE[[2]]
      FSGS_MA_DEdf <- FSGS_MA_DEdf[,!colnames(FSGS_MA_DEdf)%in%
                                     c("GSE117571__Foxc1.2_primPodo_MA","GSE106828__Tsc2_primPodo_MA" ,
                                       "GSE18358__DenDrashWt1_gloms_MA","GSE969__antiGBM_kidney_MA",
                                       "GSE47606__Dendrin_gloms_MA")]
      FSGS_MA_pvalDF <- FSGS_MA_DE[[1]]
      FSGS_MA_pvalDF <- FSGS_MA_pvalDF[,!colnames(FSGS_MA_pvalDF)%in%
                                         c("GSE117571__Foxc1.2_primPodo_MA","GSE106828__Tsc2_primPodo_MA" ,
                                           "GSE18358__DenDrashWt1_gloms_MA","GSE969__antiGBM_kidney_MA",
                                           "GSE47606__Dendrin_gloms_MA")] 
      ### load bulk
      bulkDE <-  readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/bulk_DE15.12.20.rda")
      lfc_bulkDF <- bulkDE[[2]]
      lfc_bulkDF <- lfc_bulkDF[,!colnames(lfc_bulkDF)%in%
                                 c("GSE117331__NPHS2cre_gloms_bulk" )]
      rank_bulkDF <- bulkDE[[1]]
      rank_bulkDF <- rank_bulkDF[,!colnames( rank_bulkDF ) %in% 
                                   c( "GSE117331__NPHS2cre_gloms_bulk" )]
      
      gg2 <- Reduce(union, lapply(list(rownames(FSGS_MA_DEdf),
                                       rownames(older),
                                       rownames(lfc_bulkDF),
                                       rownames(sc_deDF) ), toupper ) )
      
      FSGS_XX <- cbind( older[match(gg2, toupper(rownames(older)) ),],
                        FSGS_MA_DEdf[match(gg2, toupper(rownames(FSGS_MA_DEdf))),] ,
                        lfc_bulkDF[match(gg2, toupper(rownames(lfc_bulkDF))),] ,
                        sc_deDF[match(gg2, toupper(rownames(sc_deDF)) ), ] )
      
      ### combine ranks
    
      # colnames(older_rank) <- c( "KFO.Nphs2.hetdel_bulk", "KFO.adriomycin_bulk" , "GSE127235_eNOS.strep._podo_scRNAseq")
      FSGS_rank <-  cbind( older_rank[match(gg2, toupper(rownames(older_rank))),],
                           FSGS_MA_pvalDF[match(gg2, toupper(rownames(FSGS_MA_pvalDF))),] , 
                           rank_bulkDF[match(gg2, toupper(rownames(rank_bulkDF))),] ,
                           sc_rankDF[match(gg2, toupper(rownames(sc_rankDF))), ] )
      FSGS_XX <- as.matrix(FSGS_XX)
      FSGS_rank <- as.matrix(FSGS_rank)
      rownames(FSGS_XX) <- rownames(FSGS_rank) <- gg2
      FSGS_XX[FSGS_XX <= -Inf] <- min(FSGS_XX[is.finite(FSGS_XX)], na.rm = T)
      FSGS_XX[FSGS_XX >= Inf] <- max(FSGS_XX[is.finite(FSGS_XX)], na.rm = T)
      
      saveRDS( list(lfc=FSGS_XX, rank=FSGS_rank), file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/CrossValidation/38fsgs_data.rda")
      
      #### mean Rank 
      # choose genes 
      FSGS_rankM <- FSGS_rank[rowSums(!is.na(FSGS_rank))>round(ncol(FSGS_rank)*0.75),]
      FSGS_rankMean <- rowMeans( FSGS_rankM, na.rm = T )
      # FSGS_rankMean <- FSGS_rankMean[order(FSGS_rankMean)]
      # choose same direction of expression in at least 75% of datasets
      FSGS_lfcM <- FSGS_XX[rownames(FSGS_rankM),]
      FSGS_lfcM <- ifelse( rowSums(FSGS_lfcM<0, na.rm = T)/rowSums(!is.na(FSGS_lfcM)) > 0.75 | 
                             rowSums(FSGS_lfcM>0, na.rm = T)/rowSums(!is.na(FSGS_lfcM)) > 0.75 , T, F )
      
      topFSGS26 <- FSGS_rankMean[ FSGS_lfcM ]
      topFSGS26 <- as.data.frame(topFSGS26)
      topFSGS26$meanLFC <- rowMeans( FSGS_XX[rownames(FSGS_rankM),][FSGS_lfcM,] , na.rm = T )
      
      # write.table(FSGS26_0.8 , sep = "\t","/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/FSGS23_markers0.8_21stud.tsv")
      
      FSGS26_0.8_top100 <- topFSGS26[order(topFSGS26$topFSGS26),][1:100,]
      FSGS26_0.8_top50 <- topFSGS26[order(topFSGS26$topFSGS26),][1:50,]
      FSGS26_0.8_top20 <- topFSGS26[order(topFSGS26$topFSGS26),][1:20,]
    }
    
    
  }
  
  ### train a classifier 
  {
    library(glm)
    mydata <- scPodoWt1[ rownames(scPodoWt1) %in% rownames( FSGS26_0.8_top100),]
    
    mydata <- cbind.data.frame( gt= as.factor ( c(rep("wt",587),rep("mut", ncol(scPodoWt1)-587)) ) ,
                                t(mydata))
    scPodo_markLogit <- glm( gt ~ . , data = mydata, family = binomial(link = "logit"), maxit = 100)
    scPodo_markLogit <- as.data.frame( coef(summary(scPodo_markLogit)) )
    scPodo_markLogit$padj <- p.adjust(scPodo_markLogit[,4], method = "fdr" )
    scPodo_markLogit <- scPodo_markLogit[which(scPodo_markLogit$padj<0.1),]
    XX <- rownames(scPodo_markLogit[which(scPodo_markLogit$padj<0.1),])
    
    ## up
    mydata <- scPodoWt1[ rownames(scPodoWt1) %in% rownames( FSGS26_0.8[which(FSGS26_0.8$meanLFC>0),]),]
    
    mydata <- cbind.data.frame( gt= as.factor ( c(rep("wt",587),rep("mut", ncol(scPodoWt1)-587)) ) ,
                                t(mydata))
    scPodo_markLogit <- glm( gt ~ . , data = mydata, family = binomial(link = "logit"), maxit = 100)
    scPodo_markLogit <- as.data.frame( coef(summary(scPodo_markLogit)) )
    scPodo_markLogit$padj <- p.adjust(scPodo_markLogit[,4], method = "fdr" )
    scPodo_markLogit <- scPodo_markLogit[which(scPodo_markLogit$padj<0.1),]
    XXup <- rownames(scPodo_markLogit[which(scPodo_markLogit$padj<0.1),])[-1]
    
    ## down
    mydata <- scPodoWt1[ rownames(scPodoWt1) %in% rownames( FSGS26_0.8[which(FSGS26_0.8$meanLFC<0),]),]
    
    mydata <- cbind.data.frame( gt= as.factor ( c(rep("wt",587),rep("mut", ncol(scPodoWt1)-587)) ) ,
                                t(mydata))
    scPodo_markLogit <- glm( gt ~ . , data = mydata, family = binomial(link = "logit"), maxit = 100)
    scPodo_markLogit <- as.data.frame( coef(summary(scPodo_markLogit)) )
    scPodo_markLogit$padj <- p.adjust(scPodo_markLogit[,4], method = "fdr" )
    scPodo_markLogit <- scPodo_markLogit[which(scPodo_markLogit$padj<0.1),]
    XXdown <- rownames(scPodo_markLogit[which(scPodo_markLogit$padj<0.1),])[-1]
    
  }
  
}


#### validate the classic disease score on Wt1 studies  #### 
  {
  # read marker genes
  FSGS26_0.8 <- read.table("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/FSGS23_markers0.8_21stud.tsv")
  FSGS26_0.8_top100 <- FSGS26_0.8[order(FSGS26_0.8$topFSGS26),][1:100,]
  FSGS26_0.8_top50 <- FSGS26_0.8[order(FSGS26_0.8$topFSGS26),][1:50,]
  FSGS26_0.8_top20 <- FSGS26_0.8[order(FSGS26_0.8$topFSGS26),][1:20,]
  FSGS26_0.8_top200 <- FSGS26_0.8[order(FSGS26_0.8$topFSGS26),][1:200,]
  
  ## wt1 sc RNAseq
  {
    scPodoWt1 <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/scPodo_libNorm.rda")
    rownames(scPodoWt1) <- toupper(rownames(scPodoWt1))
    
    scPodoWt1_disScore <- disease_score(
      scRNAseq = scPodoWt1 , wt_cells=1:587 , 
      marker_gene_list=rownames(FSGS26_0.8_top20) , 
      marker_gene_LFC=FSGS26_0.8_top20$meanLFC , 
      mainPlot="disease score based on 100 Universal markers (23stud)"
    )
    
    scPodoWt1_disScore20 <- cbind.data.frame( score=scPodoWt1_disScore20, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    scPodoWt1_disScore50 <- cbind.data.frame( score=scPodoWt1_disScore50, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    scPodoWt1_disScore100 <- cbind.data.frame( score=scPodoWt1_disScore100, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    scPodoWt1_disScoreLogit <- cbind.data.frame( score=scPodoWt1_disScoreLogit, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    
    scPodoWt1_disScore20x <- cbind.data.frame( score=scPodoWt1_disScore20x, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    scPodoWt1_disScore50x <- cbind.data.frame( score=scPodoWt1_disScore50x, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    scPodoWt1_disScore100x <- cbind.data.frame( score=scPodoWt1_disScore100x, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
    
    library( ggplot2 )
    library(ggthemes)
    g1 <- ggplot2::ggplot(scPodoWt1_disScore20, aes( x=score, color=gt)) +
      scale_color_colorblind() +  geom_density(size=1.5) + theme_bw() + theme(legend.position = "none") +
      geom_vline(data=ddply(scPodoWt1_disScore20, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
    g2 <-ggplot2::ggplot(scPodoWt1_disScore50, aes( x=score, color=gt)) + scale_color_colorblind() + 
      geom_density(size=1.5) + theme_bw() + theme(legend.position = "none") +
      geom_vline(data=ddply(scPodoWt1_disScore50, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
    g3 <-ggplot2::ggplot(scPodoWt1_disScore100, aes( x=score, color=gt)) + 
      scale_color_colorblind() + geom_density(size=1.5) + theme_bw() + theme(legend.position = "none")
    geom_vline(data=ddply(scPodoWt1_disScore100, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
    g4 <- ggplot2::ggplot(scPodoWt1_disScoreLogit, aes( x=score, color=gt)) + 
      scale_color_colorblind() + geom_density(size=1.5) + theme_bw() +
      geom_vline(data=ddply(scPodoWt1_disScoreLogit, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
    
    cowplot::plot_grid (g1, g2, g3, g4, nrow = 2, labels=c("top 20 markers","top 50 markers","top 100 markers","logit 15 markers"))     
    
  }
  
  ### Wt1 MA denDrash
  {
    Wt1DenDrash_expr <- read.table( fill = T, header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA/GSE18358_series_matrix.txt")
    Wt1DenDrash_DE <- read.table( fill = T, header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DATA/MA/DE/GSE18358__DenDrashWt1_gloms.txt")
    Wt1DenDrash_expr$Gene.symbol <- Wt1DenDrash_DE$Gene.symbol[ match( Wt1DenDrash_expr$ID_REF , Wt1DenDrash_DE$ID)]
    Wt1DenDrash_expr <- Wt1DenDrash_expr[,-1]
    Wt1DenDrash_expr <- aggregate( .~Gene.symbol , data= Wt1DenDrash_expr, mean)
    Wt1DenDrash_expr$Gene.symbol <- toupper(Wt1DenDrash_expr$Gene.symbol)
    Wt1DenDrash_expr<- Wt1DenDrash_expr[-1,] # remove empty gene
    rownames(Wt1DenDrash_expr) <- Wt1DenDrash_expr$Gene.symbol
    Wt1DenDrash_expr <- Wt1DenDrash_expr[,-1]
    
    MAdendrash_disScore20 <- disease_score(
      scRNAseq = Wt1DenDrash_expr , wt_cells=1:5 , 
      marker_gene_list=rownames(FSGS26_0.8_top20) , 
      marker_gene_LFC=FSGS26_0.8_top20$meanLFC , 
      mainPlot="disease score based on 20 Universal markers (23stud)"
    )
    
    MAdendrash_disScore20 <- cbind.data.frame( score=MAdendrash_disScore20, gt=c( rep("wild.type", 5), rep("Wt1dd", 5)))
    MAdendrash_disScore50 <- cbind.data.frame( score=MAdendrash_disScore50, gt=c( rep("wild.type", 5), rep("Wt1dd", 5)))
    MAdendrash_disScore100 <- cbind.data.frame( score=MAdendrash_disScore100, gt=c( rep("wild.type", 5), rep("Wt1dd", 5)))
    MAdendrash_disScoreLogit <- cbind.data.frame( score=MAdendrash_disScoreLogit, gt=c( rep("wild.type", 5), rep("Wt1dd", 5)))
    
    g1 <- ggplot2::ggplot( MAdendrash_disScore20, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    g2 <- ggplot2::ggplot( MAdendrash_disScore50, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    g3 <- ggplot2::ggplot( MAdendrash_disScore100, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    g4 <- ggplot2::ggplot( MAdendrash_disScoreLogit, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    
    cowplot::plot_grid (g1, g2, g3, g4, nrow = 2, labels=c("top 20 markers","top 50 markers","top 100 markers","15 logit markers"))     
    
    
  }
  
  ### Wt1 het.del bulk (KFO)
  {
    bulkWt1hd_expr <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/expr_IntExon_glvl.rda")
    bulkWt1hd_expr$Gene.symbol <- tx2gene$external_gene_name[ match( rownames(bulkWt1hd_expr) , tx2gene$ensembl_gene_id)]
    bulkWt1hd_expr <- aggregate( .~Gene.symbol , data= bulkWt1hd_expr, mean)
    rownames(bulkWt1hd_expr) <- toupper( bulkWt1hd_expr$Gene.symbol )
    # bulkWt1hd_expr <- bulkWt1hd_expr[,-1]
    # bulkWt1hd_expr <- log( bulkWt1hd_expr[,-1] +1)
    
    bulkWt1hd_disScore20 <- disease_score(
      scRNAseq = bulkWt1hd_expr, wt_cells=c(2,3,5,7,8,12) , 
      marker_gene_list=rownames(FSGS26_0.8_top20) , 
      marker_gene_LFC=FSGS26_0.8_top20$meanLFC , 
      mainPlot="disease score based on 20 Universal markers (23stud)"
    )
    
    
    
    bulkWt1hd_disScore20 <- cbind.data.frame( score=bulkWt1hd_disScore20, gt=c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12"))
    bulkWt1hd_disScore50 <- cbind.data.frame( score=bulkWt1hd_disScore50, gt=c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12"))
    bulkWt1hd_disScore100 <- cbind.data.frame( score=bulkWt1hd_disScore100, gt=c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12"))
    bulkWt1hd_disScoreLogit <- cbind.data.frame( score=bulkWt1hd_disScoreLogit, gt=c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12"))
    
    g1 <- ggplot2::ggplot( bulkWt1hd_disScore20, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    g2 <- ggplot2::ggplot( bulkWt1hd_disScore50, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    g3 <- ggplot2::ggplot( bulkWt1hd_disScore100, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    g4 <- ggplot2::ggplot( bulkWt1hd_disScoreLogit, aes( y=score, x=gt, color=gt)) +
      scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
      stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
    
    cowplot::plot_grid (g1, g2, g3,    g4, nrow = 2, labels=c("top 20 markers","top 50 markers","top 100 markers","logit 15 markers"))     
    
  }
  
  
}

#### rank based disease score #### 
  {
  
  ### prepare the data
  # load expression matrix
  exprMatrix <- scPodoWt1
  # exprMatrix[exprMatrix==0] <- NA
  exprMatrixMA <- Wt1DenDrash_expr
  # exprMatrixMA[exprMatrixMA==0] <- NA
  exprMatrixBULK <- bulkWt1hd_expr
  # exprMatrixBULK[exprMatrixBULK==0] <- NA
  
  
  # calculate disease score using ranks
  cells_rankings_scoreUP <- apply(exprMatrix , 2, rank, ties.method = "last")
  cells_rankings_scoreDOWN <- apply(exprMatrix*-1 , 2, rank, ties.method = "last")
  cells_rankings_score <- scale(colMeans(cells_rankings_scoreUP[rownames(cells_rankings_scoreUP) %in% rownames( FSGS26_0.8_top50[FSGS26_0.8_top50$meanLFC>0,]),], na.rm = F) +
                                  colMeans(cells_rankings_scoreDOWN[rownames(cells_rankings_scoreDOWN) %in% rownames( FSGS26_0.8_top50[FSGS26_0.8_top50$meanLFC<0,]),], na.rm = F) )
  
  cells_rankingsMA_scoreUP <- apply(exprMatrixMA , 2, rank, ties.method = "last")
  cells_rankingsMA_scoreDOWN <- apply(exprMatrixMA*-1 , 2, rank, ties.method = "last")
  cells_rankingsMA_score <- scale(colMeans(cells_rankingsMA_scoreUP[rownames(cells_rankingsMA_scoreUP) %in% rownames( FSGS26_0.8_top50[FSGS26_0.8_top50$meanLFC>0,]),], na.rm = F) + 
                                    colMeans(cells_rankingsMA_scoreDOWN[rownames(cells_rankingsMA_scoreDOWN) %in% rownames( FSGS26_0.8_top50[FSGS26_0.8_top50$meanLFC<0,]),], na.rm = F) )
  
  cells_rankingsBULK_scoreUP <- apply(exprMatrixBULK , 2, rank, ties.method = "last")
  cells_rankingsBULK_scoreDOWN <- apply(exprMatrixBULK *-1 , 2, rank, ties.method = "last")
  cells_rankingsBULK_score <- scale(colMeans(cells_rankingsBULK_scoreUP[rownames(cells_rankingsBULK_scoreUP) %in% rownames( FSGS26_0.8_top50[FSGS26_0.8_top50$meanLFC>0,]),], na.rm = F) + 
                                      colMeans(cells_rankingsBULK_scoreDOWN[rownames(cells_rankingsBULK_scoreDOWN) %in% rownames( FSGS26_0.8_top50[FSGS26_0.8_top50$meanLFC<0,]),], na.rm = F) )
  
  ### plot separation
  AUCell_top50_CC <-   cbind.data.frame( score=cells_rankings_score, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
  AUCell_top50ma_CC <-   cbind.data.frame( score= cells_rankingsMA_score , gt=c( rep("wild.type", 5), rep("Wt1dd", 5)) )
  AUCell_top50bulk_CC <-   cbind.data.frame( score=cells_rankingsBULK_score , gt= c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12"))
  
  
  # AUCell_Logit_up <-   cbind.data.frame( score=cells_AUC_mat[4,], gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
  # AUCell_Logit_down <-   cbind.data.frame( score=cells_AUC_mat[5,], gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
  
  g1 <-ggplot2::ggplot(AUCell_top50_CC, aes( x=score, color=gt)) + 
    scale_color_colorblind() + geom_density(size=1.5) + theme_bw() + theme(legend.position = "none") +
    geom_vline(data=ddply(AUCell_top50_CC, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
  g2 <- ggplot2::ggplot( AUCell_top50ma_CC, aes( y=score, x=gt, color=gt)) +
    scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
  g3 <- ggplot2::ggplot( AUCell_top50bulk_CC, aes( y=score, x=gt, color=gt)) +
    scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
  
  cowplot::plot_grid (g1, g2, g3, nrow = 1, labels=c("top 50 markers, scRNAseq","top 50 markers, MA", "top 50 markers, bulk" )) 

  
  ### plot relation between the score and number of detected cells in scRNAseq
  plot( AUCell_top50_CC$score , colSums(exprMatrix>0), xlab="scaled rank-based score",
        ylab = "number of detected genes")
  mtext("Spearman correlation -0.013")
  
}

#### AUCell rank based score ####
  {
  library(AUCell)
  library(GSEABase)
    FSGS26_0.8_top50_UP <- rownames( FSGS26_0.8_top50[which(FSGS26_0.8_top50$meanLFC > 0) ,] )
    FSGS26_0.8_top50_DOWN <- rownames( FSGS26_0.8_top50[which(FSGS26_0.8_top50$meanLFC < 0) ,] )
    FSGS26_0.8_top200_UP <- rownames( FSGS26_0.8_top200[which(FSGS26_0.8_top200$meanLFC > 0) ,] )
    FSGS26_0.8_top200_DOWN <- rownames( FSGS26_0.8_top200[which(FSGS26_0.8_top200$meanLFC < 0) ,] )
    FSGS26_0.8_top20_UP <- rownames( FSGS26_0.8_top20[which(FSGS26_0.8_top20$meanLFC > 0) ,] )
    FSGS26_0.8_top20_DOWN <- rownames( FSGS26_0.8_top20[which(FSGS26_0.8_top20$meanLFC < 0) ,] )
    FSGS26_0.8_top100_UP <- rownames( FSGS26_0.8_top100[which(FSGS26_0.8_top100$meanLFC > 0) ,] )
    FSGS26_0.8_top100_DOWN <- rownames( FSGS26_0.8_top100[which(FSGS26_0.8_top100$meanLFC < 0) ,] )
    
    geneSets <-  GeneSetCollection( c(
      GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50)"),
      GeneSet( ( FSGS26_0.8_top50_UP), setName="top50 UP markers") , 
      GeneSet( ( FSGS26_0.8_top50_DOWN), setName="top50 DOWN markers"),
      GeneSet( ( FSGS26_0.8_top20_UP), setName="top20 UP markers") , 
      GeneSet( ( FSGS26_0.8_top20_DOWN), setName="top20 DOWN markers"),
      GeneSet( ( FSGS26_0.8_top200_UP), setName="top200 UP markers") , 
      GeneSet( ( FSGS26_0.8_top200_DOWN), setName="top200 DOWN markers") , 
      GeneSet( ( FSGS26_0.8_top100_UP), setName="top100 UP markers") , 
      GeneSet( ( FSGS26_0.8_top100_DOWN), setName="top100 DOWN markers"))) 
    
  # 1. Build gene-expression rankings for each cell  
  cells_rankings <- AUCell_buildRankings( exprMatrix, nCores=1, plotStats=TRUE)
  cells_rankingsMA <- AUCell_buildRankings( as.matrix( exprMatrixMA ), nCores=1, plotStats=F)
  cells_rankingsBULK <- AUCell_buildRankings( as.matrix( exprMatrixBULK ), nCores=1, plotStats=F)
  
  # 2. Calculate enrichment for the gene signatures (AUC)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = 500 )
  cells_AUCma <- AUCell_calcAUC(geneSets, cells_rankingsMA )
  cells_AUCbulk <- AUCell_calcAUC(geneSets, cells_rankingsBULK )
  
  cells_AUC_mat <- getAUC( cells_AUC)
  cells_AUCma_mat <- getAUC( cells_AUCma )
  cells_AUCbulk_mat <- getAUC( cells_AUCbulk )
  
  # make a plot 
  # XX <-   cbind.data.frame( score=( cells_AUC_mat[4,]-cells_AUC_mat[5,] ), gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
  XX1 <-   cbind.data.frame( score=( cells_AUC_mat[2,]*12-cells_AUC_mat[3,]*19 ), gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
  XX2 <-   cbind.data.frame( score=( cells_AUCma_mat[2,]*23-cells_AUCma_mat[3,]*26 ), gt=c( rep("wild.type", 5), rep("Wt1dd",5)))
  XX3 <-   cbind.data.frame( score=( cells_AUCbulk_mat[2,]*24-cells_AUCbulk_mat[3,]*26 ), gt= c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12") )
  
  g1 <- ggplot2::ggplot(XX1, aes( x=score, color=gt)) + 
    scale_color_colorblind() + geom_density(size=1.5) + theme_bw() + theme(legend.position = "none") +
    geom_vline(data=ddply(XX1, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
  g2 <- ggplot2::ggplot( XX2, aes( y=score, x=gt, color=gt)) +
    scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
  g3 <- ggplot2::ggplot( XX3, aes( y=score, x=gt, color=gt)) +
    scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
  
  cowplot::plot_grid (g1, g2, g3, nrow = 1, labels=c("top 50 markers, scRNAseq","top 50 markers, MA", "top 50 markers, bulk" )) 
  
  }


#### correlation based disease score #### 
  {
  
  x <- FSGS26_0.8[which(rownames(FSGS26_0.8) %in% rownames(scPodoWt1) ),]
  exprReg_sc <-  cbind.data.frame( marks= x$meanLFC ,
                                    log2( (scPodoWt1[ match( rownames(x),rownames(scPodoWt1) ),]+1) /
                                      (rowMeans( scPodoWt1[ match( rownames(x),rownames(scPodoWt1) ),1:587])+1 ) ) )
  scoreReg_sc <- sapply( 2:ncol(exprReg_sc), function (x){
    cor( exprReg_sc[which(rownames(exprReg_sc) %in% rownames(FSGS26_0.8_top20)),1 ] , 
         exprReg_sc[which(rownames(exprReg_sc) %in% rownames(FSGS26_0.8_top20)),x ] ,
         method = "spearman") 
  })
  scoreReg_sc <- cbind.data.frame( score=scoreReg_sc, gt=c( rep("wild.type", 578), rep("Wt1KO", ncol(scPodoWt1)-578)))
  
    
   exprReg_ma <- ( cbind.data.frame( marks= FSGS26_0.8$meanLFC ,
                                   Wt1DenDrash_expr[ match( rownames(FSGS26_0.8),rownames(Wt1DenDrash_expr) ),] -
                                   rowMeans( Wt1DenDrash_expr[ match( rownames(FSGS26_0.8),rownames(Wt1DenDrash_expr) ),1:5]) ) )
   scoreReg_ma <- sapply( 2:ncol(exprReg_ma), function (x){
     cor( exprReg_ma[which(rownames(exprReg_ma) %in% rownames(FSGS26_0.8_top20)),1 ] , 
          exprReg_ma[which(rownames(exprReg_ma) %in% rownames(FSGS26_0.8_top20)),x ] ,
          method = "spearman") 
     })
   scoreReg_ma <-   cbind.data.frame( score=scoreReg_ma, gt=c( rep("wild.type", 5), rep("Wt1dd",5)))

   
   exprReg_bulk <- cbind.data.frame( marks= FSGS26_0.8$meanLFC ,
                                     log2( bulkWt1hd_expr[ match( rownames(FSGS26_0.8),rownames(bulkWt1hd_expr) ), ] /
                                     rowMeans( bulkWt1hd_expr[ match( rownames(FSGS26_0.8),rownames(bulkWt1hd_expr) ),c(2,3,5,7,8,12)]) ))
   scoreReg_bulk <- sapply( 2:ncol(exprReg_bulk), function (x){
     cor( exprReg_bulk[which(rownames(exprReg_bulk) %in% rownames(FSGS26_0.8_top20)),1 ] , 
          exprReg_bulk[which(rownames(exprReg_bulk) %in% rownames(FSGS26_0.8_top20)),x ] ,
          method = "spearman") 
   })
   scoreReg_bulk <- cbind.data.frame( score=scoreReg_bulk, gt= c( "Wt1hd_4","wt_4","wt_4","Wt1hd_4","wt_4","Wt1hd_4","wt_12","wt_12","Wt1hd_12","Wt1hd_12","Wt1hd_12","wt_12") )
   ### make a plot 
   g1 <-ggplot2::ggplot(scoreReg_sc, aes( x=score, color=gt)) + 
     scale_color_colorblind() + geom_density(size=1.5) + theme_bw() + theme(legend.position = "none") +
     geom_vline(data=ddply(scoreReg_sc, "gt", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=gt), linetype="dashed")
   g2 <- ggplot2::ggplot( scoreReg_ma, aes( y=score, x=gt, color=gt)) +
     scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
     stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
   g3 <- ggplot2::ggplot( scoreReg_bulk, aes( y=score, x=gt, color=gt)) +
     scale_color_colorblind() +  geom_jitter(width = 0.1) + theme_bw() + theme(legend.position = "none") +
     stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom = "crossbar", width = .2, color = "red")  
   
   cowplot::plot_grid (g1, g2, g3, nrow = 1, labels=c("top 50 markers, scRNAseq","top 50 markers, MA", "top 50 markers, bulk" )) 
   
  }

