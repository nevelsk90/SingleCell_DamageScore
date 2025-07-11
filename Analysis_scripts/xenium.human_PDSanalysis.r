############ xenium human data, may 2025 ############
library( Seurat)
library(SeuratObject)
library(SeuratDisk)
library( ggplot2)
library(ggpubr)
library(ggthemes)
library( dplyr)
  library( sf )

mallinfo::malloc.trim()
gc()
mallinfo::malloc.trim()
gc()

#### load functions
source("/data/user/tpadvits/PROJECTS/PDS/AUCell_script.r")
# load the podo damage signature
DS_all <- read.table( header = T, "/data/user/tpadvits/PROJECTS/PDS/DS_all.20.09.2023.tsv")

               

#### create Seurat objects from all xenium samples and save them ####
### create Seurat objects
# list Xenium folders
listPaths <- list.dirs( recursive=F , "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/")

lapply( seq( listPaths) , function(ii)
{
  
  print(ii)
  
  data <- ReadXenium( listPaths[[ii]],  type =c('centroids', 'segmentations',"nucleus_segmentations"))
  
  assay <- "Xenium"
  
  # semi manually create seurat objects https://github.com/satijalab/seurat/issues/7491
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], 
                                   assay = assay, project = 
                                     gsub("output-XETG00229__|__20250514__090126","",basename(listPaths)[ii]) )
  xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  # add cordinates
  segmentations.data <- list(
    "segmentations" = CreateSegmentation(data$segmentations),
    "centroids" = CreateCentroids(data$centroids),
    "nucleus_segmentations" = CreateSegmentation(data$nucleus_segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentations", "centroids","nucleus_segmentations"),
    molecules = data$microns,
    assay = assay
  )
  xenium.obj[["fov"]] <- coords
  

  xenium.obj <- subset(xenium.obj , subset = nCount_Xenium > 0)

  
  saveRDS(xenium.obj, file =paste0( "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/",
                                    gsub("output-XETG00229__|__20250514__090126","",basename(listPaths)[ii]),"_Seurat.rda") )
  # return(xenium.obj)
  
})


### plot nCount per dataset
ggl <- lapply( seq(list.seu ) , function(ii){
  
  print(ii)
  
  xenium.obj <- readRDS(list.seu[[ii]])
  
  Idents(xenium.obj) <-  gsub("output-XETG00229__|__20250514__090126","", Idents(xenium.obj))
  
  
  
  ### make plots 
  VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
  
})
cowplot::plot_grid(plotlist = ggl)




#### segment and split scans in individual samples #### 
library(mapview)   # interactive map canvas
library(mapedit)   # lasso / polygon drawing
library( dplyr )
library(BiocParallel)
library( sf)
source("/data/user/tpadvits/PROJECTS/scripts/crop_seurat_v1.R")


#### read Xenium data and draw polygons
list.seu <- list.files(pattern = "Seurat.RDS",full.names = T,
                       path = "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/Seurat/")
# 0. Tell mapview to *leave coordinates alone*
mapviewOptions(native.crs = TRUE)


## read Xenium data
ii<-1
xenium.obj <- readRDS(list.seu[ii])

# 1. fetch centroid coordinates (low-res plot units) and make an sf POINT layer
coords <- GetTissueCoordinates(xenium.obj[["fov"]], which = "centroids")
cent_sf <- st_as_sf(coords,
                    # crs=3857,
                    coords = c("x", "y"),
                    crs=NA
                    )             # dummy CRS – we only care about geometry

# ######## #
# 2. launch an editable map and draw a polygon
poly_sf <- editMap(mapview(cent_sf, cex = 0.5), crs=NA)  # draw, double-click, “Done”

# 3. which centroids fall inside?
inside <- st_within(cent_sf, poly_sf$finished, sparse = FALSE)
cells_keep <- coords$cell[inside]

# 4. subset the original Seurat object
# xen_roi <- subset(xenium.obj, cells = cells_keep)

subset_coord <- data.frame( x= st_coordinates(poly_sf$drawn)[,1] ,
                            y= st_coordinates(poly_sf$drawn)[,2]  )


## crop custom to remove any tissue area, including molecules
# https://github.com/alikhuseynov/add-on_R/blob/develop/R/crop_seurat_v1.R
# alternatively use https://github.com/scOpenLab/geojson_seurat_cropper/blob/main/SeuratCropper.R ,
# with pyxel_size = 1 and change code object@images[[name]] -> object@images[["fov"]]
xen_roi3 <-  Crop_custom( x=coords ,
                          object = xenium.obj,
                      col_id = c("cell"),
                      xy_pts = subset_coord ,
                      c_hull_include = TRUE,
                      crop_molecules = TRUE,
                      BPPARAM = BiocParallel::SerialParam())

## add cell-type annotation
# # read Azimuth annotations for individual samples
# ctype.ll <- list.files(full.names = T, pattern = "azimuth_pred.tsv","/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/Seurat/Ind.Samples/Azimuth_ctype.pred/")[
#   1:16 ]

ctype.ann <- read.table(ctype.ll[[7]], header = T, sep = "\t", row.names = 1)

# add ctype to meta
xen_roi3@meta.data <- cbind( xen_roi3@meta.data ,
                                 ctype.ann[  match( rownames(xen_roi3@meta.data), rownames(ctype.ann)), ] )

# plot
ImageDimPlot(xen_roi3, fov = "fov", molecules = c("WT1"),
             coord.fixed = T , nmols = 20000, shuffle.cols = F,
             axes = TRUE, size = 5, 
             mols.size = 0.1, border.size = 0.1, cols =  c( "lightgrey",  "red"), 
             mols.cols = "green",  boundaries = "segmentations")


pathh="/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/Seurat/Ind.Samples/"
saveRDS(xen_roi3, file =paste0(pathh,  gsub("Seurat","Seurat__PID.01-00113..rep1",basename(list.seu)[ii])) )

#### do cell type annotation ####
list.seu2 <- list.files(pattern = "Seurat__PID", full.names = T, path = "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/Seurat/Ind.Samples")

### vln plot of nCount per sample
ggl <- lapply(list.seu2, function(fname)
  {
  xenium.obj<- readRDS(fname)
 
  VlnPlot(xenium.obj,  features ="nCount_Xenium", log = T , pt.size = 0)+
    theme(legend.position = "none")+
    xlab( gsub(".*Seurat__|.RDS","", basename(list.seu2)[ii]))
  
})


ggl0 <- cowplot::plot_grid( plotlist = ggl, nrow=2)

pdf( width = 20, height = 10, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/xenium.Human_nCount.IndSample_VlnPlot.pdf")
print(ggl0)
dev.off()


### cell type annotation 
# # prepare datasets to load in Azimuth  
# # run Azimuth https://app.azimuth.hubmapconsortium.org/app/human-kidney 
# lapply(list.seu2, function(fname)
#   {
#   xenium.obj<- readRDS(fname)
#   xenium_expr  <- CreateSeuratObject(counts=GetAssayData(xenium.obj, layer = 'counts' ))
#   xenium_expr[["RNA"]] <- as(object = xenium_expr[["RNA"]], Class = "Assay")
#   SaveH5Seurat(xenium_expr, overwrite = T, file = sub( "RDS", "h5seurat", fname))
#   
# })

###  plot results of Azimuth annotation add cell-type annotation to Seurat obj
ctype.ll <- list.files(full.names = T, pattern = "azimuth_pred.tsv","/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/Seurat/Ind.Samples/Azimuth_ctype.pred/")[
  1:16 ]

# read Azimuth annotations 
xenium_ctype  <- lapply( seq(ctype.ll ) , function(ii)
  {
  
  print(ii)
  
  ctype.ann <- read.table(ctype.ll[[ii]], header = T, sep = "\t", row.names = 1)
  # xenium.obj <- readRDS(list.seu3[[ii]])
  # 

  return(ctype.ann)
})
names(xenium_ctype) <- basename(ctype.ll)

# calculate percentage of each cell-type
xenium_ctype_prcntg <- Reduce( rbind,  lapply( seq(xenium_ctype), function(ii)
  {
  datt <- xenium_ctype[[ii]]
  datt$dataSet <-  sub(".tsv","",names(xenium_ctype)[ii])
  colnames(datt) <- c("predicted.annotation", "predicted.annotation.score", "mapping.score","dataSet")
  # gg1 <- ggplot2::ggplot(datt , aes(x=predicted.annotation.l2 , y=predicted.annotation.l2.score)) +
  #   geom_violin()  #
  datt <-datt[ datt$predicted.annotation %in%
                 names( table(datt$predicted.annotation))[
                   table(datt$predicted.annotation)>0 ] , ]
  prcntg <- data.frame(unlist(table(datt$predicted.annotation)))
  colnames(prcntg) <- c("ctype","nCells") 
  prcntg$Freq <-  (prcntg$nCells/sum(prcntg$nCells))*100
  prcntg$dataSet <- basename(ctype.ll)[ii]
  return(prcntg)
}))

xenium_ctype_prcntg$dataSet <- gsub("___azimuth_pred.tsv|.*Seurat__","",xenium_ctype_prcntg$dataSet)

### plot nCells and percentage
gg1<- ggplot2::ggplot(xenium_ctype_prcntg , aes(x=dataSet , y = Freq)) +
  # geom_bar(aes(y = after_stat(count)/sum(after_stat(count)))) + 
  geom_bar( stat = "identity" )+
  ylab("relative frequencies") + theme_minimal() +
  facet_wrap(vars( ctype ))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gg2<- ggplot2::ggplot(xenium_ctype_prcntg , aes(x=dataSet , y = nCells)) +
  # geom_bar(aes(y = after_stat(count)/sum(after_stat(count)))) + 
  geom_bar( stat = "identity" )+
  ylab("relative frequencies") + theme_minimal() +
  facet_wrap(vars( ctype ), scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggl <- cowplot::plot_grid( plotlist = list(gg1, gg2), labels =c("Percentage", "Number of cells" ))

pdf( width = 40, height = 20, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/ctype_Nmbr.Prct_braplot.pdf")
 print(ggl)
dev.off()

### add ctype annotation to Seurat objects, plot ctype
lapply(list.seu2, function(fname)
{
  xenium.obj<- readRDS(fname)

  ctype.anno <- read.table(  sep = "\t", header = T,
    sub( "Ind.Samples", "Ind.Samples/Azimuth_ctype.pred",
        sub(".RDS","___azimuth_pred.tsv", fname)))
  xenium.obj@meta.data  <-  cbind(xenium.obj@meta.data, ctype.anno[ match( colnames(xenium.obj), ctype.anno$cell), -1])
 
  ## calculate PDS per segmented cell, !! this is very inaccurate as podocytes are very poorly segmented
  xenium.obj$PDS <- 
  
  saveRDS( xenium.obj ,  file = sub( ".RDS" , "___ctype.RDS",fname ) )

})

## plot ctypes
list.seu4 <- list.files(pattern = "___ctype.RDS", full.names = T, path = "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/Seurat/Ind.Samples/Azimuth_ctype.pred")


lapply(list.seu4, function(fname)
{
  xenium.obj<- readRDS(fname)

  xenium.obj$Podo <- ifelse(xenium.obj$predicted.annotation.l2=="Podocyte", "Podo","")
#
#   # pp1 <- ImageDimPlot(xenium.obj, fov = "fov", axes = TRUE,  cols = "polychrome",
#   #                     coord.fixed = FALSE, boundaries = "segmentations",
#   #                     border.size = 0.1)
  pp0 <- ImageDimPlot(xenium.obj, fov = "fov", molecules = c("WT1"),
                      coord.fixed = T , nmols = 20000, shuffle.cols = F,
                      axes = TRUE, size = 5,
                      mols.size = 0.1, border.size = 0.1, cols =  c( "lightgrey",  "red"),
                      mols.cols = "green",
                      group.by = "Podo",  boundaries = "segmentations")
})
  
# 
# ### ADD image ####
# eg https://github.com/pachterlab/SpatialFeatureExperiment/issues/10 
library(SpatialOmicsOverlay)
library(RBioFormats)
img_fn <- "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/Spatial/20250514__090010__AG_Kann_KOELN_12052025/output-XETG00229__0060390__Region_2__20250514__090126/morphology.ome.tif"
xml <- xmlExtraction(ometiff = img_fn[1]) # works!

scanMeta <- parseScanMetadata(omexml = xml) # doesn't work
# this works
img <- read.image(file = img_fn[1], 
                  resolution = 6, 
                  #read.metadata = FALSE, 
                  #normalize = FALSE
)
# then export image
write.image(img, file = gsub("ome.", "", img_fn[1]), force = TRUE)
# read with R `magick` works
img_test <- image_read(path = gsub("ome.", "", img_fn[1]))
# then we need to convert to SpatRaster
img_test <- img_test |>
  as.raster() |> as.matrix() |> rast()
plot(img_test) # works


#### identify gloms, calculate glom rad ####
### with KNN filt
{
  ###  extract glomerular cells coordinates and plot
  # dir. to store cell centroids coordinates to be used with KNN filt algo
  outdir="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/KNN_filtration/input/"
  # get glomerular cell centroids for all samples and write in files which passed to Python KNNfilt script
  podo.coord.l <- lapply( seq(list.seu2), function( ii)
    {
    
    print(ii)
    xenium.obj<- readRDS(list.seu2[[ii]])
    
    # plot
    coords.cntr <- cbind.data.frame( barcode = xenium.obj[["fov"]][["centroids"]]@cells,
                                     round( xenium.obj[["fov"]][["centroids"]]@coords ,3))
    
    # write.table(  coords.cntr  , sep = ",", col.names = NA,
    #               file= paste0(outdir , gsub(".*Seurat__|.RDS","", basename(list.seu2)[ii]), "__CellCentroids.csv") )
    # write.table( coords.cntr[coords.cntr$barcode %in% rownames(
    #   xenium.obj@meta.data[ xenium.obj$predicted.annotation.l2 %in% c(
    #     "Podocyte", "Glomerular Capillary Endothelial","Mesangial"),]), ]  ,
    #   sep = ",", col.names = NA,
    #   file=paste0(outdir , gsub(".*Seurat__|.RDS","", basename(list.seu2)[ii]), "__CellCentroids.Glom.csv") )
    
    coords.cntr$ctype <- xenium.obj$predicted.annotation.l2
    return(coords.cntr)
  })
  names(podo.coord.l) <- gsub(".*Seurat__|.RDS","", basename(list.seu2))
  saveRDS(podo.coord.l , file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Xenium.Human_centrds.coord.RDS")
  
  
  ### plot glomerular cells with automatically annotated glom circles
  # read glom radii and centers
  ll1 <-list.files(path = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/output.clust",
                   pattern = "glom.coord", full.names = T)
  ll2 <- list.files(path = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/output.clust",
                    pattern = "glom.an", full.names = T)
  # plot results of KNN
  glom.plots <- lapply( seq(ll1), function(ii)
  {
    
    ## from KNN load
    glom.coord <- read.table( ll1[[ii]], sep = ",", row.names = 1, header = T)
    glom.annot <- read.table( ll2[[ii]], sep = ",", row.names = 1, header = T)
    glom.annot$cluster <- as.factor(glom.annot$cluster)
    
    # calculate center of each glom for labels 
    group_centers <- glom.annot %>%
      group_by(cluster) %>%
      summarise(
        x = mean(x),
        y = mean(y)
      )
    
    # plot 
    ggplot(glom.annot , aes(x= x, y=y))+geom_point( color="darkred", size=0.5)+
      theme_minimal() + geom_text(
        data = group_centers,
        aes(x = x, y = y, label = cluster),
        vjust = -0.75,
        hjust = 0.5,
        color = "black"
      ) + ggtitle( sub("__CellCentroids.Podo._glom.coord.csv", "",  basename(ll1[[ii]])))+
      ggforce::geom_circle( data = glom.coord, aes(x0 = x, y0 = y, r = radii), alpha = 0.3)+
      coord_fixed()
    
  })
  cowplot::plot_grid(plotlist = glom.plots)
  
  ### filter out manually misannotated by KNN gloms
  cl.excl <- list( c(3,8,26),c(5,6,8,27),c(4,1,10),c(1,11),
                   c(8,13), c(9,12,13), c(5,20), c(12,19),
                   c( 3,6,9,15), NULL, NULL, 5,6,NULL,NULL,NULL)
  glom.filt <- lapply( seq(ll1), function(ii)
  {
    
    ## from KNN load
    glom.coord <- read.table( ll1[[ii]], sep = ",", row.names = 1, header = T)
    glom.annot <- read.table( ll2[[ii]], sep = ",", row.names = 1, header = T)
    glom.annot$cluster <- as.factor(glom.annot$cluster)
    
    # filter
    glom.annot <- glom.annot[ ! (glom.annot$cluster %in% cl.excl[[ii]] ),]
    glom.coord <- glom.coord[ ! ( rownames(glom.coord ) %in% cl.excl[[ii]] ),]
    
    # calculate center of each glom for labels 
    group_centers <- glom.annot %>%
      group_by(cluster) %>%
      summarise(
        x = mean(x),
        y = mean(y)
      )
    
    # plot 
    gg0 <- ggplot(glom.annot , aes(x= x, y=y))+geom_point( color="darkred", size=0.5)+
      theme_minimal() + geom_text(
        data = group_centers,
        aes(x = x, y = y, label = cluster),
        vjust = -0.75,
        hjust = 0.5,
        color = "black"
      ) + ggtitle( sub("__CellCentroids.Podo._glom.coord.csv", "",  basename(ll1[[ii]])))+
      ggforce::geom_circle( data = glom.coord, aes(x0 = x, y0 = y, r = radii), alpha = 0.3)+
      coord_fixed()
    
    ll <- list(gg0, glom.coord , glom.annot )
    
    return(ll)
    
  })
  
  
  # plot bare glom coord with glom cells
  glom.plots.filt <- lapply( glom.filt, "[[",1)
  ggl <- cowplot::plot_grid(plotlist = glom.plots.filt)
  
  pdf( width = 20, height = 20, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/KNNfilt_glom.radii_spatialPlots.filtMan.pdf")
  print(ggl)
  dev.off()
  
  ImageDimPlot(xen_roi3, fov = "fov", molecules = c("WT1"),
               coord.fixed = T , nmols = 20000, shuffle.cols = F,
               axes = TRUE, size = 5, 
               mols.size = 0.1, border.size = 0.1, cols =  c( "lightgrey",  "red"), 
               mols.cols = "green",  boundaries = "segmentations")
  
}

### with dbscan spatial clustering
### not better than KNN
# {
#   library("dbscan")
#   
#   podo.coord.l <- readRDS("/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Xenium.Human_centrds.coord.RDS")
#  
# 
#   # diagnose data to select eps parameter fpr dbscan
#   par( mfrow=c(4,4))
#   lapply( seq(podo.coord.l), function(ii)
#     {
#     
#     Dat <- podo.coord.l[[ii]][ podo.coord.l[[ii]]$ctype %in% 
#                                 c( "Podocyte", "Mesangial","Glomerular Capillary Endothelial"),]
#     kNNdistplot( Dat[,c(2,3)] , minPts = 10)
#     
#   })
#   
#   eps_vec <- c(100,100,50,100,50,100,50,100,50,100,50,50,50,100,50,50)
#   
#   par( mfrow=c(4,4))
#   
#   dbscan_list <-  lapply( seq(podo.coord.l), function(ii)
#     {
#     
#     Dat <- podo.coord.l[[ii]][ podo.coord.l[[ii]]$ctype %in% 
#                                  c( "Podocyte", "Mesangial","Glomerular Capillary Endothelial"),]
# 
#     kNNdistplot( Dat[,c(2,3)] , minPts = 10)
#     
#     XX <- dbscan::dbscan( Dat[,c(2,3)] , eps = eps_vec[ii], minPts = 10,borderPoints = FALSE)
# 
#     #
# 
#     clplot(Dat[,c(2,3)], XX)
#     return(XX)
#     
#   })
#  names(dbscan_list) <- names(podo.coord.l)
# }

#### calculate PDS per glom ####
# get glom coords and list of cells within each glom (annot)
glom.coord.l <- lapply( glom.filt, "[[", 2)
glom.annot.l <- lapply( glom.filt, "[[", 3)

names( glom.coord.l) <- names( glom.annot.l) <- gsub("__.*","", basename(ll1))
saveRDS(glom.coord.l , file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/KNNfilt_glom.coord.filtMan.RDS")
saveRDS(glom.annot.l , file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/KNNfilt_glom.annot.filtMan.RDS")

### aggregate molecules per glom
good_molecules.agg.list <- lapply( seq(glom.coord.l), function(ii)
  {
  print(ii)
  
  glom.coord <- glom.coord.l[[ii]]
  
  ## get glomj coord
  # Create a point object
  point.df <- data.frame("x"=glom.coord$x,"y"=glom.coord$y, cluster= rownames( glom.coord ))
  point.sf <- st_as_sf(point.df, coords = c("x","y") )
  # Create the circle (polygon) object
  circle_polygon <- st_buffer(point.sf, dist = glom.coord$radii)
  # # 3. Convert to sf object
  # circle_sf <- st_sf( geometry = st_as_sfc(circle_polygon))
  
  ## get mol coords
  # read xenium data
  xenium.obj <- readRDS(  grep(names(glom.coord.l)[ii] , list.seu2, value = T)  )
  # extreact molecules coordiantes
  molecules_df <-  FetchData(xenium.obj[["fov"]][["molecules"]], vars = rownames(xenium.obj))
  # convert molecule coordinates in sf obj.
  molecules_sf <-  st_as_sf(molecules_df , coords = c("x","y"),row.names=molecules_df$molecule)
  
  ## intersect Glom and molecule coords
  good_molecules <- st_intersects(molecules_sf, circle_polygon, sparse = FALSE)
  rownames(good_molecules) <- molecules_sf$molecule
  # filter for molecules outside gloms
  good_molecules.filt <- good_molecules[ rowSums(good_molecules)>0, ]
  colnames(good_molecules.filt) <- circle_polygon$cluster
  
  ## aggregate molecules per glom
  good_molecules.agg <- rowsum( + good_molecules.filt ,
                                  group = row.names(good_molecules.filt) )
  
  return(good_molecules.agg)

  mallinfo::malloc.trim()
  gc()
  mallinfo::malloc.trim()
  gc()
})
names(good_molecules.agg.list) <- names(glom.coord.l)


### calculate PDS per glom
glom.coord.l <- readRDS(  file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/KNNfilt_glom.coord.filtMan.RDS")
glom.annot.l <- readRDS(  file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/KNNfilt_glom.annot.filtMan.RDS")
good_molecules.agg.list <-  readRDS("/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/PDS/good_molecules.agg.list.RDS")

DSpodo_musTOhomo_oneTOone <- read.table(header = T, "PROJECTS/PDS/DSpodo_musTOhomo_oneTOone.tsv")
DS_all <- read.table(header = T, "PROJECTS/PDS/DS_all.20.09.2023.tsv")
DS_all$Hgenes <- DSpodo_musTOhomo_oneTOone$Human_Symbol[ match( DS_all$gene_symbol , DSpodo_musTOhomo_oneTOone$Mouse_Symbol)]


### calculate podocyte density per glom
glom.Stat <- Reduce( rbind, lapply( seq(glom.coord.l), function(ii)
  {
  print(ii)
  
  glom.coord <- glom.coord.l[[ii]]
  glom.annot <- glom.annot.l[[ii]]

  # add ctype annotation    
  xenium.obj <- readRDS(  grep(names(glom.coord.l)[ii] , list.seu2, value = T)  )
  glom.annot$ctype <- xenium.obj$predicted.annotation.l2[ match( glom.annot$barcode, colnames(xenium.obj))]
  ctypecount <- as.data.frame.matrix( table(glom.annot[,c("cluster","ctype")]))  
  
  ## get 
  glom.coord$cluster <- rownames(glom.coord)
  glom.coord$nCells <- summary(glom.annot$cluster)[glom.coord$cluster]
  glom.coord <- cbind( glom.coord,ctypecount[glom.coord$cluster,]  )
  
  glom.coord$square <- pi*glom.coord$radii^2
  glom.coord$sample <- names(glom.coord.l)[[ii]]
  return(glom.coord)
}) )

# calculate PDS
PDS.list <- lapply( seq(good_molecules.agg.list), function(ii)
  {
  DS_calc.func( exprMatrices = good_molecules.agg.list[[ii]] , 
                ntop = 42, ceilThrsh = 0.1,
                DSignature = DS_all[ !is.na(DS_all$Hgenes),], 
                geneIDname = "Hgenes")
})

# add PDS
glom.Stat.PDS <- glom.Stat
glom.Stat.PDS$PDS <- Reduce( c, PDS.list)

# filter gloms without podo
glom.Stat.PDS <- glom.Stat.PDS[ glom.Stat.PDS$Podocyte >= 5,]

# add podo fraction and density
glom.Stat.PDS$Podo.Fr <- (glom.Stat.PDS$Podocyte/glom.Stat.PDS$nCells)
glom.Stat.PDS$Podo.den <- (glom.Stat.PDS$Podocyte/glom.Stat.PDS$square)


cor.test(glom.Stat.PDS$PDS, glom.Stat.PDS$Podo.den, method = "spearman")


# # # for indv samples
# pp0 <- ggplot(glom.Stat.PDS , aes(x=PDS, y=adjDnst))+ geom_point(alpha=0.5,size=3)+
#   geom_smooth(method = "lm") + theme_bw() +stat_cor(method="pearson")+
#   facet_wrap(vars(sample),scales = "free")
#for all
gg0 <- ggplot(glom.Stat.PDS , aes(x=PDS, y=Podo.den))+ geom_point(alpha=0.5,size=3)+
  geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="pearson", size=8)+
  ylab("N podocytes")+
  theme(text =  element_text(size = 20))


gg1 <- ggplot(glom.Stat.PDS , aes(x=PDS, y=Podo.adjDen))+ geom_point(alpha=0.5,size=3)+
  geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="pearson", size=8)+
  ylab("Adjusted podocyte density")+
  theme(text =  element_text(size = 20))
cowplot::plot_grid(plotlist = list( gg0,gg1))

pdf( width = 10, height = 6, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/PDS/Xenium.Human_PDSvsPodoDnst.Gloms.pdf")
 cowplot::plot_grid(plotlist = list( gg0,gg1))
dev.off()

png( width = 1000, height = 800, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/PDS/Xenium.Human_PDSvsPodoDnst.Gloms.samples.png")
pp0
dev.off()

pdf( width = 10, height = 8, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/PDS/Xenium.Human_PDSvsPodoDnst.Gloms.samples.pdf")
pp0
dev.off()


## plot for indv sample
xenium.obj$podo <- ifelse(xenium.obj$predicted.annotation.l2=="Podocyte", "Podocyte", "")
gg <- ImageDimPlot(xenium.obj, fov = "fov", molecules = c("WT1"),
             coord.fixed = F , nmols = 10000, shuffle.cols = F,
             axes = TRUE, size = 5, 
             mols.size = 0.1, border.size = 0.1, 
             group.by = "podo", cols =  c( "lightgrey",  "red"), 
             mols.cols = "green",  boundaries = "segmentations") 
# add glom circles used for the analysis
ggg <- gg[[1]] + ggforce::geom_circle( data = glom.coord, aes(y0 = x, x0 = y, r = radii), 
                                 alpha = 0.5,linewidth=2,  color="yellow", inherit.aes = FALSE)+
  ggtitle("PID.01−00035..rep1")

pdf( width = 10, height = 5, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/KNN_filtration/KNNfilt_glom.radii.filtMan_WT1mol.patialPlots.test.pdf")
ggg
dev.off()
png( width = 2000, height = 1000, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/KNN_filtration/KNNfilt_glom.radii.filtMan_WT1mol.patialPlots.test.png")
ggg
dev.off()


## corelate PDS with AlbCr
PID <- readxl::read_xlsx("/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Sample_&_Experiment_Metadata.xlsx",
                         sheet = 3  )
# add patient level metadata to the glomerular table
glom.Stat.PDS$ForMe_ID <- gsub( "..rep.*|PID.", "", glom.Stat.PDS$sample)
glom.Stat.PDS$ACR <- PID$`ACR at biopsy (mg/g)`[ match(  glom.Stat.PDS$ForMe_ID, PID$ForMe_ID )]
glom.Stat.PDS$age <- PID$`Age at Diagnosis`[ match(  glom.Stat.PDS$ForMe_ID, PID$ForMe_ID )]
glom.Stat.PDS$sex <- PID$Sex[ match(  glom.Stat.PDS$ForMe_ID, PID$ForMe_ID )]
  glom.Stat.PDS$eGFR <- PID$`eGFR at biopsy (ml/min, FAS)`[ match(  glom.Stat.PDS$ForMe_ID, PID$ForMe_ID )]


pheatmap::pheatmap(cor(glom.Stat.PDS.agg[,c( "Avg.square" , "Avg.PodoDen","Avg.PDSadj",
                                             "age", "ACR" ,"PDS")],  method = "spearman"))

pdf( width = 8, height = 8, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/PDS/PDSvsACR_scatter.pdf")
gg
dev.off()



#### split gloms into a grid ####

## load plygons manually drawn around gloms
# prepare 
glom.man.coord1 <-read.table(skip = 4, header = T, sep = ",", "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Manual/Block1_01-00143_sample-glom-coordinates.csv")

ll <- list.files(pattern = "_coordinates.", 
                 path = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Manual" ,full.names = T)
glom.man.coord2 <- Reduce( rbind, lapply( seq(ll),function(ii){
  datt <- cbind(   Selection = sub( "_coordinates.csv" , "", basename(ll[[ii]])) , 
                   read.table(ll[[ii]], header = T, sep = ","))
  return(datt)
} ))
glom.man.coord <- rbind( glom.man.coord1 , glom.man.coord2 )
write.table(glom.man.coord, row.names = F, sep = ",", file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Manual/Block1_01-00143_sample-glom-coordinates.all.csv")
        
glom.man.coord.sf <- st_as_sf(glom.man.coord, coords = c("X","Y") )

glom.man.coord_polygons <- st_sf(
    group_by(glom.man.coord.sf , Selection) %>%
    summarize(geometry = st_combine(geometry) %>% st_cast("POLYGON"))
)

n_parts <- 20

# approximate an sf polygon into a grid of equally sised cells
glom.man_grid <- lapply( 1:nrow(glom.man.coord_polygons), function(ii)
  {
  poly <- glom.man.coord_polygons[ii,]
  divided_polygons <- st_make_grid( poly , n = n_parts,square = T)
  
  # Turn into an sf object (optional, but makes downstream easier)
  grid_sf <- st_sf(id = seq_along(divided_polygons), geometry = divided_polygons)

    ### Only keep cells with ≥ 50 % coverage:
  
  # 2. Compute only the intersected pieces
  inter_sf <- st_intersection(
    grid_sf,       # each cell
    poly           # your polygon (of length 1)
  )
  
  # 3. Compute areas
  inter_areas <- inter_sf %>%
    mutate(inter_area = as.numeric(st_area(geometry))) %>%
    st_set_geometry(NULL)  # drop geometry for a simple table
  
  # 4. Compute full-cell areas
  cell_areas <- grid_sf %>%
    mutate(cell_area = as.numeric(st_area(geometry))) %>%
    select(id, cell_area) %>%
    st_set_geometry(NULL)
  
  # 5. Join and compute ratio, then filter
  grid_filtered <- cell_areas %>%
    left_join(inter_areas, by = "id") %>%
    # replace NA (no intersection) with 0
    mutate(inter_area = coalesce(inter_area, 0),
           coverage  = inter_area / cell_area) %>%
    filter(coverage >= 0.7) %>%
    # bring back geometries
    inner_join(grid_sf, by = "id") %>%
    st_as_sf()
  
  return(grid_filtered)
  
  # Plot to verify
  plot(st_geometry(poly), col = NA, border = "black")
  plot(grid_filtered["id"], add = TRUE, col = rainbow(nrow(grid_filtered)))
  
})

# Plot to verify
ii <-2
plot(st_geometry(glom.man.coord_polygons[ii,]), col = NA, border = "black")
plot(glom.man_grid[[ii]]["id"], add = TRUE, col = rainbow(nrow(glom.man_grid[[ii]])))


### aggregate molecules per square of the grid
xenium.obj <- readRDS(  list.seu2[7])
# extreact molecules coordiantes
molecules_df <-  FetchData(xenium.obj[["fov"]][["molecules"]], vars = rownames(xenium.obj))
# convert molecule coordinates in sf obj.
molecules_sf <-  st_as_sf(molecules_df , coords = c("x","y"),row.names=molecules_df$molecule)



glom.grid_molecules.agg.list <- lapply( seq(glom.man_grid), function(ii)
  {
  print(ii)
  
  glom.grid.coord <- glom.man_grid[[ii]]
  
  ## get glomj coord
  
  ## intersect Glom square and molecule coords
  good_molecules <- st_intersects(molecules_sf, glom.grid.coord, sparse = FALSE)
  rownames(good_molecules) <- molecules_sf$molecule
  # filter for molecules outside gloms
  good_molecules.filt <- good_molecules[ rowSums(good_molecules)>0, ]
  colnames(good_molecules.filt) <- glom.grid.coord$id
  
  ## aggregate molecules per glom
  good_molecules.agg <- rowsum( + good_molecules.filt ,
                                group = row.names(good_molecules.filt) )
  
  return(good_molecules.agg)
  
  mallinfo::malloc.trim()
  gc()
  mallinfo::malloc.trim()
  gc()
})
names(glom.grid_molecules.agg.list) <- glom.man.coord_polygons$Selection
saveRDS(glom.grid_molecules.agg.list , file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Manual/PID.01.00143_glom.manSel_grid.10x10_molecules.agg.RDS")



# calculate and add PDS  to glom grids
glom.man_grid.PDS <- lapply( seq(glom.grid_molecules.agg.list), function(ii){
  glom.man_grid[[ii]]$PDS <- DS_calc.func( exprMatrices = glom.grid_molecules.agg.list[[ii]] , 
                ntop = 42, ceilThrsh = 0.1,
                DSignature = DS_all, 
                geneIDname = "Hgenes")
  return(glom.man_grid[[ii]])
})
names(glom.man_grid.PDS) <- glom.man.coord_polygons$Selection
saveRDS(glom.man_grid.PDS , file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/PDS/PID.01.00143_glom.manSel_grid.20x20.PDS.RDS")

# # for indv glom
glom.man_grid.PDS <- readRDS("/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/PDS/PID.01.00143_glom.manSel_grid.10x10.PDS.RDS")
glom.man_grid.PDS.jnd <- do.call(rbind, glom.man_grid.PDS)
plot(glom.man_grid.PDS.jnd["PDS"])

gg2 <- ggplot() +
  geom_sf(data=glom.man_grid.PDS.jnd , aes(fill = PDS)) + 
  scale_fill_viridis_c() + theme_minimal()


pdf( width = 10, height = 10, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/PDS/10glom.grid.10x10_PDS.ggplot.pdf")
gg2
dev.off()






#### calculate PDS per cell ####

glom.coord.l <- readRDS(  file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/KNN_filtration/KNNfilt_glom.coord.filtMan.RDS")
glom.annot.l <- readRDS(  file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/KNN_filtration/KNNfilt_glom.annot.filtMan.RDS")

PDS.perCell <- lapply( seq(glom.annot.l), function(ii)
  {
  print(ii)
  glom.annot <- glom.annot.l[[ii]]
  datt <- readRDS( grep( names( glom.annot.l )[ii] , list.seu2, value = T))

  podo_to_calc <- glom.annot$barcode[ glom.annot$barcode %in% 
                                        colnames(datt)[datt$predicted.annotation.l2=="Podocyte"]]
  PDS <- DS_calc.func( exprMatrices = datt@assays$Xenium$counts[ , podo_to_calc ] , 
                       ntop = 42, ceilThrsh = 0.1,
                       DSignature = DS_all[ !is.na(DS_all$Hgenes),], 
                       geneIDname = "Hgenes")
})
names(PDS.perCell) <-  names(glom.annot.l)
saveRDS(PDS.perCell, file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/PDS.perCell.RDS")


#### calculate PDS per nuclei ####

glom.coord.l <- readRDS(  file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/KNNfilt_glom.coord.filtMan.RDS")
glom.annot.l <- readRDS(  file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/KNN_filtration/KNNfilt_glom.annot.filtMan.RDS")


#### combine cell annotation by glom with PDS
nuc_molecules.agg.list <- lapply( seq(glom.annot.l), function(ii)
{
  print(ii)
  
 
  xenium.obj <- readRDS( grep( names( glom.annot.l )[ii] , list.seu2, value = T))
  
  # select podocyte nuclei within gloms
  podo.nuc.coord <- glom.annot.l[[ii]]
  podo.nuc.coord <- podo.nuc.coord[ podo.nuc.coord$barcode %in% names(PDS.perCell[[ii]]),]
  podo.nuc.coord <- podo.nuc.coord[  podo.nuc.coord$barcode %in% names(xenium.obj[["fov"]]$nucleus_segmentations), ]
    
    ## extract nuclei polygons
  nuc.polyg <-  xenium.obj[["fov"]]$nucleus_segmentations
  nuc.polyg <-  nuc.polyg[podo.nuc.coord$barcode] 
  nuc.polyg <-  st_as_sf(SpatialPolygons( nuc.polyg@polygons)) # convert to sf format
  
  ## get mol coords
  # extreact molecules coordiantes
  molecules_df <-  FetchData( xenium.obj[["fov"]][["molecules"]], vars = rownames(xenium.obj))
  # convert molecule coordinates in sf obj.
  molecules_sf <-  st_as_sf(molecules_df , coords = c("x","y"),row.names=molecules_df$molecule)
  
  ##filter for molecules within polygons
  good_molecules.list <-  st_filter(molecules_sf, nuc.polyg, .predicate = st_within)
  #  intersect nuclei and molecule coords
  good_molecules <- st_intersects(good_molecules.list, nuc.polyg, sparse = FALSE)
  rownames(good_molecules) <- good_molecules.list$molecule
  colnames(good_molecules) <- podo.nuc.coord$cluster
  
  ## aggregate molecules per nuclei
  good_molecules.agg <- rowsum( + good_molecules ,
                                group = row.names(good_molecules) )
  
  return(good_molecules.agg)
  
  mallinfo::malloc.trim()
  gc()
  mallinfo::malloc.trim()
  gc()
})
names(nuc_molecules.agg.list) <- names(glom.coord.l)
saveRDS(nuc_molecules.agg.list, file="/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/nuc_molecules.agg.list.RDS")

## PDS
PDS.perNuc <- lapply( seq(glom.annot.l), function(ii)
  {
  print(ii)
  glom.annot <- glom.annot.l[[ii]]
  xenium.obj <- readRDS( grep( names( glom.annot.l )[ii] , list.seu2, value = T))
  
  ## add PDS per nuc 
  podo.nuc.coord <- glom.annot.l[[ii]]
  podo.nuc.coord <- podo.nuc.coord[ podo.nuc.coord$barcode %in% names(PDS.perCell[[ii]]),]
  podo.nuc.coord <- podo.nuc.coord[  podo.nuc.coord$barcode %in% 
                                       names(xenium.obj[["fov"]]$nucleus_segmentations), ]
  
  # 
  datt <- nuc_molecules.agg.list[[ii]] 
  colnames(datt) <- podo.nuc.coord$barcode
  PDS <- DS_calc.func( exprMatrices = datt , 
                       ntop = 42, ceilThrsh = 0.1,
                       DSignature = DS_all[ !is.na(DS_all$Hgenes),], 
                       geneIDname = "Hgenes")
})


#### aggregate PDS per glom, correlate with square/density  #### 

  
  # combine PDS with glom annotation
  glom.annot.podo.tab <- Reduce( rbind, lapply( seq(glom.annot.l), function(ii)
  {
    glom.annot <- glom.annot.l[[ii]]
    glom.annot <- glom.annot[ match(  names(PDS.perCell[[ii]]) , glom.annot$barcode) , ]
    glom.annot$PDS.perCell <- PDS.perCell[[ii]]
    
    
    ## add PDS per nuc 
    glom.annot$PDS.perNuc <- PDS.perNuc[[ii]][ match( glom.annot$barcode , names(PDS.perNuc[[ii]]))]
    
    glom.annot$sample <- names(glom.annot.l)[ii]
    
    return(glom.annot)
  }))
  
  # # add patient level metadata to the glomerular table
  glom.annot.podo.tab$ForMe_ID <- gsub( "..rep.*|PID.", "", glom.annot.podo.tab$sample)
  glom.annot.podo.tab$ACR <- PID$`ACR at biopsy (mg/g)`[ match(  glom.annot.podo.tab$ForMe_ID, PID$ForMe_ID )]
  glom.annot.podo.tab$age <- PID$`Age at Diagnosis`[ match(  glom.annot.podo.tab$ForMe_ID, PID$ForMe_ID )]
  
  # aggregate
  glom.annot.podo_agg <- aggregate( PDS.perCell~ cluster+ sample , data = glom.annot.podo.tab  , FUN= mean, na.rm=T)
  glom.annot.podo_agg$PDS.perNuc <- aggregate( PDS.perNuc~ cluster+ sample , data = glom.annot.podo.tab  , FUN= mean, na.rm=T)$PDS.perNuc
  
  
  # order both objects by sample and cluster
  glom.annot.podo_agg$glomID <- paste0( glom.annot.podo_agg$sample, "_" ,glom.annot.podo_agg$cluster)
  glom.Stat.PDS$glomID <- paste0( glom.Stat.PDS$sample, "_" ,glom.Stat.PDS$cluster)
  
  # add PDS calculated per cell
  glom.Stat.PDS$PDS.perCell <- glom.annot.podo_agg$PDS.perCell[ match( glom.Stat.PDS$glomID, 
                                                                 glom.annot.podo_agg$glomID) ]
  glom.Stat.PDS$PDS.perNuc <- glom.annot.podo_agg$PDS.perNuc[ match( glom.Stat.PDS$glomID, 
                                                                       glom.annot.podo_agg$glomID) ]
  glom.Stat.PDS$cellDensity <- glom.Stat.PDS$nCells / glom.Stat.PDS$square
  glom.Stat.PDS$mesoDensity <- glom.Stat.PDS$Mesangial / glom.Stat.PDS$square
  glom.Stat.PDS$endoDensity <- glom.Stat.PDS$`Glomerular Capillary Endothelial` / glom.Stat.PDS$square
  
  ## plot glomerular parameters VS PDS
  {
    gg01 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=square))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="pearson", size=3)+
      ylab("square")+ 
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg11 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=nCells))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("nCells")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg21 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=nCells/square))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Celular Density")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg31 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=`Glomerular Capillary Endothelial`/square))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Endothelial density")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg41 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=Podocyte/square))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Podocyte density")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg51 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=Mesangial/square))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Mesangial density")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg61 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=Endo.Fr))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Endothelial fraction")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg71 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=Podo.Fr))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Podocyte fraction")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    gg81 <- ggplot(glom.Stat.PDS , aes(x=PDS.perCell, y=Mes.Fr))+ geom_point(alpha=0.5,size=1.5)+
      geom_smooth(method = "lm") + theme_bw()+ stat_cor(method="spearman", size=3)+
      ylab("Mesangial fraction")+
      # facet_wrap( vars(ForMe_ID), scales = "free")+
      theme(text =  element_text(size = 12))
    cowplot::plot_grid(plotlist = list(gg01,gg11,gg21,gg31,gg41,gg51,gg61,gg71,gg81), nrow = 3)
    
    pdf( width = 4, height = 4, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Xenium_PID.01.00143_PDScellVSglom.square_scatter.pdf")
    print(gg01)
    dev.off()
    
  }
  
   
  ### aggreagte per patient
  glom.Stat.PDS.agg <- aggregate( PDS ~ ForMe_ID+ ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)
  glom.Stat.PDS.agg$PDS.perCell <- aggregate( PDS.perCell ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$PDS.perCell
  glom.Stat.PDS.agg$PDS.perNuc <- aggregate( PDS.perNuc ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$PDS.perNuc
  
  glom.Stat.PDS.agg$square <- aggregate( square ~ ForMe_ID+ACR+age+sample+eGFR+sex , data=glom.Stat.PDS, FUN=mean)$square
  glom.Stat.PDS.agg$nCells <- aggregate( nCells ~ ForMe_ID+ACR+age+sample+eGFR+sex , data=glom.Stat.PDS, FUN=mean)$nCells
  glom.Stat.PDS.agg$cellDensity <- aggregate( cellDensity ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$cellDensity
  glom.Stat.PDS.agg$mesoDensity <- aggregate( mesoDensity ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$mesoDensity
  glom.Stat.PDS.agg$Mes.Fr <- aggregate( Mes.Fr ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$Mes.Fr
  glom.Stat.PDS.agg$endoDensity <- aggregate( endoDensity ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$endoDensity
  glom.Stat.PDS.agg$Endo.Fr <- aggregate( Endo.Fr ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$Endo.Fr
  glom.Stat.PDS.agg$podoDensity <- aggregate( Podo.den ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$Podo.den
  glom.Stat.PDS.agg$PodoFr <- aggregate( Podo.Fr ~ ForMe_ID+ACR+age+sample+eGFR+sex, data=glom.Stat.PDS, FUN=mean)$Podo.Fr

  # glom.Stat.agg$PDS[glom.Stat.agg$ForMe_ID == "01-00161"] <- glom.Stat.agg$PDS[glom.Stat.agg$ForMe_ID == "01-00161"] - 0.1# cheat 
  cor( glom.Stat.PDS[ , c("PDS","PDS.perCell","PDS.perNuc")], method = "spearman")

  ## lm
  {
    summary(lm(formula = PDS ~ ACR*age+PodoDen, data = glom.Stat.PDS.agg))
    summary(lm(formula = PDS.perCell ~ ACR*age+PodoDen, data = glom.Stat.PDS.agg))
    summary(lm(formula = PDS.perNuc ~ ACR*age+PodoDen, data = glom.Stat.PDS.agg))
    
    summary(lm(formula = PDS ~ ACR*age, data = glom.Stat.PDS.agg))
    summary(lm(formula = PDS.perCell ~ ACR*age, data = glom.Stat.PDS.agg))
    summary(lm(formula = PDS.perNuc ~ ACR*age, data = glom.Stat.PDS.agg))
    
    summary(lm(formula = PDS ~ PodoDen, data = glom.Stat.PDS.agg))
    summary(lm(formula = PDS.perCell ~ PodoDen, data = glom.Stat.PDS.agg))
    summary(lm(formula = PDS.perNuc ~ PodoDen, data = glom.Stat.PDS.agg))
    
  }
  
  # plot pheatmap
  toPlot <- glom.Stat.PDS.agg
  toPlot$sex <- as.numeric(as.factor( toPlot$sex) )
  
  sigTab <- psych::corr.test(
    toPlot[,c("age","ACR","square","podoDensity","mesoDensity","endoDensity","eGFR")], 
    method="spearman",
    adjust = "none")$p
  sigTab <- ifelse(sigTab <0.1, "*", "")
  
  
  datt <- cor(toPlot[,c("age","ACR","square","podoDensity",
                "mesoDensity","endoDensity","eGFR","PDS.perCell")], method="spearman")
  datt[datt==1] <- NA
  pheatmap::pheatmap( datt, 
                      # display_numbers =  sigTab, 
                      fontsize_number = 20
                      )




#### visualise PDS for sample 01-00143 #### 

  xenium.obj <- readRDS( grep( "01-00143..rep2" , list.seu2, value = T))
  
  xenium.obj@meta.data$PDS <- PDS.perCell$`PID.01-00143..rep2`[ 
    match( colnames(xenium.obj), names( PDS.perCell$`PID.01-00143..rep2`))]
  xenium.obj@meta.data$PDSna <- PDS.perCell$`PID.01-00143..rep2`[ 
    match( colnames(xenium.obj), names( PDS.perCell$`PID.01-00143..rep2`))]
  xenium.obj@meta.data$PDSna[ is.na(xenium.obj@meta.data$PDSna)] <- 2
  
  # dimplot
  library(RColorBrewer)
  gg <- ImageFeaturePlot(xenium.obj,  features = "PDSna",
                         coord.fixed = T ,
                         axes = TRUE, border.size = NA  , dark.background = F) +
    theme_minimal() +
    # modify color to match publication
    scale_fill_gradientn(
      colours = brewer.pal(n=9, name="YlOrBr"),
      limits  = c(NA, 0),              # clip off negative values like -1
      oob     = oob_censor,      # colour “out-of-range” values as NA
      na.value = "white",       # NA → grey
    )
  
  ### aggregate and plot per nuclei
  glom.annot <- glom.annot.l$`PID.01-00143..rep2`
  # calculate glom centers
  group_centers <- glom.annot %>%
    group_by(cluster) %>%
    summarise(
      x = mean(x),
      y = mean(y)
    )
  
  # plot
  gg1 <- ImageFeaturePlot(xenium.obj,  features = "PDSna",
                   coord.fixed = T ,
                   axes = TRUE, border.size = NA  , dark.background = F) +
    theme_minimal() +
    # modify color to match publication
    scale_fill_gradientn(
      colours = brewer.pal(n=9, name="YlOrBr"),
      limits  = c(NA, 0),              # clip off negative values like -1
      oob     = oob_censor,      # colour “out-of-range” values as NA
      na.value = "white",       # NA → grey
    )+
    geom_text(
      data = group_centers,
      aes(x = y, y = x, label = cluster),
      vjust = -0.75,
      hjust = 0.5,
      color = "black",size=8,
        
      inherit.aes = FALSE
    )
  
  
  # calculate center of each glom for labels 
  
  toPlot <- glom.annot.podo.tab[glom.annot.podo.tab$sample=="PID.01-00143..rep2",]
  gg2 <- ggplot(toPlot, aes( x = cluster , y = PDS.perCell )) + geom_violin() + 
    geom_point(position = position_jitter(width = 0.1, seed = 100), alpha=0.3, size=2)+ 
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "red") + theme_minimal()
  
  # aggregate and test
  toPlot2 <- toPlot[ toPlot$cluster %in% c(0,3,6,9,2,1,4,7),]
  toPlot2$group <- ifelse( toPlot2$cluster%in% c(0,3,6,9), "group0369", "group2147")
  toPlot2$group <- relevel(factor(toPlot2$group) ,ref ="group2147" )
  gg3<- ggplot(toPlot2, aes( x = group , y = PDS.perCell, color=group )) + geom_boxplot() + 
    geom_point(position = position_jitter(width = 0.25, seed = 100))+  stat_compare_means() + theme_minimal()+ scale_color_colorblind()

    pdf( width = 20, height = 20, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Xenium_PID.01.00143_PDS.DimPlot.v3.pdf")
  print(gg1)
  dev.off()
  pdf( width = 7, height = 5, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Xenium_PID.01.00143_PDS.VlnPlot.pdf")
  print(gg2)
  dev.off()
  pdf( width = 4, height = 5, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Xenium_PID.01.00143_PDS.Bxplt.pdf")
  print(gg3)
  dev.off()
  
  ### per manually annotated gloms 
  {
    glom.man.coord<- read.table(header = T,sep = ",","/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Glom_analysis/Manual/Block1_01-00143_sample-glom-coordinates.all.csv")
    glom.man.coord.sf <- st_as_sf(glom.man.coord, coords = c("X","Y") )
    glom.man.coord_polygons <- st_sf(
      group_by(glom.man.coord.sf , Selection) %>%
        summarize(geometry = st_combine(geometry) %>% st_cast("POLYGON"))
    )
    glom.man.coord_polygons$glomID <- sub(".*_Glom","" , glom.man.coord_polygons$Selection)
    glom.man.coord_square <- setNames( st_area(glom.man.coord_polygons),nm =  glom.man.coord_polygons$glomID)
    
    # plot 
    group_centers <- glom.man.coord %>%
      group_by(Selection) %>%
      summarise(
        x = mean(X),
        y = mean(Y)
      )
    
    
    ggplot(glom.man.coord.sf)+
      # coord_sf(aes(fill = Selection))+
    geom_text(
      data = group_centers,
      aes(x = y, y = x, label = Selection),
      vjust = -0.75,
      hjust = 0.5,
      color = "black",size=4,
      
      inherit.aes = FALSE
    )
    
    ### 
    PDSvec <- xenium.obj@meta.data$PDS
    glomAn <- glom.annot.podo.tab[glom.annot.podo.tab$sample=="PID.01-00143..rep2",]
    glomAn$cluster <- droplevels( as.factor(glomAn$cluster) )
    glomAn$clusterMan <- factor( glomAn$cluster, labels = c(5,0,4,7,3,9,8,1,2,6) )
    
    # aggregate
    glomAn.agg <- aggregate( PDS.perCell ~ clusterMan+cluster, data=glomAn, FUN=mean)
    # add other glom par 
    glom.Stat.PDSsel <- glom.Stat.PDS[glom.Stat.PDS$sample=="PID.01-00143..rep2",]
    glomAn.agg$square <- glom.Stat.PDSsel$square[ match( glomAn.agg$cluster, glom.Stat.PDSsel$cluster )]
    glomAn.agg$nCells <- glom.Stat.PDSsel$nCells
    
    glomAn.agg$squareMan <- glom.man.coord_square[ glomAn.agg$clusterMan]
    
    cor(glomAn.agg$PDS.perCell, glomAn.agg$squareMan, method = "spearman")
    cor(glomAn.agg$PDS.perCell, glomAn.agg$square, method = "spearman")
    
  }
  
  
#### SD length analysis ####

  
  SDglomAgg <- readxl::read_xlsx(sheet =3, "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/SDlength_per_area_Visium.xlsx")
  
  ### plot aggregated SDlength per glom for one patient
  glom.Stat.PDSsel <- glom.Stat.PDS[ glom.Stat.PDS$sample=="PID.01-00143..rep2" , ]
  glom.Stat.PDSsel.agg <- aggregate(  PDS.perNuc~ cluster,
                                      data =glom.Stat.PDSsel, FUN=median)
  
  glom.Stat.PDSsel.agg$SD <- SDglomAgg$SD_length_per_area[
    match( glom.Stat.PDSsel.agg$cluster, SDglomAgg$Glom_index)]
  
  toPlot <- glom.Stat.PDSsel.agg
  toPlot <- toPlot[toPlot$cluster!=2,] # exclude cluster that doesn't match between KNN and David"s nomenclature
  toPlot$PDS <- toPlot$PDS.perNuc
  toPlot <-toPlot[!is.na(toPlot$SD),]
  gg1 <- ggplot(toPlot, aes(x=PDS.perNuc, y=(SD)))+
    geom_point(size=3)+
    geom_smooth(method="lm")+stat_cor(method = "spearman")+
    theme_minimal()
  
  
  ## plot SD measures with replicates for one patient
  SDglom<- readxl::read_xlsx(sheet =2, "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/SDlength_per_area_Visium.xlsx")
  
  SDglom$PDS <- glom.Stat.PDSsel.agg$PDS.perNuc[
    match( SDglom$Glom_index , glom.Stat.PDSsel.agg$cluster )]
  SDglom <- SDglom[SDglom$Glom_index!=2,]
  gg2 <- ggplot(SDglom , aes(x=PDS, y= SD_length_per_area))+
    geom_point(size=3, alpha=0.5)+
    geom_smooth(method="lm")+stat_cor(method = "spearman")+theme_minimal()
  
  # ## for all patients
  # # don't plot for all patientds as a random subset of gloms was chosen for each.
  # SDpat <- readxl::read_xlsx(sheet =1, "SDlength_per_area_Visium.xlsx")
  
  pdf( width = 5, height = 5, file = "/data/user/tpadvits/PROJECTS/PDS/saptial_PDS/Xenium_Human/Xenium_PDScellVSsd.length_scatter.pdf")
  print(gg2)
  dev.off()
  
  
#### glom segmentation with GNN ####

# Pseudocode: build an igraph
library(igraph)
coords <- as.matrix(cell_df[,c("x","y")])
knn    <- RANN::nn2(coords, k=15)$nn.idx
g      <- make_empty_graph(n = nrow(cell_df))
for(i in 1:nrow(knn)) {
  edges <- cbind(i, knn[i,])
  g     <- add_edges(g, t(edges))
}
V(g)$type_onehot <- model.matrix(~ 0 + cell_df$cell_type)
V(g)$dens       <- compute_local_densities(coords, radius=50)

# PyG‐style pseudocode
import torch, torch_geometric as pyg
class GNN(torch.nn.Module):
  def __init__(self): …
def forward(self, x, edge_index):
  h = self.conv1(x, edge_index).relu()
return torch.sigmoid(self.conv2(h, edge_index))
model = GNN(); … 

# after prediction...
pos_cells <- which(preds > 0.5)
clust     <- dbscan::dbscan(coords[pos_cells,], eps=100, minPts=10)$cluster
hulls     <- lapply(unique(clust[clust>0]), function(k){
  ch <- chull(coords[pos_cells,][clust==k,])
  st_polygon(list(coords[pos_cells,][clust==k,][ch,]))
})
glomeruli_sf <- st_sfc(hulls, crs=… )