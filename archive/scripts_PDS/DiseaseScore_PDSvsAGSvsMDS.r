
sne_nn <- readRDS("/cellfile/datapublic/flopes/All_Data_for_Tim.rds")
saveRDS(sne_nn , file = "/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/publication_CellGenomics2026_TimPadvitski/OldYoungArDr_Seurat.SCT.rds")


sne_raw <- readRDS("/cellfile/datapublic/flopes/All_Data_unormalized_for_Tim.rds")

cct <- sne_nn$seurat_clusters[ match( colnames(sne_raw), colnames(sne_nn))]
names(cct) <- colnames(sne_raw)
sne_raw$ctype <- cct
  
frncs_podo <- subset(sne_nn, subset =seurat_clusters=="Podo" )

frncs.raws_podo <- subset(sne_i, subset =ctype=="Podo" )

saveRDS(frncs_podo , file = "/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/publication_CellGenomics2026_TimPadvitski/OldYoungArDr_Seurat_podo.raw.rds")
VlnPlot(frncs_podo, features = "nCount_RNA", group.by = "age")
VlnPlot(frncs_podo, features = "nCount_RNA", group.by = "diet")
VlnPlot(frncs_podo, features = "nCount_RNA", group.by = "condition")
VlnPlot(frncs_podo, features = "nCount_RNA", group.by = "batch")
VlnPlot(frncs_podo, features = "Wt1", group.by = "condition")

###
sne_raw <-  FindVariableFeatures(sne_raw, selection.method = "vst", nfeatures = 2000)
sne_raw <- ScaleData(sne_raw)
sne_raw <- RunPCA(sne_raw, features = VariableFeatures(object = sne_raw))
ElbowPlot(sne_raw, ndims = 50)

sne_raw <- FindNeighbors(sne_raw, dims = 1:50)
sne_raw <- FindClusters(sne_raw, resolution = 0.1)
sne_raw <- RunUMAP(sne_raw, dims = 1:50)

FeaturePlot( sne_raw, label = T, features = "Wt1", max.cutoff = 10)
DimPlot( sne_raw, label = T)
frncs.raws_podoAll <- subset(sne_raw, subset =seurat_clusters==17 )
saveRDS(frncs.raws_podoAll , file = "/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/publication_CellGenomics2026_TimPadvitski/OldYoungArDr_Seurat_podoAll.raw.rds")

### read podocyte / glomerular circadfian genes
glom.circ <- readxl::read_xlsx( "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/RNA_seq/Public_data/Mouse_data/glomerular_circadian_PRJEB67585/Glomerular circadian transcriptomics meta2d BHQ 0.1.xlsx")
glom.circ <- glom.circ[1:375,]
podo.circ <- readxl::read_xlsx( "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/RNA_seq/Public_data/Mouse_data/glomerular_circadian_PRJEB67585/Podocyte circadian transcriptomics meta2d BHQ 0.1 .xlsx")
