setwd("~/changrila2_mount/ms_epic/")

library(Seurat)

# create individual seurat object for each sample
all_sample_dirs <- list.files(path=".", pattern="-WTA$", include.dirs = T)

seurat_obj_list <- list()
for (dir in all_sample_dirs) {
  dirname <- unlist(strsplit(dir, '-'))[1]
  print(dirname)
  data <- Read10X(data.dir = paste0(dir,"/outs/filtered_feature_bc_matrix/"))
  
  # Initialize the Seurat object with the raw (non-normalized data)
  seurat_obj <- CreateSeuratObject(counts = data, project = dirname, min.cells = 3, min.features = 200)
  # seurat_obj$sample <- dirname
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj_list[[dirname]] <- seurat_obj
}

# merge datasets ----

immune.anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, dims = 1:20)
merged.integrated <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

save(immune.anchors, file="r_objects/immune.anchors.RData")
save(seurat_obj_list, file="r_objects/seurat_obj_list.RData")


# merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[2:length(seurat_obj_list)],
#                 add.cell.ids = names(seurat_obj_list), project = "MS_scRNA", do.normalize=FALSE)
# merged
# merged.list <- SplitObject(object = merged, split.by = "datasets")
# 
# reference.list <- merged.list[c("dataset1", "dataset2")] # change
# merged.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
# merged.integrated <- IntegrateData(anchorset = merged.anchors, dims = 1:50)

DefaultAssay(object = merged.integrated) <- "integrated"
merged.integrated <- ScaleData(object = merged.integrated, verbose = FALSE)
merged.integrated <- RunPCA(object = merged.integrated, npcs = 50, verbose = FALSE)

# merged.integrated <- RunTSNE(object = merged.integrated, reduction = "pca", dims = 1:50)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
merged.integrated <- FindNeighbors(object = merged.integrated)
merged.integrated <- FindClusters(object = merged.integrated, resolution = 0.2)

p1 <- DimPlot(object = merged.integrated, reduction = "tsne", group.by = "datasets")
p2 <- DimPlot(object = merged.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.2", 
              label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2)


# method 1: integrate
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# method 2: just merge (probably wrong)
ms.rna <- merge(seurat_obj_list[[1]], y = seurat_obj_list[2:length(seurat_obj_list)],
                add.cell.ids = names(seurat_obj_list), project = "MS_scRNA")
ms.rna

# remove mitochondrial genes
bcc.seurat[["percent.mt"]] <- PercentageFeatureSet(bcc.seurat, pattern = "^MT-")
VlnPlot(bcc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by="patient_id", ncol = 3)

ggplot(pbmc@meta.data, aes(nUMI, percent.mito)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.09, linetype = "dashed", colour = "red")