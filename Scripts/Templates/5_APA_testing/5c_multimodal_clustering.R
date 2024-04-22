# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
# https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub#sec4.4.9

library(Seurat)
library(dplyr)
library(ggplot2)
library(goldmine)
library(cluster)
library(factoextra)

peak_count_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_peakcount/urothelial_counts.txt'
gene_count_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_genecount/urothelial_counts.txt'
seurat_obj_file <- '~/Desktop/Ureter10_scPASU_run/Seurat_objects/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds'
apa_res_dir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_APA_testing/differentiation_stage_cellranger_peakcount//'
peakref_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3h_fragmented_peaks_to_merge/Ureter10_urothelial_final_peak_universe_updated.txt'

### APA genes between basal, intermediate and umbrella
setwd(apa_res_dir)
files <- list.files()
for (file in grep('res.txt$', files, value = TRUE)) {
  assign(gsub('_res.txt','',file),read.delim(file))
}
union_apa_peaks <- subset(basal_v_intermediate, int_sig == TRUE)$peak
union_apa_peaks <- c(union_apa_peaks,subset(intermediate_v_umbrella, int_sig == TRUE)$peak)
union_apa_peaks <- c(union_apa_peaks,subset(umbrella_v_basal, int_sig == TRUE)$peak)
union_apa_peaks <- unique(union_apa_peaks)

# Functions

add_mt_perc<-function(s.obj,sample_id,verbose=FALSE)
{
  if(verbose)
  {cat(paste0('Calculating MT gene % for ',sample_id),sep='\n')}
  
  if(sum(PercentageFeatureSet(s.obj,pattern = "^MT-"))==0)
  { #In case of mouse dataset, '^mt-' is used
    if(sum(PercentageFeatureSet(s.obj,pattern = "^mt-"))==0)
    { if(verbose)
    {cat("No mitochondrial genes found ", '\n')}
    }else
    {
      s.obj[['perc.mt']]<-PercentageFeatureSet(s.obj,pattern = "^mt-")
    }
  }else
  {
    s.obj[['perc.mt']]<-PercentageFeatureSet(s.obj,pattern = "^MT-")
  }
  
  return(s.obj)
  
}

create_seurat_obj_from_counts_data<-function(counts_dat,sample_id,verbose=FALSE)
{
  # Create Seurat object
  if(verbose){
    cat(paste0('Creating Seurat Object for sample ',sample_id),sep='\n')
  }
  
  s.obj<-CreateSeuratObject(counts=counts_dat, project = sample_id)
  
  #Add percentage of mitochondrial genes
  s.obj<-add_mt_perc(s.obj,sample_id)
  
  return(s.obj)
  
}

cca_batch_correction<-function(s.obj.list,project.name,anchors=2000,int.genes='',verbose=FALSE)
{
  if(verbose){
    cat("Performing CCA-MNN Pipeline",sep='\n')
  }
  
  # Extract the normalization method [assuming all objects have same normalization method]
  norm<-s.obj.list[[1]]@commands$NormalizeData.RNA@params$normalization.method
  
  #Determine min number of neighbors
  mink<-min(200, min(sapply(seq_along(s.obj.list),function(x) ncol(s.obj.list[[x]]) ))  )
  
  # Select features
  features<-SelectIntegrationFeatures(obj.list,nfeatures=anchors)
  
  # Find Integration anchors
  sample.anchors<-FindIntegrationAnchors(s.obj.list,dims = 1:30,k.filter=mink,reduction='cca',anchor.features=features,normalization.method=norm)
  
  # Integrate with specified number of genes
  if(int.genes=='')
  {if(verbose)
  {cat("Integrating the sample.anchors only",sep='\n')}
    s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30,normalization.method=norm)
  }else if(int.genes=='all')
  {
    all.genes <- lapply(s.obj.list, row.names) %>% Reduce(intersect, .)
    
    if(verbose)
    {cat("CCA_MNN batch correction - integrating all genes", sep='\n')
      cat("Total genes being used for integration = ", length(all.genes),'\n')
    }
    s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30,features.to.integrate=all.genes,normalization.method=norm)
  }
  
  rm(sample.anchors)
  
  s.obj.integrated@project.name<-project.name
  
  return(s.obj.integrated)
}

#### Integrate expression assay
outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5c_clustering_cellranger_peakcount/2kfeatures_bothmodality/'
batch_correction <- 'cca-mnn'
hvg <- 2000
#hvg <- length(union_apa_peaks)
#modality <- 'apapeaks'
modality <- paste0(hvg,'hvpn')
#modality <- paste0(hvg,'hvg')
fprefix <- paste0('Ureter10_urothelial_',modality)
#file_prefix <- paste0(fprefix,'_apaintegrated')
file_prefix <- paste0(fprefix,'_exprintegrated')
batch_genes <- ''
pca_dimensions <- '1:50'

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

## Load data
sobj <- readRDS(seurat_obj_file)
sub.list <- list()
sub.count <- list()

# Expression
gene_count <- read.delim(gene_count_file)
colnames(gene_count) <- gsub('\\.','-',colnames(gene_count))
gene_count <- gene_count[,colnames(sobj@assays$RNA@counts)]
gene_count <- as.sparse(gene_count)
expr_assay <- CreateAssayObject(counts = gene_count)
sub <- sobj
sub[['Expr']] <- expr_assay
sub.list<-SplitObject(sub, split.by = "orig.ident")
sample.id<-names(sub.list)
sub.counts<-lapply(X=1:length(sub.list), function(x) {return(sub.list[[x]]@assays$Expr@counts)})
obj.list<-mapply(create_seurat_obj_from_counts_data,sub.counts,sample.id,verbose=TRUE)

# APA
peak_count <- read.delim(peak_count_file)
colnames(peak_count) <- gsub('\\.','-',colnames(peak_count))
peak_count <- peak_count[,colnames(sobj@assays$RNA@counts)]
# Exclude P0
peak_count <- peak_count[-grep(':P0$',row.names(peak_count)),]
# APA peaks only
# peak_count <- peak_count[union_apa_peaks,]
peak_count <- as.sparse(peak_count)
apa_assay <- CreateAssayObject(counts = peak_count)
sub <- sobj
sub[['APA']] <- apa_assay
sub.list<-SplitObject(sub, split.by = "orig.ident")
sample.id<-names(sub.list)
sub.counts<-lapply(X=1:length(sub.list), function(x) {return(sub.list[[x]]@assays$APA@counts)})
obj.list<-mapply(create_seurat_obj_from_counts_data,sub.counts,sample.id,verbose=TRUE)

rm(sub)
rm(sub.list)
rm(sub.counts)
rm(sample.id)

# Pre-process data
cat('Performing Log normalization \n')
n<-length(obj.list)
obj.list<-lapply(X=1:n,function(x){NormalizeData(obj.list[[x]],normalization.method = 'LogNormalize',scale.factor = 10^4, verbose =TRUE)})
cat('Preparing batch correction using',batch_correction, '\n')
if(batch_correction!='none'){
  if(batch_correction=='cca-mnn')
    {assay<-'integrated'
    reduction<-'pca'
    }else if(batch_correction=='harmony')
      {assay<-'RNA'
      reduction<-'harmony'
      obj.merged<-merge(x=obj.list[[1]],y=obj.list[2:length(obj.list)])
      DefaultAssay(obj.merged)<-'RNA'
      }else
        {cat('You did not enter correct batch correction method, setting it to cca-mnn')
          batch_correction<-'cca-mnn'
          assay<-'integrated'
          reduction<-'pca'
          }
}else{obj.list<-obj}

cat('Finding Variable genes per sample \n')
obj.list<-lapply(X=1:length(obj.list),function(x){FindVariableFeatures(obj.list[[x]],selection.method = 'vst', nfeatures=hvg, verbose = TRUE)})

if(length(obj.list)>1)
{
  # Integrate data with batch correction
  
  if(batch_correction=='cca-mnn'){
    cat('Performing cca-mnn batch integration \n')
    obj.integrated<-cca_batch_correction(obj.list,project.name=file_prefix, anchors=hvg, int.genes=batch_genes, verbose=TRUE)
    
  }else if(batch_correction=='harmony')
  {cat('Performing harmony batch integration \n')
    # Variable features for merged obj
    obj.merged<-FindVariableFeatures(obj.merged,selection.method = 'vst', nfeatures=hvg, verbose = TRUE)
    
    # Scaling
    all.features<-rownames(obj.merged)
    obj.merged<-ScaleData(obj.merged,features=all.features, verbose=TRUE,assay = 'RNA')
    
    # PCA
    obj.merged<-RunPCA(obj.merged, npcs=50, verbose=TRUE,assay = 'RNA')
    
    #Run Harmony
    obj.integrated<-RunHarmony(obj.merged,
                               group.by.vars = "orig.ident",
                               plot_convergence = FALSE,
                               assay.use = 'RNA')
  }
}else
{# This allows skipping batch correction in case user wants to run the rest of the pipeline on an already batch corrected object (most likely batch corrected by another method)
  obj.integrated<-obj.list}

rm(obj.list)

## Run PCA ## (skip if ran harmony or if reduction already exists)

if(is.null(obj.integrated@reductions$pca)==TRUE)
{
  if(is.null(obj.integrated@reductions$harmony)==TRUE)
  {
    # Perform PCA since no reduction assay found
    DefaultAssay(obj.integrated)<-'integrated'
    cat('Running PCA on integrated data')
    all.features<-rownames(obj.integrated)
    obj.integrated<-ScaleData(obj.integrated,features=all.features, verbose=TRUE)
    
    # Run PCA
    obj.integrated<-RunPCA(obj.integrated, npcs=50,ndims.print = 1:15, verbose=TRUE)
    
  }
}

saveRDS(obj.integrated,file=paste0(outdir,file_prefix,"only.rds"))

## Elbow plot
options(bitmapType='cairo')
png(file=paste0(outdir,file_prefix,"_PCA_elbow_plots.png"),width = 11,height = 8.5, units='in', res=300)
ElbowPlot(obj.integrated,ndims=50,reduction=reduction)+ggtitle(paste0(file_prefix,'_ElbowPlot'))
dev.off()

# Calculate silhouette score
seurat_obj_file <- '~/Desktop/Ureter10_scPASU_run/Seurat_objects/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds'
outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5c_clustering_cellranger_peakcount/2kfeatures_bothmodality/'
gene_count_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_genecount/urothelial_counts.txt'
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}
setwd(outdir)
apa <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5c_clustering_cellranger_peakcount/2kfeatures_bothmodality/Ureter10_urothelial_2000hvpn_apaintegratedonly.rds'
expr <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5c_clustering_cellranger_peakcount/2kfeatures_bothmodality/Ureter10_urothelial_2000hvg_exprintegratedonly.rds'
sobj <- readRDS(seurat_obj_file)
apa_obj.integrated <- readRDS(apa)
expr_obj.integrated <- readRDS(expr)

gene_count <- read.delim(gene_count_file)
colnames(gene_count) <- gsub('\\.','-',colnames(gene_count))
gene_count <- gene_count[,colnames(sobj@assays$RNA@counts)]
gene_count <- as.sparse(gene_count)
expr_assay <- CreateAssayObject(counts = gene_count)

# Calculate average silhouette score of old cluster identities on old embeddings
reduction <- 'umap'
dims <- 1:2
file_prefix <- 'u10_uro_old'
file_prefix <- paste0(file_prefix,reduction,'_','oldlabels')
clusters <- sobj@meta.data[,'integrated_snn_res.0.2']
dist.matrix <- dist(x=Embeddings(object = sobj[[reduction]])[,dims])
stopifnot(identical(Embeddings(object = sobj[[reduction]])[,dims] %>% row.names(), row.names(sobj@meta.data)))
sil<-silhouette(x=as.numeric(as.factor(clusters)),dist=dist.matrix)
fviz_silhouette(sil)
ggsave(filename=paste0(file_prefix,"_silhouette.png"), width=33,height=10)
saveRDS(sil,paste0(file_prefix,"_silhouette.rds"))

z <- CreateSeuratObject(counts = expr_assay)
expr_integrated_assay <- CreateAssayObject(data = expr_obj.integrated@assays$integrated@data, key = 'exprintegrated_')
apa_integrated_assay <- CreateAssayObject(data = apa_obj.integrated@assays$integrated@data, key = 'apaintegrated_')

z@assays[['exprintegrated']] <- expr_integrated_assay
z@assays[['apaintegrated']] <- apa_integrated_assay
z@reductions$exprintegrated_pca <- expr_obj.integrated@reductions$pca
z@reductions$apaintegrated_pca <- apa_obj.integrated@reductions$pca

expr_numpc <- 1:25
apa_numpc <- 1:20
resolutions <- seq(0.1,0.5,0.1)

z <- FindMultiModalNeighbors(
  z, reduction.list = list("exprintegrated_pca", "apaintegrated_pca"), 
  dims.list = list(expr_numpc, apa_numpc), modality.weight.name = c("RNA.weight","APA.weight")
)

for (res in resolutions) {
  z <- FindClusters(z, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)
}

z <- RunUMAP(z, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
z@meta.data$integrated_snn_res.0.2 <- sobj@meta.data$integrated_snn_res.0.2[match(row.names(z@meta.data),
                                                                                  row.names(sobj@meta.data))]

for (i in 1:length(resolutions)) {
  assign(paste0('p',i),DimPlot(z, reduction = 'wnn.umap', group.by = paste0('wsnn_res.',resolutions[i]), label = TRUE, repel = TRUE, label.size = 2.5))
}
p0 <- DimPlot(z, reduction = 'wnn.umap', group.by = 'integrated_snn_res.0.2', label = TRUE, repel = TRUE, label.size = 2.5)
p0 + p1 + p2 + p3 + p4 + p5 
ggsave(paste0(outdir,'u10_uro_2khvpn_integrated_apa_PC',max(apa_numpc),'_2khvg_integrated_expr_PC',max(expr_numpc),'_clustered.png'), 
       width = 20, height = 14, units = 'in', bg = 'white')

# Calculate average silhouette score of old cluster identities on new embeddings
reduction <- 'wnn.umap'
dims <- 1:2
file_prefix <- paste0('u10_uro_PC',max(apa_numpc),'apaintPC',max(expr_numpc),'exprint_')
file_prefix <- paste0(file_prefix,reduction,'_','oldlabels')
clusters <- z@meta.data[,'integrated_snn_res.0.2']
dist.matrix <- dist(x=Embeddings(object = z[[reduction]])[,dims])
stopifnot(identical(Embeddings(object = z[[reduction]])[,dims] %>% row.names(), row.names(z@meta.data)))
sil<-silhouette(x=as.numeric(as.factor(clusters)),dist=dist.matrix)
fviz_silhouette(sil)
ggsave(filename=paste0(file_prefix,"_silhouette.png"), width=33,height=10)
saveRDS(sil,paste0(file_prefix,"_silhouette.rds"))

# Calculate average silhouette score of new cluster identities on new embeddings
for (res in resolutions) {
  file_prefix <- paste0('u10_uro_PC',max(apa_numpc),'apaintPC',max(expr_numpc),'exprint_')
  file_prefix <- paste0(file_prefix,reduction,'_wsnn_res.',res)
  clusters <- z@meta.data[,paste0('wsnn_res.',res)]
  dist.matrix <- dist(x=Embeddings(object = z[[reduction]])[,dims])
  stopifnot(identical(Embeddings(object = z[[reduction]])[,dims] %>% row.names(), row.names(z@meta.data)))
  sil<-silhouette(x=as.numeric(as.factor(clusters)),dist=dist.matrix)
  f <- fviz_silhouette(sil)
  ggsave(filename=paste0(file_prefix,"_silhouette.png"), width=33,height=10, f)
  saveRDS(sil,paste0(file_prefix,"_silhouette.rds"))
}

saveRDS(z,paste0('u10_uro_2khvpn_integrated_apa_PC',max(apa_numpc),'_2khvg_integrated_expr_PC',max(expr_numpc),'_clustered.rds'))
