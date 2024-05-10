#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(Seurat)
library(stringr)

# Local run
counts_dir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/4b_cellranger_genecount/'
counts_dir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/4b_cellranger_peakcount/'
fprefix <- 'urothelial'
seurat_obj_path <- '~/Desktop/Ureter10_scPASU_run/Seurat_objects/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds'
outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_genecount/'
outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_peakcount/'
samples <- c('U1','U2','U3','U4','U5','U6','U7','U8','U9','U10')

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

f <- paste0(counts_dir,samples,'/outs/raw_feature_bc_matrix/')

peak_counts<-lapply(f,function(x) {
  mtx <- Read10X(data.dir = x)
  mtx <- as.data.frame(mtx)
  return(mtx)
}
)

names(peak_counts) <- samples

# Modify cell barcodes in accordance with previous scRNA-seq analyses
if (is.na(seurat_obj_path)==FALSE) {
  cat('Modifying cell barcodes to match Seurat when integrated \n')
  z <- readRDS(seurat_obj_path)
  meta <- z@meta.data
  write.table(meta,paste0(outdir,fprefix,'_meta.txt'),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  meta$sample_idx <-as.numeric(sub(".*_([0-9]+)$", "\\1", row.names(meta)))
  sample_idx <- unique(meta$sample_idx)
  names(sample_idx) <- meta[match(sample_idx,meta$sample_idx),'sample']
  
  for (i in 1:length(peak_counts)) {
    for (s in names(sample_idx)) {
      if (names(peak_counts)[i]==s) {
        colnames(peak_counts[[i]]) <- sub(".*_(.*)", "\\1",colnames(peak_counts[[i]]))
        colnames(peak_counts[[i]]) <- paste0(colnames(peak_counts[[i]]),'_',sample_idx[names(sample_idx)==s])
      }
    }
  }
}

# Merge samples
names(peak_counts) <- NULL
merged_counts<-do.call(cbind,peak_counts)

stopifnot(identical(sort(colnames(merged_counts)),sort(row.names(z@meta.data))))

write.table(merged_counts,paste0(outdir,fprefix,'_counts.txt'),row.names = TRUE, col.names = TRUE, quote = FALSE, sep='\t')
