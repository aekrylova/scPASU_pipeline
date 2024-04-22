#https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# Load libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(apeglm)
library(DESeq2)
library(tibble)

outdir <- "~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_DEG_testing/differentiation_stage_cellranger_genecount/"
script_dir <- '~/Desktop/Ureter10_scPASU_run/Scripts/5_APA_testing/'
meta <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_genecount/urothelial_meta.txt'
indir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_genecount/'
counts_file <- 'urothelial_counts.txt'

source(paste0(script_dir,'scPASU_functions_for_differential_testing.R'))
if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

meta <- read.delim(meta)
meta$seurat_clusters <- paste0('cluster',meta$seurat_clusters)
counts <- read.delim(paste0(indir,counts_file))
colnames(counts) <- gsub('\\.','-',colnames(counts))

counts_by_cellstate <- stratify_matrix(merged_counts = counts, meta = meta, vars = c('seurat_clusters'),
                                     cutoff_pct = 1, min_cell_per_group = 10)

##### basal v. intermediate
inputs <- create_test_inputs(test = 'DEG', all_groups = counts_by_cellstate, ident1=c('cluster2','cluster3','cluster6'),
                             ident2=c('cluster0','cluster1','cluster5'), APA.feature.filter = FALSE, 
                             replicate = 'random', nrep = 3 ,p = 0.7)

colnames(inputs) <- c(paste0('basal_rep',1:3),paste0('intermediate_rep',1:3))

sig_genes <- DEG_DESeq2_pseudobulk(inputs = inputs, comp = c('intermediate','basal'), test.used = 'LRT',
                                   padj_thres=0.01,l2fc_frac_thres=log2(1.5),outdir)

#> The fit type with the least median absolute residual is glmGamPoi 
#> Differential testing using glmGamPoi as fit type and LRT as test 

##### intermediate v. umbrella
inputs <- create_test_inputs(test = 'DEG', all_groups = counts_by_cellstate, ident1=c('cluster0','cluster1','cluster5'),
                             ident2=c('cluster4','cluster7'), APA.feature.filter = FALSE, 
                             replicate = 'random', nrep = 3 ,p = 0.7)

colnames(inputs) <- c(paste0('intermediate_rep',1:3),paste0('umbrella_rep',1:3))

sig_genes <- DEG_DESeq2_pseudobulk(inputs = inputs, comp = c('umbrella','intermediate'), test.used = 'LRT',
                                   padj_thres=0.01,l2fc_frac_thres=log2(1.5),outdir)

#> The fit type with the least median absolute residual is glmGamPoi 
#> Differential testing using glmGamPoi as fit type and LRT as test 

##### umbrella v. basal
inputs <- create_test_inputs(test = 'DEG', all_groups = counts_by_cellstate, ident1=c('cluster4','cluster7'),
                             ident2=c('cluster2','cluster3','cluster6'), APA.feature.filter = FALSE, 
                             replicate = 'random', nrep = 3 ,p = 0.7)

colnames(inputs) <- c(paste0('umbrella_rep',1:3),paste0('basal_rep',1:3))

sig_genes <- DEG_DESeq2_pseudobulk(inputs = inputs, comp = c('basal','umbrella'), test.used = 'LRT',
                                   padj_thres=0.01,l2fc_frac_thres=log2(1.5),outdir)

#> The fit type with the least median absolute residual is glmGamPoi 
#> Differential testing using glmGamPoi as fit type and LRT as test 
