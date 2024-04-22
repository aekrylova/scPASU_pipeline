#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(reshape)
library(stringr)
library(data.table)
library(edgeR)
library(DEXSeq)
library(Seurat)

indir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5a_merged_cellranger_peakcount/'
script_dir <- '~/Desktop/Ureter10_scPASU_run/Scripts/5_APA_testing/'
counts_file <- 'urothelial_counts.txt'
meta_file <- 'urothelial_meta.txt'
peak_ref_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3h_fragmented_peaks_to_merge/Ureter10_urothelial_final_peak_universe_updated.txt'
outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_APA_testing/differentiation_stage_cellranger_peakcount2/'

source(paste0(script_dir,'scPASU_functions_for_differential_testing.R'))

# Read merged counts
merged_counts<-read.delim(paste0(indir,counts_file))
colnames(merged_counts) <- gsub('\\.','-',colnames(merged_counts))

# Read meta data [save meta data from Seurat object - ensure the set of cells match that of the count matrix]
meta<-read.table(paste0(indir,meta_file))
meta$seurat_clusters <- paste0('cluster',meta$seurat_clusters)

# Peak ref is jtu$join and final_annotation is the peak column 
peak_ref<-read.table(paste0(peak_ref_file),header=TRUE,sep='\t')

# Replace peak name column with final annotation as this is more meaningful name 
cat('Use final annotations as the peak names \n')
col<-which(colnames(peak_ref)=='final_annotation')
peak_ref$peak<-peak_ref[,col]
peak_ref<-peak_ref[,-col]

## Remove P0 peaks ##
cat('Remove P0 peaks \n')
plist<-strsplit(peak_ref$peak,split=':') %>% sapply(.,'[[',3)
rem<-which(plist=='P0')
peak_ref<-peak_ref[-rem,]

### DEXSeq test ###

apa <- stratify_matrix(merged_counts=merged_counts, meta=meta, vars=c('seurat_clusters'),cutoff_pct = 0, min_cell_per_group = 0)
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)}

###### basal v. intermediate ######

inputs <- create_test_inputs(test='APA',apa,ident1=c('cluster2','cluster3','cluster6'),ident2=c('cluster0','cluster1','cluster5'),min_cell_expr_pct=10,
                             expr.thres=1,pseudo.count=1,replicate='random',nrep=3,p=0.7)

colnames(inputs) <- c(paste0('basal_rep',1:3),paste0('intermediate_rep',1:3))

### DEXSeq test ###
DEXseq_res <- APA_DEXseq_test(ident1='basal',ident2='intermediate',inputs,min_peak=2,ncpu=4,dispersion.plot.save=TRUE,peak_ref=peak_ref,outdir=outdir)

### T Test ###
DEXseq_res <- APA_TTest(DEXseq_res,ident1='basal',ident2='intermediate')

## Check significance ##
DEXseq_res <- callsig(DEXseq_res,ident1='basal',ident2='intermediate',delta_frac_thres=0.1, padj_thres=0.01,
                      l2fc_frac_thres=log2(1.5),mean_frac_thres=0.05)

saveRDS(DEXseq_res,paste0(outdir,'basal_v_intermediate_APAtest.res.rds'))

###### intermediate v. umbrella  ######

inputs <- create_test_inputs(test='APA',apa,ident1=c('cluster0','cluster1','cluster5'),ident2=c('cluster4','cluster7'),min_cell_expr_pct=10,
                             expr.thres=1,pseudo.count=1,replicate='random',nrep=3,p=0.7)

colnames(inputs) <- c(paste0('intermediate_rep',1:3),paste0('umbrella_rep',1:3))

### DEXSeq test ###
DEXseq_res <- APA_DEXseq_test(ident1='intermediate',ident2='umbrella',inputs,min_peak=2,ncpu=4,dispersion.plot.save=TRUE,peak_ref=peak_ref,outdir=outdir)

### T Test ###
DEXseq_res <- APA_TTest(DEXseq_res,ident1='intermediate',ident2='umbrella')

## Check significance ##
DEXseq_res <- callsig(DEXseq_res,ident1='intermediate',ident2='umbrella',delta_frac_thres=0.1, padj_thres=0.01,
                      l2fc_frac_thres=log2(1.5),mean_frac_thres=0.05)

saveRDS(DEXseq_res,paste0(outdir,'intermediate_v_umbrella_APAtest.res.rds'))

###### umbrella v. basal ######

inputs <- create_test_inputs(test='APA',apa,ident1=c('cluster4','cluster7'),ident2=c('cluster2','cluster3','cluster6'),min_cell_expr_pct=10,
                             expr.thres=1,pseudo.count=1,replicate='random',nrep=3,p=0.7)

colnames(inputs) <- c(paste0('umbrella_rep',1:3),paste0('basal_rep',1:3))

### DEXSeq test ###
DEXseq_res <- APA_DEXseq_test(ident1='umbrella',ident2='basal',inputs,min_peak=2,ncpu=4,dispersion.plot.save=TRUE,peak_ref=peak_ref,outdir=outdir)

### T Test ###
DEXseq_res <- APA_TTest(DEXseq_res,ident1='umbrella',ident2='basal')

## Check significance ##
DEXseq_res <- callsig(DEXseq_res,ident1='umbrella',ident2='basal',delta_frac_thres=0.1, padj_thres=0.01,
                      l2fc_frac_thres=log2(1.5),mean_frac_thres=0.05)

saveRDS(DEXseq_res,paste0(outdir,'umbrella_v_basal_APAtest.res.rds'))

### Generate raw counts files and UCSC bedGraph tracks for significant APA genes
## Pairwise
pairwise_res <- list()
i = 1
for (comp in c('basal_v_intermediate','intermediate_v_umbrella','umbrella_v_basal')) {
  DEXseq_res <- readRDS(paste0(outdir,comp,'_APAtest.res.rds'))
  pairwise_res[[i]] <- DEXseq_res$res
  i = i + 1
}

names(pairwise_res) <- c('basal_v_intermediate','intermediate_v_umbrella','umbrella_v_basal')

bedtracks_list <- list()
i = 1
j = 1
for (comp in c('basal_v_intermediate','intermediate_v_umbrella','umbrella_v_basal')) {
  ident1 <- sapply(str_split(comp,pattern='_v_',n = 2),`[`,1)
  ident2 <- sapply(str_split(comp,pattern='_v_',n = 2),`[`,2)
  bedtracks_list[[i]] <- as.data.frame(pairwise_res[[j]])[,c('chr','start','end',paste0('mean_',ident1,'_frac'),'strand','int_sig')]
  names(bedtracks_list)[i] <- paste0(comp,'-',ident1)
  i = i + 1
  bedtracks_list[[i]] <- as.data.frame(pairwise_res[[j]])[,c('chr','start','end',paste0('mean_',ident2,'_frac'),'strand','int_sig')]
  names(bedtracks_list)[i] <- paste0(comp,'-',ident2)
  i = i + 1
  j = j + 1
}

# Create bedGraph tracks
for (i in 1:length(bedtracks_list)) {
  bedGraph <- bedtracks_list[[i]] 
  ident <- names(bedtracks_list[i])
  plus <- bedGraph[which(bedGraph$strand == '+'),]
  plus <- plus[,!(names(plus) %in% c('strand','int_sig'))]
  minus <- bedGraph[which(bedGraph$strand == '-'),]
  minus <- minus[,!(names(minus) %in% c('strand','int_sig'))]
  write.table(plus,paste0(outdir,ident,'_plus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  write.table(minus,paste0(outdir,ident,'_minus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}

for (comp in c('basal_v_intermediate','intermediate_v_umbrella','umbrella_v_basal')) { 
  DEXseq_res <- readRDS(paste0(outdir,comp,'_APAtest.res.rds'))
  raw.count <- DEXseq_res$counts.raw
  write.table(raw.count,paste0(outdir,comp,'_raw.count.txt'),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  res <- pairwise_res[[comp]]
  res <- res[,-c('dexseq_sig','ttest_sig','all_tests_sig')]
  write.table(res,paste0(outdir,comp,'_res.txt'),col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}

# Consensus tracks
setwd(outdir)
files <- list.files()
for (file in grep('res.txt$', files, value = TRUE)) {
  assign(gsub('_res.txt','',file),read.delim(file))
}

union_features <- Reduce(union, list(basal_v_intermediate$peak,intermediate_v_umbrella$peak,umbrella_v_basal$peak))

union_res <- matrix(nrow = length(union_features), ncol = 2) %>% as.data.frame()
colnames(union_res) <- c('TU','peak')
union_res$peak <- union_features
union_res$TU <- strsplit(union_res$peak,split=':') %>% sapply(.,'[[',1)

# All comparisons have the same set of features
for (var in c('basal_v_intermediate','intermediate_v_umbrella','umbrella_v_basal')) {
  comp <- strsplit(var,split = '_v_') %>% unlist()
  ident1 <- comp[1]
  ident2 <- comp[2]
  old <- get(var)
  old <- old[,c('tu','peak',paste0('mean_',ident1,'_frac'),paste0('mean_',ident2,'_frac'))]
  colnames(old)[3:4] <- paste0(c(ident1,ident2),'_',var)
  notfound <- union_features[!(union_features %in% old$peak)]
  old[(nrow(old)+1):(nrow(old)+length(notfound)),'peak'] <- notfound
  old[(nrow(old)+1):(nrow(old)+length(notfound)),'peak'] <- strsplit(notfound,split=':') %>% sapply(.,'[[',1)
  old <- old[match(union_features,old$peak),]
  stopifnot(identical(union_res$peak,old$peak))
  union_res <- cbind(union_res,old[,c(3,4)])
}

consensus_res <- union_res[,1:2]

for (celltype_name in c('basal','intermediate','umbrella')) {
  celltype_res <- union_res[,c(1:2,grep(paste0('^',celltype_name),colnames(union_res)))]
  celltype_res_spl <- split(celltype_res, celltype_res$TU)
  celltype_res_new <- lapply(celltype_res_spl,function(x) {
    tu_peak <- x[,1:2]
    all_na_check_col <- colSums(is.na(x)) == nrow(x)
    # Remove NA columns 
    x.na.rm <- x[,!all_na_check_col]
    all_numeric_check_col <- sapply(x.na.rm,is.numeric)
    x.value <- x.na.rm[,all_numeric_check_col] %>% as.data.frame()
    if (ncol(x.value) > 1) {
      # There are more than one column containing numeric values, pick whichever column has the highest value in the first row
      max.col.idx <- apply(x.value[1,],1,function(x) which.max(x)) %>% unname()
      x.value <- x.value[,max.col.idx]
      stopifnot(ncol(x.value) == 1)
    } else if (ncol(x.value) == 0) {x.value <- rep(0,times = nrow(tu_peak))}
    final_res <- cbind(tu_peak,x.value)
    colnames(final_res)[3] <- 'consensus'
    return(final_res)
  })
  
  celltype_res_new <- do.call('rbind',celltype_res_new)
  celltype_res_new <- celltype_res_new[match(consensus_res$peak,celltype_res_new$peak),]
  stopifnot(identical(consensus_res$peak,celltype_res_new$peak))
  celltype_res_new <- celltype_res_new[,ncol(celltype_res_new)] %>% as.data.frame()
  colnames(celltype_res_new) <- celltype_name
  consensus_res <- cbind(consensus_res,celltype_res_new)
}

peakref_subset <- subset(peak_ref, peak %in% consensus_res$peak) %>% dplyr::select(tu,peak,chr,start,end,strand)
peakref_subset <- peakref_subset[match(consensus_res$peak,peakref_subset$peak),]
stopifnot(identical(peakref_subset$peak,consensus_res$peak))
peakref_subset <- peakref_subset %>% dplyr::select(chr,start,end,strand)
consensus_res <- cbind(consensus_res,peakref_subset)

for (celltype_name in c('basal','intermediate','umbrella')) {
  bedGraph <- consensus_res[,c('chr','start','end','strand',celltype_name)]
  plus <- bedGraph[which(bedGraph$strand == '+'),]
  plus <- plus[,!(names(plus) %in% c('strand'))]
  minus <- bedGraph[which(bedGraph$strand == '-'),]
  minus <- minus[,!(names(minus) %in% c('strand'))]
  write.table(plus,paste0(outdir,celltype_name,'_merged_plus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  write.table(minus,paste0(outdir,celltype_name,'_merged_minus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
} 
