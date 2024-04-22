library(data.table)
library(dplyr)
library(ggplot2)
library(goldmine)
library(RColorBrewer)
library(gridExtra)

outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3h_fragmented_peaks_to_merge/'
fprefix <- 'Ureter10_urothelial'
outdir <- '~/Desktop/Ureter10_scPASU_run/'
fprefix <- 'atlas.clusters.2.0.GRCh38.96_minGapwidth1'
  
ref <- fread(paste0(outdir,'Ureter10_urothelial_final_peak_universe_updated.txt'))
ref <- fread('~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96_minGapwidth1_peak_universe_updated.txt')

ref_tu_tally <- ref %>% group_by(tu) %>% tally() %>% as.data.frame()

turef <- readRDS('~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3f_merge_two_prongs/genes_newflankupdated.rds')
turef <- readRDS('~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96_minGapwidth1_genes_flankupdated_new.rds')

# AHT colors
"#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9",
"#74add1","#4575b4","#313695","#762a83","#c2a5cF","#969696","#525252","black"
"#a50026","#fdae61","#abd9e9","#313695"

nTU <- nrow(ref_tu_tally)
nPeak <- sum(ref_tu_tally$n)
# Build dataset with different distributions
p <- ref_tu_tally %>%
  ggplot( aes(x=n)) +
  geom_histogram( binwidth=1, fill="#969696", color="black", alpha=0.9, linewidth = 0.3) +
  labs(subtitle= paste0("Number of TUs = ",nTU,"\nNumber of poly(A) sites = ",nPeak,"\n"), x = "Poly(A) sites per TU", y = "Frequency") +
  theme_classic()

p + theme(text = element_text(size = 14))            
ggsave(paste0(outdir,fprefix,'_peaknum_per_TU.png'),width = 7, height = 4, unit = 'in')

### Metagene plot 
tu <- as.data.frame(turef$tu)
flank <- as.data.frame(turef$flank)

turef_flank <- tu
turef_flank$flank_start <- flank[match(turef_flank$tu,flank$tu),'start']
turef_flank$flank_end <- flank[match(turef_flank$tu,flank$tu),'end']

turef_flank$end <- ifelse(turef_flank$strand == '+',turef_flank$flank_end,turef_flank$end)
turef_flank$start <- ifelse(turef_flank$strand == '-',turef_flank$flank_start,turef_flank$start)  
turef_flank$width <- turef_flank$end - turef_flank$start

ref_mid_pt <- ref %>% mutate(middle = (pr_start + round((pr_end - pr_start)/2) - 1))

ref_mid_pt_spl <- split(ref_mid_pt,ref_mid_pt$tu)

norm_mid_pt <- lapply(ref_mid_pt_spl,function(x) {
  tu_id <- unique(x$tu)
  peaks <- x$final_annotation
  tu_coord <- subset(turef_flank, tu == tu_id)
  width <- tu_coord$width
  strand <- unique(x$strand)
  if (strand == '+') {
    tss <- tu_coord$start
    norm_mid <- (x$middle - tss)/width
  } else if (strand == '-') {
    tss <- tu_coord$end
    norm_mid <- (tss - x$middle)/width
  } 
  stopifnot(length(peaks)==length(norm_mid))
  df <- data.frame(peak = peaks, norm_mid = norm_mid, tu = tu_id, strand = strand)
  return(df)
})

norm_mid_pt_df <- do.call('rbind',norm_mid_pt)

#nTU_wirange <- subset(norm_mid_pt_df, norm_mid >= -0.1)
p <- norm_mid_pt_df %>%
  ggplot( aes(x=norm_mid)) +
  geom_histogram(binwidth=0.025, fill="#969696", color="black", alpha=0.9) +
  coord_cartesian(xlim = c(-0.1,1.0)) +
  scale_x_continuous(breaks=c(0,0.2,0.4,0.8,0.6,1.0))+
  labs(subtitle= paste0("Number of TUs = ",length(unique(norm_mid_pt_df$tu)),"\nNumber of peaks = ",nrow(norm_mid_pt_df),"\n"), x = "Poly(A) Site Position", y = "Frequency") +
  theme_classic()
  
p + theme(text = element_text(size = 24))     

p1 <- p + theme(text = element_text(size = 24))
p2 <- p + theme(text = element_text(size = 24))

p1_gtable <- ggplot_gtable(ggplot_build(p1))
p2_gtable <- ggplot_gtable(ggplot_build(p2))

p1_gtable$widths <- p2_gtable$widths

png(paste0(outdir,'Ureter10_urothelial_metagene.png'),width = 10, height = 5, unit = 'in', res = 1200)
grid::grid.draw(p1_gtable)
dev.off()

png(paste0(outdir,'atlas.clusters.2.0.GRCh38.96_minGapwidth1_metagene.png'),width = 10, height = 5, unit = 'in', res = 1200)
grid::grid.draw(p2_gtable)
dev.off()

### PR Width ###

p <- ref %>%
  ggplot( aes(x=pr_width)) +
  geom_histogram(binwidth=10, fill="#969696", color="black", alpha=0.9, linewidth = 0.3) +
  labs(title = fprefix, subtitle= paste0("Number of TUs = ",nTU,"\nNumber of peaks = ",nPeak,"\n"), x = "PR Width", y = "Frequency") +
#  scale_x_continuous(breaks = seq(0, 250, 50)) +
  theme_classic()

p + theme(text = element_text(size = 14))       
ggsave(paste0(outdir,fprefix,'_pr_width.png'),width = 6, height = 4, unit = 'in',p)


p <- ref %>%
  ggplot( aes(x=peak_width)) +
  geom_histogram(binwidth=10, fill="#969696", color="black", alpha=0.9, linewidth = 0.3) +
  labs(title = fprefix, subtitle= paste0("Number of TUs = ",nTU,"\nNumber of peaks = ",nPeak,"\n"), x = "Peak Width", y = "Frequency") +
  #  scale_x_continuous(breaks = seq(0, 250, 50)) +
  theme_classic()

p + theme(text = element_text(size = 14))       
ggsave(paste0(outdir,fprefix,'_peak_width.png'),width = 6, height = 4, unit = 'in',p)


##### polyADb
chr<-paste0('chr',1:22)
chr<-c(chr,'chrX','chrY')

polyadb_fpath_hg38<-'~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96.tsv'
pdb_hg38 <- fread(polyadb_fpath_hg38)
pdb_hg38$chrom<-paste0('chr',pdb_hg38$chrom)
pdb_hg38<-pdb_hg38[pdb_hg38$chrom %in% chr,]

# TE, terminal exon; EX, exonic; IN, intronic; DS, 1,000 nt downstream of an annotated terminal exon
# AE, anti-sense to an exon; AI, anti-sense to an intron; AU, 1,000 nt upstream in anti-sense direction of a transcription start site; IG, intergenic
source('~/Desktop/Ureter10_scPASU_run/Scripts/3_peak_ref_cleanup/scPASU_functions.R')
pdb_hg38 <- pdb_hg38[,c(1,2,3,6)]
colnames(pdb_hg38) <- c('chr','start','end','strand')
pdb_gr <- makeGRanges(pdb_hg38, strand = T)
pdb_gr<-reduce(pdb_gr,min.gapwidth=1)
pdb_hg38 <- as.data.frame(pdb_gr)
pdb_hg38$width <- NULL
colnames(pdb_hg38) <- c('chr','start','end','strand')

pdb_gr <- makeGRanges(pdb_hg38,strand=T)
pdb_gr$polya <- 'polya'
seqlevels(pdb_gr) <- chr
pdb_gr <- pdb_gr[order(pdb_gr)]
pdb_gr$peakID <- paste0('peak_',1:length(pdb_gr))

turef <- readRDS('~/Desktop/Ureter10_scPASU_run/outputs/urothelial/hg38_turef/genes.rds')
jtu <- joinTus_peaks(pdb_gr,turef)

save(jtu,pdb_hg38,turef,file='~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96_minGapwidth1_jtu.rds')

load('~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96_minGapwidth1_jtu.rds')

# Create Peak reference with all relevant columns merged after TU assignment #
cat('Add other relevant cols \n')

jtu$join<-as.data.frame(jtu$join)
# Remove multi-TU peaks
cat('Remove peaks that are assigned to multiple TUs, i.e. non-unique_peak TUs \n')
r1<-which(jtu$join$unique_peak==FALSE)

if(length(r1)!=0) {
  tu_tab<-jtu$join[-r1,] 
}else{
  tu_tab<-jtu$join
}

# Unique TU (TU contains only peaks assigned uniquely to TU) status is determined before removing peaks with unique_peak==FALSE so there can still be FALSE unique_tu peaks left

# Create peak per TU count
tu_peak<-jtu$join %>% dplyr::select(peak,tu) %>% group_by(tu) %>% tally()

# Peaks (saved from merging all peaks from MACS2 output)
peaks<-as.data.frame(jtu$polya_peaks)
found<-match(tu_tab$peak,peaks$peakID)
cat(paste(length(found),'peaks overlap uniquely with the reference TU. Only move forward with these peaks \n'))
matched_peaks<-peaks[found,]

# Sort all tables
# tu_tab<-tu_tab[mixedorder(tu_tab$peak),] Peak numbering, in this new approach, is not necessarily reflective of their relative position
# matched_peaks<-matched_peaks[mixedorder(matched_peaks$peakID),]

# Now lets transfer peak info to TU table
stopifnot(identical(tu_tab$peak,matched_peaks$peakID)) # TRUE

tu_tab$chr<-matched_peaks$seqnames
tu_tab$start<-matched_peaks$start
tu_tab$end<-matched_peaks$end
tu_tab$strand<-matched_peaks$strand
tu_tab$gene<-strsplit(tu_tab$tu_anno,split=':') %>% sapply(.,'[[',2)

# Also transfer other info for future use. This table can then serve as comprehensive features table

tu_tab$pr_chr<-matched_peaks$seqnames
tu_tab$pr_start<-matched_peaks$start
tu_tab$pr_end<-matched_peaks$end
tu_tab$pr_strand<-matched_peaks$strand
tu_tab$pr_width<-matched_peaks$end - matched_peaks$start
tu_tab$peakID<-matched_peaks$peakID

# Now add peakID to tu annotation column (multi peak TU will look like TU1:gene:P1, TU1:gene:P2 and TU1:gene:P3 while single peak TU will look like TU2:gene:P0)
cat('Sorting peaks by genomic coordinates before peak numbering \n')
tu_tab <- makeGRanges(tu_tab, strand = T)
seqlevels(tu_tab) <- chr
tu_tab <- tu_tab[order(tu_tab)]
tu_tab <- as.data.frame(tu_tab)
tu_tab$width <- NULL
colnames(tu_tab)[1] <- 'chr'

# Get TU count
tu_count<-tu_tab %>% group_by(tu) %>% tally()

# Sort
tu_count<-tu_count[order(tu_count$n),]

#Single peak per TU
single<-which(tu_count$n==1)
s<-match(tu_count$tu[single],tu_tab$tu)
single_tu<-tu_tab[s,]
cat(paste(single_tu$tu %>% unique() %>% length(),'TUs have only one peaks. These peaks are hence annotated P0 \n'))
single_tu$final_annotation<-paste0(single_tu$tu_anno,':P0')

# Multi peak TU -
multi_tu<-tu_tab[-s,]

# Split by strand
multi_tu_p<-multi_tu[multi_tu$strand=='+',]
multi_tu_m<-multi_tu[multi_tu$strand=='-',]

# Now add final annotation col
multi_tu_p2<-create_final_annotation_col(multi_tu_p)
multi_tu_m2<-create_final_annotation_col(multi_tu_m,is_minus=TRUE)

# Merge plus and minus strand
multi_tu2<-rbind(multi_tu_p2,multi_tu_m2)

# Add P0 genes back to ref
merged_tu3<-rbind(single_tu,multi_tu2)

# Save
merged_tu3<-merged_tu3 %>% as.data.frame()

write.table(merged_tu3,'~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96_minGapwidth1_peak_universe_updated.txt',
            sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

turef_flank_updated <- update_flank(merged_tu3,turef,save=FALSE)

saveRDS(turef_flank_updated,'~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96_minGapwidth1_genes_flankupdated_new.rds')
