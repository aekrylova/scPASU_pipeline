# This code is modified from Ting lab's bulk APA pipeline to repurpose it for scRNA APA analysis
# Author: Surbhi Sona

# original script from bulkAPA pipeline: '08_pas_stats_figure.r'
# hg38 polya db: atlas.clusters.2.0.GRCh38.96.bed  (downloaded from :https://polyasite.unibas.ch/atlas#2)
# hg19 db: polyAsite_v1_hg19_clusters.bed


source('~/Desktop/Scripts/scPASU_pipeline/scPASU_plotting_functions.R')

#####################################################################

### Step 1. Prep data ###

#only PR, PR upto 100 bases, polyAdb and NULL

out_dir<-'~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3f_merge_two_prongs/hexamer_plot/'
tmp_dir<- paste0(out_dir,'tmp/')
fprefix <-'Ureter10_urothelial'

if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
if (!dir.exists(tmp_dir)) {dir.create(tmp_dir, recursive = TRUE)}

# You can find these file in the Tinglab directory Michael copied from CCF
polyadb_fpath_hg38<-'~/Desktop/Ureter10_scPASU_run/atlas.clusters.2.0.GRCh38.96.bed'

# This should be your final peak reference file (check this file for format and update code if you have different column names
f<-'~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3f_merge_two_prongs/Ureter10_urothelial_4thtu_assigned_peak_universe_updated.txt'

file<-fread(f,header=TRUE)

## (A) For PR ##

# Please update column names if they differ in your table from the ones listed below

file2<-file
colnames(file2)<-gsub('\\bchr\\b','peak_chr',colnames(file2))
colnames(file2)<-gsub('\\bstart\\b','peak_start',colnames(file2))

colnames(file2)<-gsub('\\bend\\b','peak_end',colnames(file2))
#colnames(file2)<-gsub('\\bstrand\\b','peak_strand',colnames(file2))

colnames(file2)<-gsub('\\bpr_chr\\b','chr',colnames(file2))
colnames(file2)<-gsub('\\bpr_start\\b','start',colnames(file2))
colnames(file2)<-gsub('\\bpr_end\\b','end',colnames(file2))
#colnames(file2)<-gsub('\\bpr_strand\\b','strand',colnames(file2))

gr2<-makeGRangesFromDataFrame(file2,ignore.strand = FALSE)

pr_hex<- findHex(gr2,us=50)

## (B) polyA db hg38 ##

chr<-paste0('chr',1:22)
chr<-c(chr,'chrX','chrY')

pdb_hg38 <- fread(polyadb_fpath_hg38)
pdb_hg38$V1<-paste0('chr',pdb_hg38$V1)
pdb_hg38<-pdb_hg38[pdb_hg38$V1 %in% chr,]

# Create GR range
pa_hg38 <- with(pdb_hg38,GRanges(V1,IRanges(V2,V3),strand=V6))
pa_hg38 <- pa_hg38[seqnames(pa_hg38) %in% chr]
pa_hg38 <- unique(pa_hg38)

#pa_hg38_hex <- countHexes(gr=pa_hg38,us=50)
#pa_hg38_hex_up <- countHexes_up(gr=pa_hg38,us=50)

# cluster pA within 1 bases of each other - skip the reduce() based on your workflow (I used this code in my thesis difure)
pa_hg38_clus<-reduce(pa_hg38,min.gapwidth=1)
pa_hg38_clus_hex <- findHex(gr=pa_hg38_clus,us=50)

## (C) PR upto 100 bases ##

# All longer PRs are filtered out in this set

file2$pr_width %>% range()
select<-which(file2$pr_width<=100)
file2_sub<-file2[select,]
pr_gr100<-makeGRangesFromDataFrame(file2_sub,ignore.strand = FALSE)

pr_hex100<-findHex(pr_gr100,us=50)

## (D) NULL ##

dgp_hg38_100 <- drawGenomePool(query=pr_hex100,
                               n=10,
                               genome="hg38",
                               cachedir=tmp_dir,
                               sync=FALSE)

null_hg38_100 <- findHex(gr=dgp_hg38_100,us = 50)

############################################################

### Step2. Collect hex counts ###

hexlist<-list()
hexlist$pr<-pr_hex
hexlist$pr100<-pr_hex100
hexlist$pdb_hg38_clus<-pa_hg38_clus_hex
hexlist$null<-null_hg38_100


# Make tables
hextabs <- lapply(hexlist,function(x) as.data.frame(table(x$hex)))
for(h in 1:length(hextabs))
{
  hextabs[[h]]$run <- names(hexlist)[h]
}

tab <- rbindlist(hextabs)

####################################################################

### Step 3. Generate plots ###
#cols <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#969696","#666666")
cols <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#abd9e9",
          "#74add1","#4575b4","#313695","#762a83","#c2a5cF","#969696","#525252")
hexlist_names <- names(hexlist)
L <- length(hexlist_names)

sample_names <- hexlist_names

tab.per <- tab[,list(Var1=Var1,Freq=Freq/sum(Freq)),by="run"]
#tab.per[run=="background"]
tab.per$run <- factor(tab.per$run,levels=c(sample_names))
tab.per$motif <- tab.per$Var1

library(ggplot2)


#ggplot(tab,aes(x=run,y=Freq,fill=Var1)) + geom_bar(stat="identity") + handy::ggnice() + scale_fill_manual(values=cols)
p <- ggplot(tab.per,aes(x=run,y=Freq,fill=motif)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cols) +
  labs(y="Fraction of Regions",x="") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme_classic()

ggsave(paste0(out_dir,fprefix,'_hexamer_plot.png'),width = 8, height = 8, units = 'in', p, bg = 'white')

saveRDS(hexlist,paste0(out_dir,fprefix,'_hexamer_find_res.rds'),)
