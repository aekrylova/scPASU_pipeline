library(dplyr)
library(goldmine)
library(argparser)

parser<-arg_parser(name="ref2gtf.R",description="Convert final poly(A) site reference to gtf file format.")

parser<-add_argument(
  parser,
  arg='--peaks',
  short = '-b',
  type="character",
  default=NA,
  help="Enter the path to the peak file from the last filtering step. Format: path/to/peaks_file")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

args <- parse_args(parser)

peaks_file<- args$peaks
fprefix <- args$file_prefix
outdir <- args$outdir

ref <- read.delim(peaks_file)
ref_gr <- makeGRanges(ref, strand = T)
ref_gr <- ref_gr[order(ref_gr)]
ref <- as.data.frame(ref_gr)
ref$width <- NULL
colnames(ref)[1] <- 'chr'

ref_gtf <- matrix(nrow = nrow(ref), ncol = 9) %>% as.data.frame()
ref_gtf$V1 <- ref$chr
ref_gtf$V2 <- 'unknown'
ref_gtf$V3 <- 'exon'
ref_gtf$V4 <- ref$start
ref_gtf$V5 <- ref$end
ref_gtf$V6 <- '.'
ref_gtf$V7 <- ref$strand
ref_gtf$V8 <- '.'
ref_gtf$V9 <- paste0('gene_id "',ref$final_annotation,'"; transcript_id "',ref$final_annotation,'"; gene_name "',ref$final_annotation,'"; gene_biotype "',ref$final_classification,'";')

write.table(ref_gtf,paste0(outdir,fprefix,'_final_peak_universe_updated.gtf'),quote = F, sep = '\t', row.names = F, col.names = F)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
