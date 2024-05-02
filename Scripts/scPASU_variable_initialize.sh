#!/usr/bin/bash

runname="Ureter10"
compartmentname="urothelial"
sampleset="U1 U2 U3 U4 U5 U6 U7 U8 U9 U10"
polyAfilterpath="/rsrch6/scratch/mol_cgenesis/bnle/polyAfilter/"
dirpath="/rsrch6/scratch/mol_cgenesis/bnle/${runname}_scPASU_run/outputs/"
scriptpath="/rsrch6/scratch/mol_cgenesis/bnle/${runname}_scPASU_run/Scripts/${compartmentname}"
genomename="hg38"
cellrangeroutputpath="/rsrch6/scratch/mol_cgenesis/bnle/Ureter10_data/cellranger_outputs_GRCh38-2023/"
cellrangerrefpath="/rsrch6/scratch/mol_cgenesis/bnle/GRCh38/"
fastapath="/rsrch6/scratch/mol_cgenesis/bnle/GRCh38/fasta/genome.fa"
gtfpath="/rsrch6/scratch/mol_cgenesis/bnle/GRCh38/genes/genes.gtf"
subsetbampath="/rsrch6/scratch/mol_cgenesis/bnle/"
chrsname="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
picardpath=/risapps/noarch/picard/2.23.8/build/libs/
bedtoolspath=/risapps/rhel7/bedtools/2.30.0/bin/
bowtie2path=/risapps/noarch/bowtie2/2.4.1/
ucscutilitiespath=/rsrch6/scratch/mol_cgenesis/bnle/
bamtofastqpath=/risapps/rhel7/cellranger/7.1.0/lib/bin/
cellrangerpath=/risapps/rhel7/cellranger/7.1.0/ 

mkdir -p ${scriptpath}/1_data_preprocessing/ ${scriptpath}/2_peak_ref/ ${scriptpath}/3_peak_ref_cleanup/ ${scriptpath}/4_count_matrix/ ${scriptpath}/7_ucsc_uploads/
templates="Templates/1_data_preprocessing/* Templates/2_peak_ref/* Templates/3_peak_ref_cleanup/* Templates/4_count_matrix/* Templates/7_ucsc_uploads/*"

for f in $templates; do
filename=$(basename $f)
dir=$(dirname $f)
dir=$(basename $dir)
sed "s|runname|$runname|g" $f | sed "s|compartmentname|$compartmentname|g" | sed "s|sampleset|$sampleset|g" | sed "s|polyAfilterpath|$polyAfilterpath|g" | sed "s|dirpath|$dirpath|g" | sed "s|scriptpath|$scriptpath|g" | sed "s|cellrangeroutputpath|$cellrangeroutputpath|g" | sed "s|cellrangerrefpath|$cellrangerrefpath|g" | sed "s|fastapath|$fastapath|g" | sed "s|genomename|$genomename|g" | sed "s|gtfpath|$gtfpath|g" | sed "s|subsetbampath|$subsetbampath|g" | sed "s|chrsname|$chrsname|g" | sed "s|picardpath|$picardpath|g" | sed "s|bedtoolspath|$bedtoolspath|g" | sed "s|bowtie2path|$bowtie2path|g" |sed "s|ucscutilitiespath|$ucscutilitiespath|g" | sed "s|bamtofastqpath|$bamtofastqpath|g" | sed "s|cellrangerpath|$cellrangerpath|g" > ${scriptpath}/${dir}/${filename}
done

cp -r Templates/5_APA_testing/ Templates/pa_site_ref_qc_plots/ ${scriptpath}/
