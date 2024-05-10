#!/usr/bin/bash

#BSUB -W 24:00
#BSUB -n 32
#BSUB -M 64
#BSUB -R rusage[mem=64]
#BSUB -o /rsrch6/scratch/mol_cgenesis/bnle/bsub_jobs/output/test_3e_identify_fragmented_peaks.%J.out
#BSUB -e /rsrch6/scratch/mol_cgenesis/bnle/bsub_jobs/error/test_3e_identify_fragmented_peaks.%J.err
#BSUB -q e40medium
#BSUB -J test_3e_identify_fragmented_peaks
#BSUB -u bnle@mdanderson.org
#BSUB -B
#BSUB -N

source activate scPASU_env
module load bedtools/2.30.0
module load bowtie2

# Variable
dir=/rsrch6/scratch/mol_cgenesis/bnle/test_scPASU_run/outputs/
compartment=immune
dir=${dir}${compartment}/
script_dir=/rsrch6/scratch/mol_cgenesis/bnle/test_scPASU_run/Scripts/immune
script_dir=${script_dir}/3_peak_ref_cleanup/
ref_genome=/rsrch6/scratch/mol_cgenesis/bnle/GRCh38/fasta/genome.fa
in=${dir}3g_peak_classification/
out=${dir}3h_fragmented_peaks_to_merge/
bamdir=${dir}1e_merged_bam/
refdir=${dir}3f_merge_two_prongs/
ucsc_utilities_dir=/rsrch6/scratch/mol_cgenesis/bnle/

# Argument
PAS_filtering=true
binsize=1
gap_threshold=0
min_coverage=10
spliced_read_pct_thres=40
realign_peak_read_pct_thres=40
ncores=32

if [ $PAS_filtering == 'true' ]
   then fprefix=${compartment}
else
   fprefix=${compartment}_noPASfiltering
fi

peak_saf=${refdir}${fprefix}_4thtu_assigned_peak_universe_updated.saf

# Functions
check_empty_files() {
    local files=("$@")
    local is_empty=1

    for file in "${files[@]}"; do
        if [ ! -s "$file" ]; then
            is_empty=0
            break
        fi
    done

    return $is_empty
}

mkdir -p ${out}

transcripts=${in}${fprefix}_transcripts_causing_peak_fragmentation.txt
if [ $(tail -n+2 $transcripts | wc -l) -eq 0 ]; then
    echo "There are no transcripts likely causing fragmentation. Exiting script."
    exit 1
fi
echo Extract transcripts likely causing peak fragmentation...

for strand in plus minus; do
   bed_input=${out}${fprefix}_transcripts_causing_peak_fragmentation_exons_${strand}.bed
   fa_output=${out}${fprefix}_transcripts_causing_peak_fragmentation_exons_${strand}.fa

   if [ $strand == 'plus' ]; then 
     tail -n+2 $transcripts | awk -v OFS='\t' '{print $1,$2,$3,$10,0,$5}' | awk '$6 =="+" {print}' > $bed_input
   elif [ $strand == 'minus' ]; then 
     tail -n+2 $transcripts | awk -v OFS='\t' '{print $1,$2,$3,$10,0,$5}' | awk '$6 =="-" {print}' > $bed_input
   fi

   echo Extract exons of these transcripts from $ref_genome
   bedtools getfasta -fi $ref_genome -bed $bed_input -s -nameOnly -fo $fa_output
   echo Stitch the exons to obtain these transcripts
   Rscript ${script_dir}3e_stitch_exons.R -i $fa_output -o ${out}${fprefix}_transcripts_causing_peak_fragmentation_${strand}.fa
done

echo Filter bam file for spliced reads overlapping with merge candidates...

cd $out
for strand in plus minus; do
   bamfile=${bamdir}dedup_uniq_genomicAfiltered_merged_${strand}.bam
   peak_bed=${refdir}${fprefix}_4thtu_assigned_${strand}.bed
   subset_peak_bed=${fprefix}_4thtu_assigned_${strand}_subset.bed
   subset_peak_saf=${fprefix}_4thtu_assigned_${strand}_subset.saf
   subset_bam_file=${fprefix}_dedup_uniq_genomicAfiltered_${strand}_subset.bam

   if [ $strand == 'plus' ]; then 
      tail -n+2 $transcripts | awk '$5 == "+"{print}' | cut -f11 | grep -v '^[[:space:]]*$' | tr ',' '\n' | sort | uniq > ${subset_peak_bed}_temp
   elif [ $strand == 'minus' ]; then 
      tail -n+2 $transcripts | awk '$5 == "-"{print}' | cut -f11 | grep -v '^[[:space:]]*$' | tr ',' '\n' | sort | uniq > ${subset_peak_bed}_temp
   fi

   grep -w -f ${subset_peak_bed}_temp $peak_bed > $subset_peak_bed
   echo GeneID$'\t'Chr$'\t'Start$'\t'End$'\t'Strand > $subset_peak_saf
   grep -w -f ${subset_peak_bed}_temp $peak_saf >> $subset_peak_saf

   echo Count reads for merge candidates...
   samtools view -b -h -L $subset_peak_bed $bamfile > $subset_bam_file
   Rscript ${script_dir}feature_counts.R -r $subset_peak_saf -b $subset_bam_file -o $out -f ${fprefix}_${strand}_subset -c $ncores -i no
   echo Count spliced reads for merge candidates...
   spliced_only_sam_temp=${subset_bam_file%.bam}_spliced_only.sam
   samtools view -H $subset_bam_file > $spliced_only_sam_temp
   samtools view $subset_bam_file | awk '($6 ~ /N/)' >> $spliced_only_sam_temp
   samtools view -bS $spliced_only_sam_temp > ${spliced_only_sam_temp%.sam}.bam
   samtools view $subset_bam_file | awk '($6 ~ /N/)' | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${subset_bam_file%.bam}_spliced_only.fastq   
   Rscript ${script_dir}feature_counts.R -r $subset_peak_saf -b ${spliced_only_sam_temp%.sam}.bam -o $out -f ${fprefix}_${strand}_subset_spliced_reads -c $ncores -i no
   
   rm ${subset_peak_bed}_temp $spliced_only_sam_temp $subset_bam_file
done

echo Call peaks on the realignment of spliced reads from merge candidates against their likely exomes...
for strand in plus minus; do
   ### Realign spliced reads from merge candidates against transcripts likely causing peak fragmentation
   bowtie2_prefix=${fprefix}_transcripts_causing_peak_fragmentation_${strand}
   realign_prefix=${fprefix}_realign_merge_candidates_${strand}
   
   if ! check_empty_files ${bowtie2_prefix}.fa; then
       bowtie2-build ${bowtie2_prefix}.fa ${bowtie2_prefix}
       fastq_file=${fprefix}_dedup_uniq_genomicAfiltered_${strand}_subset_spliced_only.fastq
       bowtie2 --local --threads $ncores -x ${bowtie2_prefix} -X 2000 -U ${fastq_file} -S ${realign_prefix}.sam
   
       ### Preparing files for IGV 
       samtools faidx ${bowtie2_prefix}.fa
       samtools view -bS ${realign_prefix}.sam > ${realign_prefix}.bam
       samtools sort ${realign_prefix}.bam -@ $ncores -o ${realign_prefix}.sorted.bam > ${realign_prefix}.sorted.bam
       samtools index ${realign_prefix}.sorted.bam
       rm ${realign_prefix}.bam ${realign_prefix}.sam

       ### Call peaks
       bamCoverage -b ${realign_prefix}.sorted.bam -o ${realign_prefix}.bw -p $ncores -v --binSize ${binsize}
       ${ucsc_utilities_dir}/bigWigToBedGraph ${realign_prefix}.bw ${realign_prefix}.bedGraph
   else
      echo chr$'\t'start$'\t'end$'\t'score > ${realign_prefix}.bedGraph
   fi 
done

Rscript ${script_dir}3e_peak_calling_from_bedGraph.R -f $fprefix -b $out -o $out -g $gap_threshold -m $min_coverage

realign_peak=${fprefix}_spliced_reads_realign_peaks.txt
for strand in plus minus; do   
   ### Count reads for realign peaks from spliced reads
   saf=${fprefix}_spliced_reads_realign_peaks_${strand}.saf
   realign_prefix=${fprefix}_realign_merge_candidates_${strand}
   bamfile=${realign_prefix}.sorted.bam
   echo GeneID$'\t'Chr$'\t'Start$'\t'End$'\t'Strand > $saf
   if [ $strand == 'plus' ]; then
      awk -v strand="+" -v OFS='\t' '$6 == strand { print $4,$1,$2,$3,$6}' $realign_peak >> $saf
   elif [ $strand == 'minus' ]; then 
      awk -v strand="-" -v OFS='\t' '$6 == strand { print $4,$1,$2,$3,$6}' $realign_peak >> $saf
   fi
   Rscript ${script_dir}feature_counts.R -r $saf -b $bamfile -o $out -f ${fprefix}_spliced_reads_realign_peaks_${strand} -c $ncores -i no
done

echo Identify fragmented peaks to merge...
ref=${in}${fprefix}_peak_universe_classified.txt 
Rscript ${script_dir}3e_identify_fragmented_peaks.R -t $transcripts -p $realign_peak -d $out -f $fprefix -o $out -s $spliced_read_pct_thres -r $realign_peak_read_pct_thres -b $ref -n $ncores
