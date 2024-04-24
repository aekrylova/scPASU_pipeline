#!/bin/bash

samples='sampleset'

for s in $samples; do
    sed "s/samplename/$s/g" 1b_preprocess_bam_dedup_cleanup.lsf | bsub
    #sed "s/samplename/$s/g" 1b_preprocess_bam_genomicA_filtering_1stround.lsf | bsub
    #sed "s/samplename/$s/g" 1b_preprocess_bam_genomicA_filtering_2ndround.lsf | bsub
done
