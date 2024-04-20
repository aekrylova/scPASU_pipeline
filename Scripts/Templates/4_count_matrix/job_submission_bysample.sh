#!/bin/bash

samples='sampleset'

for s in $samples; do
    sed "s/samplename/$s/g" 4a_BamToFastq.lsf | bsub
    #sed "s/samplename/$s/g" 4b_cellranger-genecount.lsf | bsub
    #sed "s/samplename/$s/g" 4b_cellranger-peakcount.lsf | bsub
done

