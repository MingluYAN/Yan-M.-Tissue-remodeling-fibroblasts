#!/bin/bash
SAMPLE=$1
BAM_DIR=$2
BAM=$3

THR=0.01

macs2 callpeak -t ${BAM_DIR}/${BAM} \
 	-f BAM -g 2.8e+9 \
    -q ${THR} \
    --keep-dup all \
	--outdir default/${SAMPLE} \
    -n ${SAMPLE}_${THR}_default 2> log/${SAMPLE}_${THR}_default.log