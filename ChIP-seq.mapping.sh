
#!/bin/bash
SAMPLE=$1

/path/to/bowtie2 \
-p 4 \
-x /path/to/bowtie_index/hg19 \
-q /path/to/${SAMPLE}.trimmed.fastq.gz \
-S ChIP_${SAMPLE}.sam

grep -v "XS:" ChIP_${SAMPLE}.sam > ChIP_${SAMPLE}.unique.sam
rm ChIP_${SAMPLE}.sam
samtools view -@ 4 -h -Sb ChIP_${SAMPLE}.unique.sam | samtools sort -@ 4 > ChIP_${SAMPLE}.bam
samtools index ChIP_${SAMPLE}.bam
rm ChIP_${SAMPLE}.unique.sam

samtools view -@ 4 -o ChIP_${SAMPLE}.autosome.bam ChIP_${SAMPLE}.bam `seq 1 22 | sed 's/^/chr/g'`
samtools index ChIP_${SAMPLE}.autosome.bam

rm ChIP_${SAMPLE}.bam.bai
rm ChIP_${SAMPLE}.bam

bedtools intersect -v -abam ChIP_${SAMPLE}.autosome.bam \
-b /path/to/hg19-blacklist.v2.bed \
> ChIP_${SAMPLE}.filtered.bam
samtools index ChIP_${SAMPLE}.filtered.bam

rm ChIP_${SAMPLE}.autosome.bam.bai
rm ChIP_${SAMPLE}.autosome.bam

mkdir -p stats

samtools stats ChIP_${SAMPLE}.filtered.bam \
> stats/ChIP_${SAMPLE}.filtered.stats

bamCoverage -b ChIP_${SAMPLE}.filtered.bam -o ChIP_${SAMPLE}.bw \
--normalizeUsing CPM