#ChIP-seq analysis 
#1 env setting

cd /path/to/analysis_dir

#2 download files
mkdir -p fastqc
awk '{print $2}' fastq_name_ID_list.txt | while read line ; do
/path/to/fastqc  \
-o fastqc ${line}.fastq.gz
done
cd fastqc
multiqc .

#QC by fastp
#adaptor trimming + quality filtering + length filtering
mkdir -p /path/to/analysis_dir/fastp
cd /path/to/analysis_dir/fastp
mkdir -p fastp_html
awk '{print $2}' ../raw/fastq_name_ID_list.txt | while read line ; do
fastp -i ../raw/${line}.fastq.gz \
-o ${line}.trimmed.fastq.gz \
-3 -5 \
-W 4 \
-M 20 \
--adapter_fasta /path/to/TruSeq3-SE.fa \
-w 4 \
-q 30 \
-n 5 \
-l 25 \
--html fastp_html/${line}.fastp.html
done

cd /path/to/analysis_dir/fastp
multiqc .

mkdir -p fastqc
awk '{print $2}' ../raw/fastq_name_ID_list.txt | while read line ; do
/path/to/fastqc  \
-o fastqc ${line}.trimmed.fastq.gz
done
cd fastqc
multiqc .
