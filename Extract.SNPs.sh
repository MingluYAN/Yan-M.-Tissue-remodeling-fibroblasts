#Extraction of the RA-associated variants

awk '$2==11 && $3<128957453 && $3>127828656 && $9 < 5e-8 {print $0}' \
/path/to/RA_GWASmeta_TransEthnic_v2.txt \
> Transethnic_sig.txt

awk '{print "chr"$2":"$3"\t"$4"\t"$5}' \
Transethnic_sig.txt | sort -k 1,1 \
> Transethnic_sig_ID_allele.txt

awk '{print $2"\t"$5"\t"$6}' /path/to/\
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.no_multi.no_mono.EAS.CHR.POS.bim | \
sort -k 1,1 | join -j 1 Transethnic_sig_ID_allele.txt - \
> Transethnic_sig_ID_allele.check

awk '{print $1}' Transethnic_sig_ID_allele.txt \
> Transethnic_sig_ID_list.txt

/path/to/plink \
--threads 1 \
--memory 8000 \
--bfile /path/to/\
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.no_multi.no_mono.EAS.CHR.POS \
--extract Transethnic_sig_ID_list.txt \
--r2 'inter-chr' \
--ld-window-r2 0 \
--out EAS_r2_Transethnic_sig

/path/to/plink \
--threads 1 \
--memory 8000 \
--bfile /path/to/\
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.no_multi.no_mono.EUR.CHR.POS \
--extract Transethnic_sig_ID_list.txt \
--r2 'inter-chr' \
--ld-window-r2 0 \
--out EUR_r2_Transethnic_sig

/path/to/plink \
--threads 1 \
--memory 8000 \
--bfile /path/to/\
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.no_multi.no_mono.EAS.CHR.POS \
--show-tags Transethnic_sig_ID_list.txt \
--tag-r2 0.6 \
--tag-kb 500 \
--out EAS_r2_06_ETS1_sig_all.txt

/path/to/plink \
--threads 1 \
--memory 8000 \
/path/to/\
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.no_multi.no_mono.EUR.CHR.POS \
--show-tags Transethnic_sig_ID_list.txt \
--tag-r2 0.6 \
--tag-kb 500 \
--out EUR_r2_06_ETS1_sig_all.txt

sort EAS_r2_06_ETS1_sig_all.txt.tags > EAS_r2_06_ETS1_sig_all.txt
rm EAS_r2_06_ETS1_sig_all.txt.tags

sort EUR_r2_06_ETS1_sig_all.txt.tags > EUR_r2_06_ETS1_sig_all.txt
rm EUR_r2_06_ETS1_sig_all.txt.tags

join EAS_r2_06_ETS1_sig_all.txt EUR_r2_06_ETS1_sig_all.txt \
> EAS_EUR_r2_06_ETS1_sig_all.txt

#lift over

cat EAS_EUR_r2_06_ETS1_sig_all.txt | \
sed 's/:/\t/g' | awk '{print $1"\t"$2-1"\t"$2}' \
> EAS_EUR_r2_06_ETS1_sig_all.hg19.bed

/path/to/liftOver \
EAS_EUR_r2_06_ETS1_sig_all.hg19.bed \
/path/to/hg19ToHg38.over.chain.gz \
EAS_EUR_r2_06_ETS1_sig_all.hg38.bed EAS_EUR_r2_06_ETS1_sig_all.unmapped.bed

awk '{print $1":"$3}' EAS_EUR_r2_06_ETS1_sig_all.hg38.bed \
> EAS_EUR_r2_06_ETS1_sig_all.hg38.txt

#join with allele info in GWAS
paste -d "\t" EAS_EUR_r2_06_ETS1_sig_all.txt EAS_EUR_r2_06_ETS1_sig_all.hg38.txt \
> EAS_EUR_r2_06_ETS1_sig_all.hg19_38.txt

awk '{print $2"\t"$5"\t"$6}' /path/to/\
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.no_multi.no_mono.EAS.CHR.POS.bim | \
sort -k 1,1 | join -j 1 EAS_EUR_r2_06_ETS1_sig_all.hg19_38.txt - | \
sed '1ihg19\thg38\tA1\tA2' > EAS_EUR_r2_06_ETS1_sig_all.allele