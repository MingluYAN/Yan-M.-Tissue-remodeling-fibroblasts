#Code for Fig Ex.3f
library(tidyverse)
library(RColorBrewer)

SNP<-read_tsv("EAS_EUR_r2_06_ETS1_sig_all.hg19_38.txt",col_names=F)
TOP_LD<-read_tsv("RA_GWAS_top_5_dash_LD06.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)
SIG<-read_tsv("RA_GWAS_significant.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)
SIG_LD<-read_tsv("RA_GWAS_sig_plus_LD06.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

colnames(SNP)<-c("hg19","hg38")

OL_WITH_RA_CTL_FIBRO<-vector()
OL_WITH_RA_TNF_FIBRO<-vector()
OL_WITH_RA_CD4<-vector()
OL_WITH_RA_Mono<-vector()
OL_WITH_RA_CTL_FIBRO_2<-vector()
OL_WITH_RA_TNF_FIBRO_2<-vector()
OL_WITH_RA_NK<-vector()
OL_WITH_OA_CD4<-vector()
OL_WITH_OA_B<-vector()
OL_WITH_OA_TNF_FIBRO<-vector()
OL_WITH_OA_CTL_FIBRO<-vector()
OL_WITH_OA_NK<-vector()
OL_WITH_CCRE_PLS<-vector()
OL_WITH_CCRE_pELS<-vector()
OL_WITH_CCRE_dELS<-vector()
OL_WITH_CCRE_DNase_H3K4me3<-vector()
OL_WITH_CCRE_CTCF_only<-vector()

OL_WITH_RA_CTL_FIBRO<-read_tsv("RA_GWAS_sig_plus_LD06.ChIP_CTL_1_1e-2.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

OL_WITH_RA_TNF_FIBRO<-read_tsv("RA_GWAS_sig_plus_LD06.ChIP_72h_1_1e-2.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

OL_WITH_CCRE_PLS<-read_tsv("RA_GWAS_sig_plus_LD06.cCRE_PLS.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

OL_WITH_CCRE_pELS<-read_tsv("RA_GWAS_sig_plus_LD06.cCRE_pELS.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

OL_WITH_CCRE_dELS<-read_tsv("RA_GWAS_sig_plus_LD06.cCRE_dELS.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

OL_WITH_CCRE_DNase_H3K4me3<-read_tsv("RA_GWAS_sig_plus_LD06.cCRE_DNase_H3K4me3.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

OL_WITH_CCRE_CTCF_only<-read_tsv("RA_GWAS_sig_plus_LD06.cCRE_CTCF_only.OL.bed",col_names=F)%>%
mutate(NAME=str_c(X1,":",X3))%>%
pull(NAME)

data<-SNP %>% mutate(TOP=ifelse(hg19=="chr11:128496952",1,0))%>%
mutate(TOP_LD_0.6=ifelse(hg19%in%TOP_LD,1,0))%>%
mutate(SIG_IN_GWAS=ifelse(hg19%in%SIG,1,0))%>%
mutate(SIG_IN_GWAS_LD_0.6=ifelse(hg19%in%SIG_LD,1,0))%>%
mutate(RA_CTL_FIBRO_H3K27ac=ifelse(hg19%in%OL_WITH_RA_CTL_FIBRO,2,0))%>%
mutate(RA_TNF_FIBRO_H3K27ac=ifelse(hg19%in%OL_WITH_RA_TNF_FIBRO,2,0))%>%
mutate(cCRE_PLS=ifelse(hg38%in%OL_WITH_CCRE_PLS,3,0))%>%
mutate(cCRE_pELS=ifelse(hg38%in%OL_WITH_CCRE_pELS,3,0))%>%
mutate(cCRE_dELS=ifelse(hg38%in%OL_WITH_CCRE_dELS,3,0))%>%
mutate(cCRE_DNase_H3K4me3=ifelse(hg38%in%OL_WITH_CCRE_DNase_H3K4me3,3,0))%>%
mutate(cCRE_CTCF_only=ifelse(hg38%in%OL_WITH_CCRE_CTCF_only,3,0))

LIST<-colnames(data)
LIST<-LIST[3:length(LIST)]

data$hg19<-factor(data$hg19,levels=unique(data$hg19))
data$hg38<-factor(data$hg38,levels=unique(data$hg38))

data<-data %>%
filter(SIG_IN_GWAS==1)%>%
dplyr::select(all_of(c("hg19","hg38","RA_CTL_FIBRO_H3K27ac","RA_TNF_FIBRO_H3K27ac","cCRE_PLS",
"cCRE_pELS","cCRE_dELS","cCRE_DNase_H3K4me3",
"cCRE_CTCF_only")))%>%
pivot_longer(names_to="Category",values_to="value",-c("hg19","hg38"))

data$Category<-factor(data$Category,levels=rev(LIST))
data$col<-factor(data$value)
data<-mutate(data,value=ifelse(value>0,1,0))

p<-ggplot(data,aes(x=hg19,y=Category,size=value,colour=col))+
geom_point()+
theme_classic()+
labs(x="",y="") +
theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5,size=12),
    axis.text.y = element_text(size=12))+
scale_colour_manual(values=c("grey","red","blue"))
ggsave("hg19_GWAS_summary_fig.pdf",p,width=12,height=4)
ggsave("hg19_GWAS_summary_fig.png",p,width=12,height=4)


VAR<-read.table("5e-8_all.txt",header=T)
VAR$SNP<-factor(VAR$SNP,levels=VAR$SNP)
p<-ggplot(VAR,aes(x=SNP,y=-log10(P)))+
geom_point(shape=23,size=5,colour="#000033",fill="#FF0099")+
theme_classic()+
labs(x="",y="") +
ylim(0,11)+geom_hline(yintercept=-log10(5e-8),linetype="dashed",colour="#000000")+
theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5,size=12),
    axis.text.y = element_text(size=12))
ggsave("hg19_GWAS_P.pdf",p,width=12,height=4)
ggsave("hg19_GWAS_P.png",p,width=12,height=4)
