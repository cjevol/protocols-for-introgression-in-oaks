########################################################################################################################
# Part 1 ABBA-BABA test 
# This analysis was performed using the python script 'ABBABABAwindows.py' ,'calculate_abba_baba.r' and Fst values for the same windows were calculated using the script 'popgenWindows.py'
# The script is sourced at the following URL
# https://github.com/simonhmartin/genomics_general
# https://github.com/palc/tutorials-1
########################################################################################################################

## The D-statistic was quantified at whole-genome level
# step1 vcf file to genotype file 
tabix Quercus_all267_allchr_flt2SNP.vcf.gz
python3 parseVCFs.py -i Quercus_all267_allchr_flt2SNP.vcf.gz --threads 70 | bgzip -c > Quercus_all267_allchr_flt2SNP.geno2.gz
# step2 Genome wide allele frequencies
python3 freq.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P LFS -p TBS -p TB -p Outgroup  --popsFile all10pop_8qa_3qu.txt  --threads 70 --target derived -o Quercus_all267_allchr_.LFS_TBS_TB.Fre.tsv.gz
# step3 Zscore D fd caculate
Rscript calculate_abba_baba.r Quercus_all267_allchr_.LFS_TBS_TB.Fre.tsv.gz genome_LFS_TBS_TB_abba_baba.txt LFS TBS TB Outgroup chr_lengths.txt 

### To identify specific genomic regions subjected to introgression by nonoverlapping windows
## step1 Zscore D fd fdm caculate
#100kb
python3 ABBABABAwindows.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P1 LFS -P2 TBS -P3 TB -O Outgroup --popsFile all10pop_8qa_3qu.txt -w 100000 -m 200 --T 70 -f phased --writeFailedWindows -o Quercus_all267_allchr.ABBABABA_LFS_TBS_TB_out.w100km200.csv.gz 
#50kb
python3 ABBABABAwindows.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P1 LFS -P2 TBS -P3 TB -O Outgroup --popsFile all10pop_8qa_3qu.txt -w 50000 -m 100 --T 70 -f phased --writeFailedWindows -o Quercus_all267_allchr.ABBABABA_LFS_TBS_TB_out.w50km100.csv.gz 
#10kb
python3 ABBABABAwindows.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P1 LFS -P2 TBS -P3 TB -O Outgroup --popsFile all10pop_8qa_3qu.txt -w 10000 -m 20 --T 70 -f phased --writeFailedWindows -o Quercus_all267_allchr.ABBABABA_LFS_TBS_TB_out.w10km20.csv.gz 
## step2 Fst caculate
#100kb
python3 popgenWindows.py -w 100000 -p LFS -p TBS -p TB -o Outgroup --popsFile all10pop_8qa_3qu.txt -g Quercus_all267_allchr_flt2SNP.geno2.gz -f phased -T 70 --analysis {popFreq,popDist,popPairDist} -o Quercus_allchr_LFS_TBS_TB_allstatistic_w100k.csv.gz
#50kb
python3 popgenWindows.py -w 50000 -p LFS -p TBS -p TB -o Outgroup --popsFile all10pop_8qa_3qu.txt -g Quercus_all267_allchr_flt2SNP.geno2.gz -f phased -T 70 --analysis {popFreq,popDist,popPairDist} -o Quercus_allchr_LFS_TBS_TB_allstatistic_w50k.csv.gz
#10kb
python3 popgenWindows.py -w 10000 -p LFS -p TBS -p TB -o Outgroup --popsFile all10pop_8qa_3qu.txt -g Quercus_all267_allchr_flt2SNP.geno2.gz -f phased -T 70 --analysis {popFreq,popDist,popPairDist} -o Quercus_allchr_LFS_TBS_TB_allstatistic_w10k.csv.gz


########################################################################################################################
# Part 2 Mantel and partial Mantel regression tests by R package 'vegan'
########################################################################################################################

### The correlation between shared introgression and environmental/geographic variables
library(vegan)
library(dplyr)
geodistancematrix<- read.table("D:/mantel/equel2indsample/matrix/geo_distance10pop_matrix.txt",header = TRUE) %>% as.matrix
all91envir<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/all10pop_91environ_PCA_matrix.txt",header = TRUE) %>% as.matrix
Fst <- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/reorder_all10qa_fst_scalar_matrix.txt",header = TRUE) %>% as.matrix
PSIG<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/fd/fdpgi_all10qv_10k.txt",header = TRUE) %>% as.matrix
longitude<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/longitude_matrix.txt",header = TRUE) %>% as.matrix
latitude<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/latitude_matrix.txt",header = TRUE) %>% as.matrix
pre<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/all10pop_pre_PCA_matrix.txt",header = TRUE) %>% as.matrix
tem<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/all10pop_tem_PCA_matrix.txt",header = TRUE) %>% as.matrix
srad<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/all10pop_srad_PCA_matrix.txt",header = TRUE) %>% as.matrix
wind<- read.table("D:/mantel/equel2indsample/refilter_fd/matrix/PCA/all10pop_wind_PCA_matrix.txt",header = TRUE) %>% as.matrix
# Mantel test
mantel(all91envir, PSIG)
mantel(geodistancematrix, PSIG)
mantel(longitude, PSIG)
mantel(latitude, PSIG)
mantel(pre, PSIG)
mantel(tem, PSIG)
mantel(srad, PSIG)
mantel(wind, PSIG)
mantel(all91envir, geodistancematrix)
# partial Mantel test
mantel.partial( geodistancematrix,PSIG, Fst)
mantel.partial( geodistancematrix,PSIG, all91envir)
mantel.partial(environ8matrix, PSIG, geodistancematrix)
mantel.partial(environ8matrix, PSIG, Fst)
mantel.partial(all91envir, PSIG, geodistancematrix)
mantel.partial(all91envir, PSIG, Fst)
mantel.partial(longitude, PSIG, geodistancematrix)
mantel.partial(longitude, PSIG, Fst)
mantel.partial(longitude, PSIG, all91envir)
mantel.partial(latitude, PSIG, geodistancematrix)
mantel.partial(latitude, PSIG, Fst)
mantel.partial(pre, PSIG, geodistancematrix)
mantel.partial(pre, PSIG, Fst)
mantel.partial(tem, PSIG, geodistancematrix)
mantel.partial(tem, PSIG, Fst)
mantel.partial(srad, PSIG, geodistancematrix)
mantel.partial(srad, PSIG, Fst)
mantel.partial(wind, PSIG, geodistancematrix)
mantel.partial(wind, PSIG, Fst)

### The correlation between genetic distance (FST) and environmental/geographic variables
mantel(geodistancematrix, Fst)
mantel(all91envir, Fst)
mantel(longitude, Fst)
mantel(latitude, Fst)
mantel(pre, Fst)
mantel(tem, Fst)
mantel(srad, Fst)
mantel(wind, Fst)
# partial Mantel test
mantel.partial(geodistancematrix,Fst, all91envir)
mantel.partial(all91envir, Fst, geodistancematrix)
mantel.partial(longitude, Fst, geodistancematrix)
mantel.partial(latitude, Fst, geodistancematrix)
mantel.partial(pre, Fst, geodistancematrix)
mantel.partial(tem, Fst, geodistancematrix)
mantel.partial(srad, Fst, geodistancematrix)
mantel.partial(wind, Fst, geodistancematrix)


########################################################################################################################
# Part 3 RNA data processing,mapping to refrence genome using GSNAP,read count using featureCounts
########################################################################################################################

### mapping to refrence genome
# step1 build index
gmap_build -D /usr/local/share -d Qa /data/frr/quercus/genome/genome.fasta
cat Quce_anno.gtf | gtf_splicesites > /data/frr/Quercus/RNA_seqcallsnp/0_mapping/splice/Qa.splicesites
cat Qa.splicesites | iit_store -o Qa.gtf.index
# step2 run gsnap
for i in {1..3};do gsnap --gunzip -N 1 -s Qa.gtf.index.iit -A sam -d Qa -t 60 --split-output TBS${i} /data/frr/RNAseq_0804/0_clean_data/TBS${i}_BRRL210002171-1A_1.fq.gz /data/frr/RNAseq_0804/0_clean_data/TBS${i}_BRRL210002171-1A_2.fq.gz

### read count
# step1 sort produce bam file
for i in {1..3};do less /data/frr/Quercus/RNA_seqcallsnp/0_mapping/TBS${i}.concordant_uniq |samtools view -@ 20 -b - | samtools sort -@ 20 - > /data/frr/Quercus/RNA_seqcallsnp/0_mapping/concordant_bam/TBS${i}_srt.bam;done
# step2 run featureCounts
/data/frr/RNAseq_0804/0_software/subread-2.0.1-Linux-x86_64/bin/featureCounts -T 30 -p -t exon -g gene_id -a Quce_anno.gtf -o allfeaturecount.txt /data/frr/Quercus/RNA_seqcallsnp/0_mapping/concordant_bam/*.bam

### logFC caculate by R package 'edgeR'
library(limma)
library("edgeR")
library(dplyr)
countmatrix_keep <- read.delim('D:/quercus/resequence/RNA/gsnap_count/allfeaturecount_idmodified.txt', row.names=1, stringsAsFactors=FALSE, sep='\t',header=T);
HCS_KM <- select(countmatrix_keep,1,2,3,4,5,6)
#DGEList
KM_group <- factor (c(rep('HCS',3),rep('KM',3)))
KM_DEG_list_data <- DGEList(counts = HCS_KM, group = KM_group)
head(KM_DEG_list_data$counts)
KM_DEG_list_data$samples
#filter
keep <- filterByExpr(KM_DEG_list_data)
KM_DEG_list_data_filtered <- KM_DEG_list_data[keep, , keep.lib.sizes=FALSE]
head(KM_DEG_list_data_filtered$counts)
KM_DEG_list_data_filtered$samples
#normalization
KM_DEG_list_data_norm <- calcNormFactors(KM_DEG_list_data_filtered,method = 'TMM')
png("D:/quercus/resequence/RNA/gsnap_count/HCS_KM_MDS.png",width=3500,height=2150,res=20*20)
plotMDS(KM_DEG_list_data_norm,col=as.numeric(KM_DEG_list_data_norm$samples$group))
legend("bottomright",as.character(unique(KM_DEG_list_data_norm$samples$group)), col=1:4, pch=20)
dev.off()
design <- model.matrix(~0+KM_group)
disper <- estimateDisp(KM_DEG_list_data_norm,design)
png("D:/quercus/resequence/RNA/gsnap_count/HCS_KM_BCV.png",width=3500,height=2150,res=20*20)
plotBCV(disper)
dev.off()
fit <- glmFit(disper,design) 
KM_HCS <- glmLRT(fit, contrast=c(1,-1)) 
deKM_HCS <- decideTestsDGE(KM_HCS, adjust.method="BH", p.value = 0.05)
summary(deKM_HCS)
KM_HCS_result <- topTags(KM_HCS, n = Inf, adjust.method = 'BH', sort.by = 'PValue') %>% as.data.frame()
KM_HCS_result <- data.frame(geneid=rownames(KM_HCS_result),KM_HCS_result)
KM_HCS_result[which(KM_HCS_result$logFC >= 1 & KM_HCS_result$FDR <=0.05),'sig'] <- 'Up'
KM_HCS_result[which(KM_HCS_result$logFC <= -1 & KM_HCS_result$FDR <=0.05),'sig'] <- 'Down'
KM_HCS_result[which( KM_HCS_result$FDR >0.05),'sig'] <- 'NotSig'
write.table(KM_HCS_result,file="D:/quercus/resequence/RNA/gsnap_count/edgeR_expression/HCS_KM.txt",row.name=F,quote=FALSE,sep="\t")

