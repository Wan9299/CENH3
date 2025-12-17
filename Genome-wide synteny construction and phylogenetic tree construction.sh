
### SVGAP pipeline

## 0. align
mkdir ./0.alignment
minimap2 -t 24 -x asm10 -c --cs=long ./genome/ref.genome.fa ./genome/query.genome.fa > ./0.alignment/refvsquery.paf


## 1. chainnet
mkdir ./1.chainnet
perl ./SVGAP/1_Convert2Axt.pl \
-ali minimap2 \
-input ./0.alignment \
-wk ./1.chainnet

## 2 synnet
perl /share/Pub/user/wan/sf/SVGAP/2_ChainNetSyn.pl \
--gd ./genome/ --ad 1.chainnet/Target_Ref/ --lst ./genome/g.lst --sing --syn ./2.synnet

## 3 filter
perl ./SVGAP/3_SynNetFilter.pl \
--synnet ./2.synnet/ --gd ./genome/ --chain 1.chainnet/Target_Ref

## 4 CleanSV
perl ./SVGAP/4_PairwiseSV.pl \
--syndir ./CleanSynNet --gd ./genome/ --outdir ./4.CleanSV

## 5 CombinedSV
perl ./SVGAP/5_Combined.pl \
--sv ./4.CleanSV --refname ref.genome --outdir ./5.CombinedSV

## 6 Genotyping Deletions (≥50 bp)
perl ./SVGAP/6_DELGenoptyping.pl \
--del ./5.CombinedSV/All.DELs.50bplarge.bed.combined.sorted.txt \
--singmaf ./filtered_sing/ --refseq ./genome/ref.genome --refsize ./genome/ref.genome.sizes \
--refname ref.genome --outdir ./6.VCFtemp

## 7 Genotyping Insertions (≥50 bp)
perl ./SVGAP/7_INSGenotyping.pl \
--ins ./5.CombinedSV/All.INTs.50bplarge.bed.combined.sorted.txt \
--query ./genome/ \
--sing ./filtered_sing/ \
--stretcher ./SVGAP/pub/Stretcher/ \
--refseq ./genome/ref.genome \
--refsize ./genome/ref.genome.sizes \
--refname ref.genome \
--outdir 7.INS_GT_temp

## 8 Genotyping InDels del
perl ./SVGAP/8_InDelGenotyping.pl \
--indel ./5.CombinedSV/All.DELs.50bpsmall.bed.combined.sorted.txt \
--refseq ./genome/ref.genome \
--refsize ./genome/ref.genome.sizes \
--refname ref.genome \
--query ./genome/ \
--indeltype del \
--singmaf ./filtered_sing/ \
--outdir 8.del.50bpsmall.temp


## 8 Genotyping InDels ins
perl ./SVGAP/8_InDelGenotyping.pl \
--indel ./5.CombinedSV/All.INTs.50bpsmall.bed.combined.sorted.txt \
--refseq ./genome/ref.genome \
--refsize ./genome/ref.genome.sizes \
--refname ref.genome \
--query ./genome/ \
--indeltype del \
--singmaf ./filtered_sing/ \
--outdir 8.ins.50bpsmall.temp

# 9 Genotyping SNPs (require big memory)
perl ./SVGAP/9_SNPgenotyping.pl \
--refseq ./genome/ref.genome \
--refsize ./genome/ref.genome.sizes \
--refname ref.genome \
--singmaf ./filtered_sing/ \
--outdir 9.SNPVCF




### call snp from cenh3 chip and input data of Sa01 - Sa07 using GATK pipeline

## 0 merge data
cat cenh3-chip.good.r1.fq.gz input.good.r1.fq.gz > chip.good.r1.fq.gz
cat cenh3-chip.good.r2.fq.gz input.good.r2.fq.gz > chip.good.r2.fq.gz

## 1 mapping
bwa-mem2 mem -t 10 -M -a \
-R "@RG\tID:S${n}\tSM:S${n}\tPL:ILLUMINA\tLB:lb" \
genome.index \
chip.good.r1.fq.gz \
chip.good.r2.fq.gz | \
samtools sort -@ 2 > chip.bam && \
samtools index -@ 12 chip.bam chip.bam.bai

## 2 remove duplicates
picard -Xmx200g -XX:ParallelGCThreads=4 MarkDuplicates \
I=chip.bam \
O=chip.picard.bam \
CREATE_INDEX=true \
REMOVE_DUPLICATES=true \
M="$g".picard/A"$n".map"$g".marked_dup_metrics.log

## 3 generate GVCF
gatk --java-options "-Xmx40g -Djava.io.tmpdir=tmp" HaplotypeCaller \
--intervals "$g"."$Chrs" \
--native-pair-hmm-threads 1 \
-R /CEPH/penglongwan/plant.genome/qiezi/"$g".fa \
--emit-ref-confidence GVCF \
-I "$g".picard/"$n".map"$g".picard.bam \
-O gatkHtC."$g"/"$Chrs"/"$n"."$Chrs".gvcf.gz

## 4 build DB and GenotypeGVCF
gatk --java-options "-Xmx100g -Djava.io.tmpdir=tmp" GenomicsDBImport \
--genomicsdb-workspace-path ./gatkDB \
--sample-name-map gvcf.map.list \
--reader-threads 4 \
--batch-size 50

gatk --java-options "-Xmx40g -Djava.io.tmpdir=tmp" GenotypeGVCFs \
-R genome.fa \
-V gendb://gatkDB \
-O GtG.vcf.gz


## 5 call SNP

gatk --java-options "-Xmx80g -Djava.io.tmpdir=tmp" SelectVariants \
-R genome.fa \
-V GtG.vcf.gz \
--select-type SNP \
-O raw.SNP.vcf.gz

gatk --java-options "-Xmx80g -Djava.io.tmpdir=tmp" VariantFiltration \
-R genome.fa \
-V raw.SNP.vcf.gz \
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name SNPfail \
-O filter.SNP.vcf.gz

gatk --java-options "-Xmx80g -Djava.io.tmpdir=tmp" SelectVariants \
-R genome.fa \
-V filter.SNP.vcf.gz \
--exclude-filtered \
-O filtered.SNP.vcf.gz










