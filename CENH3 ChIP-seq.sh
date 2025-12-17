## CENH3 ChIP-seq data process
fastp \
-i chip.r1.fq.gz \
-o chip.good.r1.fq.gz \
-I chip.r2.fq.gz \
-O chip.good.r2.fq.gz


bowtie2 \
-x genome.index \
-1 chip.good.r1.fq.gz \
-2 chip.good.r2.fq.gz | \
samtools sort > chip.bam && \
samtools index chip.bam chip.bam.bai && \
bamCoverage -p 16 -b chip.bam -o chip.bw

bamCompare -p 12 --normalizeUsing CPM --scaleFactorsMethod None \
-b1 cenh3.chip.bam \
-b2 input.bam \
-o log2ratio.bw
















