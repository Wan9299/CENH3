
## genome asm
hifiasm -l 2 -n 4 -o asm \
--h1 hic.r1.fq.gz \
--h2 hic.r2.fq.gz \
hifi.fasta.gz

run-asm-pipeline.sh \
--editor-coarse-resolution 100000 \
--editor-fine-resolution 1000 \
contigs.fasta \
merged_nodups.txt

## tandem repeat
Rscript TRASH.R \
-p 24 -m 10000 -i 200 \
-f genome.fa \
-o genome.trash

## TE anno
EDTA.pl --genome genome.fasta \
--species others \
--u 2.35e-8 \
--anno 1 \
--sensitive 1 \
--evaluate 1

## gene anno

cd-hit -i pep.fa -o pep.cdhit90.fa -c 0.9 -n 5 -M 16000 -d 0 -T 8

miniprot -t 24 --gff -K 20M genome.fa pep.cdhit90.fa > miniprot.gff

fastp -w 16 \
-i rna.r1.fq.gz \
-I rna.r2.fq.gz \
-o out.good.r1.fq.gz \
-O out.good.r2.fq.gz \

hisat2 --new-summary --dta --max-intronlen 20000 \
-x softmask.genome.index \
-1 out.good.r1.fq.gz \
-2 out.good.r2.fq.gz | \
samtools sort -@ 4 > out.bam && \
samtools index -@ 20 out.bam out.bam.bai

braker.pl --gff3 \
--genome=softmask.genome.fa \
--species=species \
--prot_seq=pep.cdhit90.fa \
--bam=all.bam

gtf_genome_to_cdna_fasta.pl \
braker/GeneMark-ETP/rnaseq/stringtie/transcripts_all.sorted.gff \
softmask.genome.fa > transcripts.all.fasta

gtf_to_alignment_gff3.pl transcripts_all.sorted.gff > transcripts.all.gff3

TransDecoder.LongOrfs -t transcripts.all.fasta

TransDecoder.Predict -t transcripts.all.fasta

cdna_alignment_orf_to_genome_orf.pl \
transcripts.all.fasta.transdecoder.gff3 \
transcripts.all.gff3 \
transcripts.all.fasta > transcripts.all.transdecoder.gff3

EVidenceModeler \
--weights list \
--sample_id sample \
--genome genome.fa \
--gene_predictions gene_predictions.gff \
--protein_alignments miniprot.match.gff \
--transcript_alignments trans_evm.gff \
--segmentSize 5000000 \
--overlapSize 5000

Launch_PASA_pipeline.pl -c config -C -R -g genome.fa -t transcripts.all.fasta --ALIGNERS gmap

Load_Current_Gene_Annotations.dbi -c config -g genome.fa -P evm.gff

Launch_PASA_pipeline.pl -c config -A -g genome.fa -t transcripts.all.fasta -L --annots_gff evm.gff

## evaluation
busco -m genome -i genome.fa -l embryophyta_odb10
busco -m proteins -i Protein.pep -l embryophyta_odb10