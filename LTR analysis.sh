

## calculate intact LTR insert time
cd genome.mod.EDTA.raw/LTR
cat genome.mod.pass.list | \
awk 'BEGIN { FS=OFS="\t" } !/^#/ {print $1, $(NF-3), $NF}' | \
sed -e 's/:/\t/g' -e 's/\.\./\t/g' | \
awk 'BEGIN { FS=OFS="\t" } \
{ if ($4 == "-") { temp = $2; $2 = $3; $3 = temp } print }' | \
cut -f 1,2,3,5 | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4 -d 0 -o max \
> genome.LTR.time.bed

awk -F'\t' '{OFS="\t"; $4=$4/1000000; print $0}' \
genome.LTR.time.bed > genome.LTR.time.mya.bed

ucsc/bedGraphToBigWig \
genome.LTR.time.mya.bed \
genome.chr.len.list \
genome.LTR.time.mya.bw



## intact LTR tree and class
awk '$3~"LTR_retrotransposon"' \
genome.mod.EDTA.intact.gff3 | \
grep "Method=structural" > genome.intact.LTR.gff3

awk -F';' '{print $1}' genome.intact.LTR.gff3 | \
sed 's/ID=//' | awk '$7!="?"' | \
awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' > genome.intact.LTR.bed

awk -F';' '{print $1}' genome.intact.LTR.gff3 | \
sed 's/ID=//' | awk '$7=="?"' | \
awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t+"}' >> genome.intact.LTR.bed

bedtools getfasta -fi "$g" \
-bed genome.intact.LTR.bed -fo - -name -s | \
awk -F'::' '{print $1}' | seqtk seq -l 60 > genome.intact.LTR.fa

TEsorter genome.intact.LTR.fa -db rexdb-plant -p 20 > TEsorter.log  2> TEsorter.err.log

grep -P "\-RT\t" genome.intact.LTR.fa.rexdb-plant.dom.tsv > genome.intact.LTR-RT.dom.tsv
grep -P "Ty1-RT\t" genome.intact.LTR-RT.dom.tsv > genome.intact.Ty1-RT.dom.tsv
grep -P "Ty3-RT\t" genome.intact.LTR-RT.dom.tsv > genome.intact.Ty3-RT.dom.tsv

for i in LTR Ty1 Ty3 ; do
module load annotation
cut -f 1 genome.intact."$i"-RT.dom.tsv | \
seqtk subseq genome.intact.LTR.fa.rexdb-plant.dom.faa - \
> genome.intact."$i"-RT.raw.fa

## 简化ID
awk -F':' '{print $1}' genome.intact."$i"-RT.raw.fa | \
sed 's/\Class_I\/LTR\///' | \
sed -e "s/>/>genome-/" -e 's/|/_/g' -e 's/\//_/g' -e 's/\*//g' > genome.intact."$i"-RT.fa
mafft --auto genome.intact."$i"-RT.fa > genome.intact."$i"-RT.aln.fa
done



## sequence similarity identification
snakemake --cores 48 make_figures \
-s ./StainedGlass/workflow/Snakefile \
--config sample=genome \
fasta=seq.fa \
cooler_window=100 window=10000 mm_f=10000 tempdir=temp
