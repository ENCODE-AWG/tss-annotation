#!/bin/bash

accession=$1
key=$2
genomeRef=~/Lab/Reference/Human/hg38/chromInfo.txt

if [ ! -f $accession.bam ]
then
    wget https://www.encodeproject.org/files/$accession/@@download/$accession.bam
fi

echo "Sorting and indexing bam ..."
samtools sort $accession.bam -o tmp.bam
samtools index tmp.bam

echo "Seperating strands ..."
samtools view -b -F 16 tmp.bam -o pos.bam
samtools view -b -f 16 tmp.bam -o neg.bam

echo "Converting to bigWigs ..."
~/bin/bamToBigWig pos.bam ~/Lab/Reference/Human/hg38/chromInfo.txt $key.full.plus.bigWig
~/bin/bamToBigWig neg.bam ~/Lab/Reference/Human/hg38/chromInfo.txt $key.full.minus.bigWig

###Extracting 5'ends

echo "Extracting 5'ends ..."
bedtools bamtobed -i pos.bam | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_'$key'" \
    "\t" $5 "\t" $6}' > p.bed
bedtools bamtobed -i neg.bam | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_'$key'" \
    "\t" $5 "\t" $6}' > n.bed

awk '{print $1 "\t" $2 "\t" $2+1 "\t" $4 "\t" $5}' p.bed > tp5.bed
awk '{print $1 "\t" $3-1 "\t" $3 "\t" $4 "\t" $5}' n.bed > tn5.bed

echo "Converting to bedGraphs ..."
bedtools genomecov -bg -i tp5.bed -g $genomeRef | sort -k1,1 -k2,2n > p5.bg
bedtools genomecov -bg -i tn5.bed -g $genomeRef | sort -k1,1 -k2,2n > n5.bg

echo "Converting to bigWigs ..."
~/bin/bedGraphToBigWig p5.bg $genomeRef $key.5ends.plus.bigWig
~/bin/bedGraphToBigWig n5.bg $genomeRef $key.5ends.minus.bigWig

mv tp5.bed $key.5ends.plus.bed
mv tn5.bed $key.5ends.minus.bed
mv p.bed $key.full.plus.bed
mv n.bed $key.full.minus.bed

cp $key.*.bigWig /home/moorej3/Public-HTML/TSS-Annotations/Track-Hub/hg38/data
