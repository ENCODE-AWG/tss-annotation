
cellType=$1
mode=$2
boundary=$3

dataDir=/data/tss_annotations
scriptDir=~/Scripts/tss-annotation/cCRE-Intersections/
ccres=$dataDir/Chromatin-Datasets/cCREs/hg38/$cellType/$cellType-cCREs.bed
rankedList=$dataDir/Transcription-Datasets/$mode/hg38/$cellType/$cellType-$mode-*Peak-Quantifications.txt
peakList=$dataDir/Transcription-Datasets/$mode/hg38/All-Biosamples-$mode-*Peaks.bed


if [ $mode == "RAMPAGE" ]
then
    awk '{if (NR != 1 && $NF > 0) print $0}' $rankedList | sort -k3,3rg > tmp.ranked-list
    col=3
elif [ $mode == "CAGE" ]
then
    awk '{if (NR != 1 && $NF > 0) print $0}' $rankedList | sort -k2,2rg > tmp.ranked-list
    col=2
fi

if [ $boundary == "summit" ]
then
    peakList=$dataDir/Transcription-Datasets/$mode/hg38/All-Biosamples-$mode-Summits.bed
fi

grep -v "Low-DNase" $ccres > tmp.active-ccres
totalPeaks=$(wc -l tmp.ranked-list | awk '{print $1}')

rm -f log.*
for i in `seq 500 500 $totalPeaks`
do
    echo "Processing" $i "peaks..."
    head -n $i tmp.ranked-list | tail -n 500 > tmp.mini-list
    x=$(tail -n 1 tmp.mini-list | awk '{print $'$col'}')
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.mini-list $peakList > tmp.mini-bed
    bedtools intersect -u -a tmp.mini-bed -b tmp.active-ccres | wc -l | 
	awk '{print "'$i'" "\t" "'$x'" "\t" $1 "\t" $1/500*100}' >> log.tss
    bedtools intersect -u -a tmp.active-ccres -b tmp.mini-bed > tmp.intersect
    python $scriptDir/count-ccre-groups.py tmp.intersect | \
	awk '{print "'$i'" "\t" "'$x'" $0}' >> log.ccre 
done

mv log.tss $mode-$cellType-TSS-Overlap.txt
mv log.ccre $mode-$cellType-cCRE-Overlap.txt

Rscript $scriptDir/plot-ccre-intersections.R $cellType $mode

rm -f tmp.* log.*
