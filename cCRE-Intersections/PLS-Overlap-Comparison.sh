cellType=$1

dataDir=~/GitHub/tss-annotation/Annotations
scriptDir=~/GitHub/tss-annotation/cCRE-Intersections/
ccres=$dataDir/Chromatin-Datasets/cCREs/hg38/$cellType/$cellType-cCREs.bed

rampageBed=$dataDir/Transcription-Datasets/RAMPAGE/hg38/$cellType/$cellType-RAMPAGE-*Peaks.RPM-1.bed
cageBed=$dataDir/Transcription-Datasets/CAGE/hg38/$cellType/$cellType-CAGE-*Peaks.TPM-1.bed
pacbioBed=$dataDir/Transcription-Datasets/PacBio/hg38/$cellType/$cellType-PacBio-Sequel*-TSS-Counts.bed

grep "PLS" $ccres > tmp.pls

bedtools intersect -u -a tmp.pls -b $rampageBed > tmp.rampage
bedtools intersect -u -a tmp.pls -b $cageBed > tmp.cage
bedtools intersect -u -a tmp.pls -b $pacbioBed > tmp.pacbio

awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.rampage tmp.cage > tmp.RC
awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.rampage tmp.pacbio > tmp.RP
awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.cage tmp.pacbio > tmp.CP

rcpCount=$(awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.RC tmp.RP | wc -l | awk '{print $1}')

rcCount=$(wc -l tmp.RC | awk '{print $1-'$rcpCount'}')
rpCount=$(wc -l tmp.RP | awk '{print $1-'$rcpCount'}')
cpCount=$(wc -l tmp.CP | awk '{print $1-'$rcpCount'}')

rCount=$(wc -l tmp.rampage | awk '{print $1-'$rcpCount'-'$rcCount'-'$rpCount'}')
pCount=$(wc -l tmp.pacbio | awk '{print $1-'$rcpCount'-'$cpCount'-'$rpCount'}')
cCount=$(wc -l tmp.cage | awk '{print $1-'$rcpCount'-'$rcCount'-'$cpCount'}')

cat tmp.rampage tmp.cage tmp.pacbio | sort -u > tmp.overlap
noCount=$(awk 'FNR==NR {x[$0];next} !($0 in x)' tmp.overlap tmp.pls | wc -l | awk '{print $1}')

echo -e "All-three" "\t" "RAMPAGE-CAGE" "\t" "RAMPAGE-PacBio" "\t" "CAGE-PacBio" \
    "\t" "RAMPAGE-only" "\t" "CAGE-only" "\t" "PacBio-only" "\t" "None"
echo -e $rcpCount "\t" $rcCount "\t" $rpCount "\t" $cpCount "\t" $rCount \
    "\t" $cCount "\t" $pCount "\t" $noCount

bedtools intersect -u -a tmp.pls -b $rampageBed | awk '{print $0 "\t" "R"}' > tmp.bed
bedtools intersect -v -a tmp.pls -b $rampageBed | awk '{print $0 "\t" "-"}' >> tmp.bed
mv tmp.bed tmp.pls

bedtools intersect -u -a tmp.pls -b $cageBed | awk '{print $0 "C"}' > tmp.bed
bedtools intersect -v -a tmp.pls -b $cageBed | awk '{print $0 "-"}' >> tmp.bed
mv tmp.bed tmp.pls

bedtools intersect -u -a tmp.pls -b $pacbioBed | awk '{print $0 "P"}' > tmp.bed
bedtools intersect -v -a tmp.pls -b $pacbioBed | awk '{print $0 "-"}' >> tmp.bed
mv tmp.bed $cellType-PLS-Overlap-Groups.bed

rm tmp.*
