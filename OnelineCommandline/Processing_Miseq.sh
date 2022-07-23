#!/usr/bin/bash

inputone=$1
inputtwo=$2
IssacDirectory=$3
Processed_Reference=$4
Manta=$5
Reference=$6
Strelka2=$7
IndexedTargetRegion=$8
TargetRegion=$9
FREEC=${10}
CurrentDirectory=$(pwd)


TemporDirectory=$(echo $inputone|awk -F"/" '{print $(NF-1)}')



DirectoryName=$(echo $TemporDirectory|awk -F"_" '{print $1}' ) 


mkdir $DirectoryName

ln -s $inputone $DirectoryName/lane1_read1.fastq.gz
ln -s $inputtwo $DirectoryName/lane1_read2.fastq.gz


$IssacDirectory/isaac-align -r $Processed_Reference -b $DirectoryName -m 100 --base-calls-format fastq-gz -o ./Result_$DirectoryName

mkdir manta_$DirectoryName

$Manta/configManta.py --bam ./Result_$DirectoryName/Projects/default/default/sorted.bam --referenceFasta $Reference --runDir manta_$DirectoryName

python manta_$DirectoryName/runWorkflow.py -m local

mkdir Strelka2_$DirectoryName


$Strelka2/configureStrelkaGermlineWorkflow.py --bam ./Result_$DirectoryName/Projects/default/default/sorted.bam --referenceFasta $Reference --targeted --callRegions=$IndexedTargetRegion --indelCandidates=manta_$DirectoryName/results/variants/candidateSmallIndels.vcf.gz  --runDir ./Strelka2_$DirectoryName

python ./Strelka2_$DirectoryName/runWorkflow.py -m local


mkdir FREEC_$DirectoryName


echo '''
[general]

chrLenFile = /storage/home/leefall2/tools/FREEC-11.0/cancer_panel.len
ploidy = 2
window = 0
breakPointThreshold = 1.2
chrFiles = /storage/home/reference/UCSC/hg19/Sequence/Chromosomes/
outputDir = FREEC_'''$DirectoryName'''
sambamba=/storage/home/leefall2/tools/sambamba_v0.6.6
SambambaThreads=20

maxThreads=20


[sample]

mateFile = '''$CurrentDirectory'''/Result_'''$DirectoryName'''/Projects/default/default/sorted.bam
inputFormat = bam

[control]



[target]

captureRegions = '''$TargetRegion''' '''   > config_FREEC.txt

 

echo $FREEC

$FREEC/freec -conf config_FREEC.txt





