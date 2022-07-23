#!/usr/bin/bash

#cnvkit.py batch -m amplicon -t target.split.bed /storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_001/*.bam

touch MT

for sFile in /storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_001/*.bam;
do
sPrefix1=$(echo $sFile|awk -F"/" '{print $(NF)}')

#echo $sPrefix1
sPrefix=$(echo $sPrefix1|awk -F"." '{print $1}')
#echo "cnvkit.py coverage $sFile target.split.bed -p 10 -o /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/CNVkit/Target_001/"$sPrefix.cnn
cnvkit.py coverage $sFile Gene_target.split.bed -p 10 -o /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/CNVkit/Target_001/$sPrefix.cnn

done


