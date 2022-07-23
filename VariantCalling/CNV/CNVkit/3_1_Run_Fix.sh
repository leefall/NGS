#!/usr/bin/bash
touch /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/CNVkit/MT_1
for sFile in /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/CNVkit/Target_001/*.cnn;
do
#echo $sFile
cnvkit.py fix $sFile /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/CNVkit/MT_1 /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/CNVkit/Target_001/ref-tas.cnn --no-edge
done

