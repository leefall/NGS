#!/usr/bin/bash


sTarget=Target_001
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' SNU*.txt > ../$sTarget\.txt
ls SNU*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt

sTarget=Target_002
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' SNU*.txt > ../$sTarget\.txt
ls SNU*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt


sTarget=Target_003
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' SNU*.txt > ../$sTarget\.txt
ls SNU*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt

sTarget=Target_004
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' SNU*.txt > ../$sTarget\.txt
ls SNU*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt

sTarget=Target_005
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' SNU*.txt > ../$sTarget\.txt
ls SNU*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt

sTarget=Target_010
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' SNU*.txt > ../$sTarget\.txt
ls SNU*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt

sTarget=Control
cd /storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Input/$sTarget

paste -d '\t' C*.txt > ../$sTarget\.txt
ls C*.txt|awk -F"." '{print $1}' > ../$sTarget\_ID.txt


