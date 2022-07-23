#!/usr/bin/env python
import glob, os

#lBamControlFile=glob.glob("/mnt/QNAP/leefall2/Package/Cancer_core_Control/*.bam")
#
#fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Panelcn_MOPS/Control_File_list.txt","w")
#
#for sFile in lBamControlFile:
#	fout.write("{0}\n".format(sFile))
#
lCase1File=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_001/*.bam")
lCase2File=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_002/*.bam")

lCase3File=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_003/*.bam")
lCase4File=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_004/*.bam")
lCase5File=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_005/*.bam")
lCase10File=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_010/*.bam")

lCaseFiles=lCase1File+lCase2File+lCase3File+lCase4File+lCase5File+lCase10File

fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Panelcn_MOPS/Case_File_list.txt","w")

for sFile in lCaseFiles:
	fout.write("{0}\n".format(sFile))


