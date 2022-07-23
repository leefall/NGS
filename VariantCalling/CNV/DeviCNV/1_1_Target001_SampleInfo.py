#!/usr/bin/env python
import glob


dSexDict=dict()
lSexFiles=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/CNVkit/sex/SNU*.sex")

for sFile in lSexFiles:
	fp=open(sFile)
	sTemporID=sFile.split("/")[-1]
	sID=sTemporID.split(".")[0]
	fp.readline()
	sSecondLine=fp.readline()
	lSecondLine=sSecondLine.split("\t")
	dSexDict[sID]=lSecondLine[1]
	fp.close()
	
	
	


#print(dSexDict)


lTargetFiles=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_001/SNU*.bam")
fout=open("sampleInfo/Target001.txt","w")
fout.write("Sample\tSex\n")

for sFile in lTargetFiles:
	sTemporID=sFile.split("/")[-1]
	sID=sTemporID.split(".")[0]
	
	fout.write("{0}\t{1}\n".format(sID,dSexDict[sID]))
	
	





