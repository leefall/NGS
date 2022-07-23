#!/usr/bin/env python
import glob

lResultsFiles=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ExomeDepth/Result/*.txt")

#lResultsFiles=lResultsFiles[0:1]

fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/Result/Total_Result_ExomeDepth_CaseOnly.txt","w")
dIDCNV=dict()


for sFile in lResultsFiles:
	fp=open(sFile)
	sID=sFile.split("/")[-1]
	sID=sID.split(".")[0]
	fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		
		(sCNV,sGenes)=(t[2],t[-1])
		
		
		if sCNV=="deletion":
			sOutCNV="Del"
		elif sCNV=="duplication":
			sOutCNV="Amp"
		
		
		sGenenumber=sGenes.split(",")[0]
		sGene=sGenenumber.split("_")[0]
		
		
		
		fout.write("{0}\n".format("\t".join([sID,sGene,sOutCNV])))
	
	






