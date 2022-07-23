#!/usr/bin/env python
import glob, os

Target_Gene=["CDK4","CDK6","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","MET","MYC"]
lResult=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Panelcn_MOPS/Control_1/ResultTable/*.txt")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Total_Result/Total_Result_Control_1_Panelcn_MOPS.txt","w")
for sFile in lResult:
	fp=open(sFile)
	fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sSample, sGene, sLowQual, sCN)=(t[0],t[2],t[-2],t[-1])
		if sGene in Target_Gene:
			sSample=sSample.split(".")[0]
			nCN=sCN.split("N")[-1]
			nCN=int(nCN)
			if sLowQual=='':
				if ((nCN==1) or (nCN==0)):
					fout.write("{0}\t{1}\t{2}\n".format(sSample,sGene,"Del"))
				elif ((nCN==3) or (nCN==4)):
					fout.write("{0}\t{1}\t{2}\n".format(sSample,sGene,"Amp"))
				elif (nCN==2):
					pass
				else:
					print(sLine)
	
	

