#!/usr/bin/env python
import glob
import sys
Target_Gene=["AKT1","BRAF","BRCA1","BRCA2","CDK6","CTNNB1","DDR2","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","FLT3","IDH1","IDH2","JAK1","JAK2","KIT","MAP2K1","MAP2K2","MET","MTOR","MYC","NPM1","PDGFRA","PIK3CA","RB1","STK11","TP53"]

fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/Result/Total_Result_CovCopCan_batch.txt","w")

lResult=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/Result/CovCopCan/Intersected/*.vcf")

dIDCNV=dict()


for sFile in lResult:
	fp=open(sFile)
	sID=sFile.split("/")[-1]
	sID=sID.split("_")[1]
	sID=sID.split(".")[0]
	dIDCNV[sID]=dict()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		sGene=t[-2]
		sGene=sGene.split("=")[1]
		if t[4]=="<DUP>":
			dIDCNV[sID][sGene]="Amp"
		elif t[4]=="<DEL>":
			dIDCNV[sID][sGene]="Del"
		else:
			print("Other CNV!!!!!!!!!!!!")
			print(sLine)
			sys.exit()
		
		
		
	

for sID in dIDCNV.keys():
	for sGene in dIDCNV[sID].keys():
		fout.write("{0}\t{1}\t{2}\n".format(sID,sGene,dIDCNV[sID][sGene]))





