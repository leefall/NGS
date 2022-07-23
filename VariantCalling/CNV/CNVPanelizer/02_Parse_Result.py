#!/usr/bin/env python
import glob
import os


fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Total_Result/Total_Result_CNV_Panelizer.txt","w")
Target_Gene=["CDK4","CDK6","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","MET","MYC"]
#lFilelist=["/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/CNVPanelizer/Batch2_with_Whole_Control_except_Case.txt"]
lFilelist=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/CNVPanelizer/*with_Whole_Control_except_Case.txt")
for sFile in lFilelist:
	fp=open(sFile)
	dIDCNV=dict()
	sFline=fp.readline()
	sFline=sFline.strip()
	t=sFline.split("\t")
	n=1
	for sCategory in t:
		if "ReliableStatus" in sCategory:
			sID=sCategory.split(".")[0]
			dIDCNV[n]=sID
			
		n+=1

	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		sGene=t[0]
		if sGene in Target_Gene:
			for i in range(0,len(t[1:])):
				if i in dIDCNV:
					sCell=t[i]
					if sCell=="Amplification":
						fout.write("{0}\t{1}\t{2}\n".format(dIDCNV[i],sGene,"Amp"))
					elif sCell=="Deletion":
						fout.write("{0}\t{1}\t{2}\n".format(dIDCNV[i],sGene,"Del"))
					else:
						pass

	fp.close()
		






