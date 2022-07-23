#!/usr/bin/env python
import glob
import os


sDirlist=["001","002","003","004","005","010"]
Target_Gene=["CDK4","CDK6","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","MET","MYC"]
#fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/VisCap/Total_Result.txt","w")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Total_Result/Total_Result_Viscap.txt","w")

for sDir in sDirlist:
	lFiles=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/VisCap/iter_Target_"+sDir+"_batch/iter_Target_"+sDir+"_batch_run1/*.xls")
	lFiles.remove("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/VisCap/iter_Target_"+sDir+"_batch/iter_Target_"+sDir+"_batch_run1/VisCap_run_info.xls")
	for sFile in lFiles:
		fp=open(sFile)
		fp.readline()
		for sLine in fp.readlines():
		#	print(t)
			t=sLine.split("\t")
			#print(t)
			
			if t[5] in Target_Gene:
				
				if t[2]=="Gain":
					fout.write("{0}\t{1}\t{2}\n".format(t[0],t[5],"Amp"))
				else:
					fout.write("{0}\t{1}\t{2}\n".format(t[0],t[5],"Del"))
		
		fp.close()

