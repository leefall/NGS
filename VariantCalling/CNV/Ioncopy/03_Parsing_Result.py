#!/usr/bin/env python
import glob

Target_Gene=["CDK4","CDK6","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","MET","MYC"]
flist=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Result/WithControl/*Bonferroni.txt")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Total_Result/Total_Result_Ioncopy_with_Control.txt","w")

for sFile in flist:
	fp=open(sFile)
	fp.readline()
	for sLine in fp.readlines():
		t=sLine.split("\t")
		if t[2] in Target_Gene:
			if t[-2]=="GAIN":
				fout.write("{0}\t{1}\t{2}\n".format(t[1],t[2],"Amp"))
			elif t[-2]=="LOSS":
				fout.write("{0}\t{1}\t{2}\n".format(t[1],t[2],"Del"))
		

fout.close()
fp.close()

flist=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Unfilter/Ioncopy/Result/BatchOnly/*Bonferroni.txt")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Total_Result/Total_Result_Ioncopy_BatchOnly.txt","w")

for sFile in flist:
	fp=open(sFile)
	fp.readline()
	for sLine in fp.readlines():
		t=sLine.split("\t")
		if t[2] in Target_Gene:
			if t[-2]=="GAIN":
				fout.write("{0}\t{1}\t{2}\n".format(t[1],t[2],"Amp"))
			elif t[-2]=="LOSS":
				fout.write("{0}\t{1}\t{2}\n".format(t[1],t[2],"Del"))
		


