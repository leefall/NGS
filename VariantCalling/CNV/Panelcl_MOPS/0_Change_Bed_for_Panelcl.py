#!/usr/bin/env python

fp=open("/storage/home/leefall2/mypro/DeviCNV/Target_Bed.bed")
fp.readline()

fout=open("//storage/home/leefall2/mypro/DeviCNV/Panelcl_Target_Bed.bed","w")

for sLine in fp.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	(sChr,nStart,nEnd,sAmplicon,sGeneID)=(t[0],t[1],t[2],t[3],t[-1])
	sGene=sGeneID.split("=")[-1]
	fout.write("{0}\n".format("\t".join([sChr, nStart, nEnd, sGene+"."+sAmplicon+"."+sChr+"."+nStart+"."+nEnd])))



