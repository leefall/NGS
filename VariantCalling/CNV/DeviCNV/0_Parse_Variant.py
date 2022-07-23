#!/usr/env/bin python
fp=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CovCopCan_Design_CancerPanel.csv")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/ForAnnotation_Desing.txt","w")
#Example input
#1       948921  948921  T       C       comments: rs15842, a SNP in 5' UTR of ISG15
#1       1404001 1404001 G       T       comments: rs149123833, a SNP in 3' UTR of ATAD3C
#1       5935162 5935162 A       T    

fp.readline()
for sLine in fp.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	nBase=int(t[3])-int(t[2])
	sRefVariant=t[-1][nBase//2]
	
	fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(t[1],str(int(t[2])+nBase//2),str(int(t[2])+nBase//2),sRefVariant,sRefVariant,t[-3],t[-2]))
	
	



