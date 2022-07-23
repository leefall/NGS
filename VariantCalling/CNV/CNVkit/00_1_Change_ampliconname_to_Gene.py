#!/usr/bin/env python

dAmpliconGene=dict()
fp=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Cancer_Panel/Cancer_panel_sequencing_target.bed")

fp.readline()


for sLine in fp.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	(sAmplicon, sGeneName)=(t[3],t[-1])
	sGene=sGeneName.split("=")[-1]
	dAmpliconGene[sAmplicon]=sGene

fp.close()



fout=open("Gene_target.split.bed","w")
fp=open("target.split.bed","r")
for sLine in fp.readlines():
#	sLine=sLine.split("\t")
	sLine=sLine.strip()
	t=sLine.split("\t")
	if "," in t[-1]:
		sOneAmplicon=t[-1].split(",")[0]
		t[-1]=dAmpliconGene[sOneAmplicon]
	else:
		t[-1]=dAmpliconGene[t[-1]]
	print(t)
	fout.write("{0}\n".format("\t".join(t)))

