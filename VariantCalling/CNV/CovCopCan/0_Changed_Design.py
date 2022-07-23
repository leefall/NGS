#!/usr/bin/env python
dAmplSeqeunce=dict()

fSequence=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/Sequence_per_CancerPanel_Deisgn.bed")

for sLine in fSequence.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	dAmplSeqeunce[t[0]]=t[1]
	
	
	
fBed=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/Cancer_panel_sequencing_target.bed")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CovCopCan_Design_CancerPanel.csv","w")
fBed.readline()
fout.write("#Pool\tChromosome\tStart\tEnd\tAmplicon\tGene\tSequence\n")
for sLine in fBed.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	(sChr,sStart,sEnd,sAmplicon,sGeneID)=(t[0],t[1],t[2],t[3],t[7])
	sGene=sGeneID.split("=")[1]
	sGene=sGene.replace("_1","")
	sGene=sGene.replace("_2","")
	
	sSequence=dAmplSeqeunce[sAmplicon]
	#2	chr1	16254614	16254715	AMPL7154512776	SPEN-E11	AACCAAGATCGTACATATTATGAGAGTGTTCGAACTCCAGGCACTTATCCTGAGGATTCCAGGCGGGACTATCCAGCTCGAGGGAGAGAGTTTTATTCAGAA
	fout.write("1\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(\
	sChr,sStart,sEnd,sAmplicon,sGene,sSequence))
	
	
	
	





