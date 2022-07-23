#!/usr/bin/env python
sFile="/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ReadCount/ChangedID/Target001.xls"

#sFirstLine=sFile
dGeneCount=dict()
dAmpliconInformction=dict()
fSequence=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/Header_Caner_Panel_Seqeuncing_GC.bed")

for sLine in fSequence.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	sAmplicon=t[3]
	
	(sChr, sStart, sEnd, sTemporGene, sGC)=(t[0],t[1],t[2],t[7],t[8])
	#sTemporGene=sTemporGene.
	sGene=sTemporGene.split("=")[1]
	sGene=sGene.replace("_1","")
	sGene=sGene.replace("_2","")
	sGene=sGene.replace("_HS","-HS")
	if sGene in dGeneCount.keys():
		dGeneCount[sGene]+=1
		sGeneNumber=sGene+"_"+str(dGeneCount[sGene])
	else:
		dGeneCount[sGene]=1
		sGeneNumber=sGene+"_"+str(dGeneCount[sGene])
	
	
	
	dAmpliconInformction[sAmplicon]=[sChr,sStart,sEnd,sGeneNumber,sGC]
	#dAmplSeqeunce[t[0]]=t[1]
	
	
	
#fBed=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/Cancer_panel_sequencing_target.bed")

#fBed.readline()
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ExomeDepth/Target001.txt","w")
fout.write("seqnames\tstart\tend\twidth\tstrand\tGC\t")
fp=open(sFile)
sFirstLine=fp.readline()
sFirstLine=sFirstLine.strip()
lFirstLine=sFirstLine.split("\t")
lIDline=lFirstLine[2:]
fout.write("{0}".format("\t".join(lIDline)))


#fout.write
fout.write("\t")
fout.write("names\tchromosome\n")
#seqnames start   end width strand        GC Exome1 Exome2 Exome3 Exome4          names chromosome
#chr1 12012 12058    47      * 0.6170213      0      0      0      0 DDX11L10-201_1          1
#chr1 12181 12228    48      * 0.5000000      0      0      0      0 DDX11L10-201_2          1

dAmplSeqeunce=dict()

n=1

for sLine in fp.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	sAmplicon=t[1]
	#(sChr,sStart,sEnd,sAmplicon,sGeneID)=(t[0],t[1],t[2],t[3],t[7])
	#sGene=sGeneID.split("=")[1]
	#sGene=sGene.replace("_1","")
	#sGene=sGene.replace("_2","")
	
	[sChr,sStart,sEnd,sGeneNumber,sGC]=dAmpliconInformction[sAmplicon]
	
	noChr=sChr.replace("chr","")
	sWidth=str(int(sEnd)-int(sStart))
	
	fout.write("{0}\t".format("\t".join([sChr,sStart,sEnd,sWidth,"*",sGC])))
	
	fout.write("{0}\t".format("\t".join(t[2:])))
	fout.write(sGeneNumber)
	fout.write("\t")
	fout.write(noChr)
	fout.write("\n")
	
	
	#2	chr1	16254614	16254715	AMPL7154512776	SPEN-E11	AACCAAGATCGTACATATTATGAGAGTGTTCGAACTCCAGGCACTTATCCTGAGGATTCCAGGCGGGACTATCCAGCTCGAGGGAGAGAGTTTTATTCAGAA
	#fout.write("1\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(\
	#sChr,sStart,sEnd,sAmplicon,sGene,sSequence))
	
	
	
	





