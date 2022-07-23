#!/usr/bin/env python
#dBedDict=dict()
dExonDict=dict()
dNMID=dict()
fExon=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/ANNOVAR_Tempora/Result_GenneAnno.txt.exonic_variant_function")
for sLine in fExon.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	
	#
	#MTOR:NM_004958:exon53:c.C7185A:p.N2395K
	sExonInfo=t[2]
	lExonInfo=sExonInfo.split(":")
	(sNMID,sExon)=(lExonInfo[1],lExonInfo[2])
	sAmplicon=t[-2]
	sGene=t[-1]
	
	dExonDict[sAmplicon]=sExon
	if sGene in dNMID.keys():
		pass
	else:
		dNMID[sGene]=sNMID
	
	
	
	
fExon.close()
	
fGenome=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/ANNOVAR_Tempora/Result_GenneAnno.txt.variant_function")

for sLine in fGenome.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	sAmplicon=t[-2]
	sClass=t[0]
	
	if sAmplicon in dExonDict.keys():
		pass
	else:
		dExonDict[sAmplicon]=sClass
	
	


fGenome.close()

fp=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CovCopCan_Design_CancerPanel.csv")
fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/DeviCNV_CancerPanel_probeInformation.txt","w")
##DeviCNV Example
#Amplicon_ID     Chr     Amplicon_Start  Amplicon_End    Insert_Start    Insert_End      Gene    Transcript      Exon    Pool
#MUT.1   6       49399400        49399520        49399400        49399520        MUT     NM_000255       Exon13  Pool1
#MUT.2   6       49399420        49399540        49399420        49399540        MUT     NM_000255       Exon13  Pool1
#MUT.3   6       49399440        49399560        49399440        49399560        MUT     NM_000255       Exon13  Pool1
#

fp.readline()
fout.write("Amplicon_ID\tChr\tAmplicon_Start\tAmplicon_End\tInsert_Start\tInsert_End\tGene\tTranscript\tExon\tPool\n")
for sLine in fp.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	(sAmpliconID,sChr,sStart,sEnd,sGene)=(t[4],t[1],t[2],t[3],t[5])
	
	sExon=dExonDict[sAmpliconID]
	sNMID=dNMID[sGene]
	
	
	fout.write("{0}\n".format("\t".join([sAmpliconID,sChr,sStart,sEnd,sStart,sEnd,sGene,sNMID,sExon,"Pool1"])))
	
	
	
	
	
	






