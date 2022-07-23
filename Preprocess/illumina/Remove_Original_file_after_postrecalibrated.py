#!/usr/bin/env python
import glob
import os
dTargetTCGA=dict()
dNotanalysis=dict()
sSetwhat=set()
def getSampleName():
	dUUID_Sample=dict()
	
	fp=open("Paired_TCGAID_table.txt","r")
	fp.readline()
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sUUID, sTCGAID)=(t[0],t[-1])
		dUUID_Sample[sUUID]=sTCGAID
		
	return dUUID_Sample

fFinal=open("/mnt/alpha/leefall2/TCGA_LGG/QC/Sample_for_Remove.txt")

for sLine in fFinal.readlines():
	sLine=sLine.strip()

	dTargetTCGA[sLine]=0
	

os.chdir("/gaia2/home/leefall2/mypro/TCGA_LGG/")
dUUID_Sample=getSampleName()
lDirlist=next(os.walk('.'))[1]
dIDdict=dict()

for sDir in lDirlist:
	if sDir in dUUID_Sample.keys():
		dSample=dUUID_Sample[sDir]
		#print(dSample)
		sTCGAbarcode=dSample[0:12]
		if dSample[13:15]=="01":
			sTissue="T"
		else:
			sTissue="N"
		sTCGASample=sTCGAbarcode+"-"+sTissue
		if sTCGASample in dTargetTCGA:
			#print(dSample)
			
			
			
			lBam=glob.glob(sDir+"/*.bam*")
			sCurrent=os.getcwd()
			for i in lBam:
				if not ".bai" in i:
					sBam=i
			#if dSample=="TCGA-A2-A3KC-01A":
#				print(lBam)
			try:
				sPathBam=sCurrent+"/"+sBam
				if sTCGASample in dIDdict.keys():
					dIDdict[sTCGASample]+=1
				else:
					dIDdict[sTCGASample]=1

				print("removing "+sPathBam)
				os.system("rm "+sPathBam)
			except NameError:
				pass
#		

