#!/usr/bin/env python
import glob
#dFinal=set()
dTargetTCGA=dict()
sSetWhat=set()

#fp

#lPostrecalBam1=glob.glob("/mnt/alpha/leefall2/TCGA_HNSC/rawbam/Paired/Paired/Preprocess/postrecal_TCGA-*.bam")
lPostrecalBam=glob.glob("/mnt/Beta/leefall2/TCGA_LGG/postrecalbam/postrecal_TCGA-*bam")
#lPostrecalBam=lPostrecalBam1+lPostrecalBam2
#lPostrecalBam=lPostrecalBam1
fFinal=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/allmc3Barcode.txt")
fFinal.readline()



for sLine in fFinal.readlines():
	sLine=sLine.strip()
	sTCGAID=sLine.split("\t")[0]
	#sFinal.add(t[0])
	dTargetTCGA[sTCGAID+"-T"]=0
	dTargetTCGA[sTCGAID+"-N"]=0

for sPathFile in lPostrecalBam:
	sFile=sPathFile.split("/")[-1]
	sTCGAFile=sFile.split("_")[1]
	sTCGASample=sTCGAFile.split(".")[0]
	if sTCGASample in dTargetTCGA.keys():
		dTargetTCGA[sTCGASample]=1
	else:
		sSetWhat.add(sTCGASample)
	

sNotYet=set()

for sKey in dTargetTCGA.keys():
	if dTargetTCGA[sKey]==0:
		sNotYet.add(sKey)
print("Total Number of Sample for Process")
print(len(dTargetTCGA))
#print(dTargetTCGA)
print("Not Processed Yet")
print(len(sNotYet))
print(sNotYet)
print("What the")
print(len(sSetWhat))
print(list(sSetWhat)[0:10])
#print(sTNBC&sFinal)

