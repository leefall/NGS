#!/usr/bin/env python
import os
import glob
def Parse_Barcode_Dictionary():
	
	#lFilterlist=["SNU055","SNU041","SNU054","SNU053","SNU104,SNU101,SNU064,SNU140,SNU163,SNU162,SNU161,SNU294,SNU247,SNU281,SNU295,SNU283,SNU251,SNU248,SNU293,SNU246,SNU252,SNU243,SNU343,SNU359,SNU341,SNU344
	lFilterlist=[]
	fFilter=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Uniformity_Filter_ID.txt")
	
	
	for sLine in fFilter.readlines():
		sLine=sLine.strip()
		lFilterlist.append(sLine)
	
	
	dBarcode_Dict=dict()
	fp=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/Sample_Barcode.txt","r")
	for sLine in fp.readlines():
		sLine=sLine.strip()
		(sTarget,sBarcode,sSampleID)=sLine.split("\t")
		sTarget=sTarget.replace("_","")
		sTarget_Barcode=sTarget+"|"+sBarcode
		
		if "SNU" in sSampleID:
			if not sSampleID in lFilterlist:
				dBarcode_Dict[sTarget_Barcode]=sSampleID
	return(dBarcode_Dict)



def ChangedID(sFile,dIDdict):
	fp=open(sFile)
	fout=open("ChangedID/"+sFile,"w")
	sFile=sFile.split("/")[-1]
	sTargetInfo=sFile.split(".")[0]
	nOutlist=[0,1]
	sLine=fp.readline()
	sLine=sLine.strip()
	lLine=sLine.split("\t")
	nIndex=2
	lTemporOut=[]
	fout.write(lLine[0]+"\t"+lLine[1]+"\t")
	#lTemporLine=lLine[2:]
	for sCate in lLine[2:]:
		sCate=sCate.split("_")[1]
		sTargetID=sTargetInfo+"|"+sCate
		if sTargetID in dIDdict:
			nOutlist.append(nIndex)
			lTemporOut.append(dIDdict[sTargetID])
			#fout.write(dIDdict[sTargetID])
		else:
			print(sTargetID)
			pass
		
		nIndex+=1
	fout.write("{0}\n".format("\t".join(lTemporOut)))
	
	for sLine in fp.readlines():
		nIndex=0
		sLine=sLine.strip()
		lLine=sLine.split("\t")
		lTemporOut=[]
		for sCate in lLine:
			if nIndex in nOutlist:
				lTemporOut.append(sCate)
			else:
				pass
			nIndex+=1
		fout.write("{0}\n".format("\t".join(lTemporOut)))
		
		
		
	
	
	
	#print(nOutlist)
	#print(len(nOutlist))
	
	
	
	
	
	
	
	



if __name__=="__main__":
	
	
	dIDdict=Parse_Barcode_Dictionary()
	print(len(dIDdict))
	os.chdir("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ReadCount")
	lXlsfiles=glob.glob("*.xls")
	for sFile in lXlsfiles:
		ChangedID(sFile,dIDdict)
	#ChangedID("Target001.xls",dIDdict)
	
	
	


