#!/usr/bin/env python
import os
import glob
def Parse_Control():
	lControlBam=glob.glob("/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Control_2/*.bam")
	
	lControllist=[]
	
	
	for sFile in lControlBam:
		sID=sFile.split("/")[-1]
		sID=sID.split(".")[0]
		lControllist.append(sID)
	
	
	#lFilterlist=["SNU055","SNU041","SNU054","SNU053","SNU104,SNU101,SNU064,SNU140,SNU163,SNU162,SNU161,SNU294,SNU247,SNU281,SNU295,SNU283,SNU251,SNU248,SNU293,SNU246,SNU252,SNU243,SNU343,SNU359,SNU341,SNU344
	dControlIDRead=dict()
	#fFilter=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Uniformity_Filter_ID.txt")
	
	for sID in lControllist:
		dControlIDRead[sID]=[]
		fp=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Control/"+sID+".txt")
		for sLine in fp.readlines():
			sLine=sLine.strip()
			sRead=sLine.split("\t")[4]
			dControlIDRead[sID].append(sRead)
		
		

	return(lControllist,dControlIDRead)



def AddControlID(sFile,lControllist,dControlIDRead):
	fp=open(sFile)
	fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ReadCount/ChangedID/WithControl/"+sFile,"w")
	sFile=sFile.split("/")[-1]
	sTargetInfo=sFile.split(".")[0]
	#nOutlist=[0,1]
	sLine=fp.readline()
	sLine=sLine.strip()
	lLine=sLine.split("\t")
	fout.write("{0}\t".format("\t".join(lLine)))
	fout.write("{0}\n".format("\t".join(lControllist)))
	
	
	nIndex=0
	#lTemporOut=[]
	#fout.write(lLine[0]+"\t"+lLine[1]+"\t")
	#lTemporLine=lLine[2:]
	for sLine in fp.readlines():
		lTemporOut=[]
		sLine=sLine.strip()
		fout.write(sLine)
		fout.write("\t")
		for sControlID in lControllist:
			lTemporOut.append(dControlIDRead[sControlID][nIndex])
		fout.write("{0}\n".format("\t".join(lTemporOut)))
		nIndex+=1
		
	
	
	
	
	#print(nOutlist)
	#print(len(nOutlist))
	
	
	
	
	
	
	
	



if __name__=="__main__":
	
	
	(lControllist,dControlIDRead)=Parse_Control()
	#print(len(dIDdict))
	#print(lControllist,dControlIDRead)
	os.chdir("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ReadCount/ChangedID")
	lXlsfiles=glob.glob("*.xls")
	for sFile in lXlsfiles:
		AddControlID(sFile,lControllist,dControlIDRead)
	#ChangedID("Target001.xls",dIDdict)
	
	
	


