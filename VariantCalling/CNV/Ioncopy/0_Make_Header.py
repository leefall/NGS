#!/usr/bin/env python
import glob


def Header(sFile):
	dIonHeader=dict()
	dInhouseHeader=dict()
	lAMPID=[]
	fp=open(sFile,"r")
	fInhouse=open("../Real_Header.txt","w")
	fIon=open("../Ioncopy_Header.txt","w")
	for sLine in fp.xreadlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sAmpID, sGeneID, nStartPosition, nDepth,nEndPosition)=(t[3],t[7].split("=")[-1],str(t[1]),int(t[-1]),str(t[2]))
		#print (sAmpID, sGeneID, nPosition, str(nDepth))
		sIonID=sGeneID+"_"+sAmpID+"_"+nStartPosition
		sInhouseID=str(nStartPosition)+"\t"+str(nEndPosition)+"\t"+sGeneID
		
		
		#fIon.write("{0}\n".format(sIonID))
		
		
		if not sIonID in dIonHeader.keys():
			dInhouseHeader[sIonID]=sInhouseID
		else:
			pass
		
		if not sIonID in lAMPID:
			lAMPID.append(sIonID)
		else:
			pass
	#fout=open("../Average/Average"+sFile,"w")
	for sAMPID in lAMPID:
		fIon.write("{0}\n".format(sAMPID))
		fInhouse.write("{0}\n".format(dInhouseHeader[sAMPID]))
		#fout.write("{0}\n".format(float(nTheDepth)/nLength))
		
		
		

if __name__=="__main__":
	sFilelist=glob.glob("IonXpress_*.txt")
	Header(sFilelist[0])
	
	
	
	


