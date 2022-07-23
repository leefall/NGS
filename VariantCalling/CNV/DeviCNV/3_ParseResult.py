#!/usr/bin/env python
import glob
import sys





def ParseMQ20Result(fout,sFile):
	
	
	fp=open(sFile)
	
	for sLine in fp.readlines():
		if sLine[0]=="#":
			pass
		else:
			sLine=sLine.strip()
			t=sLine.split("\t")
			(sScore,sCNV,sID,sIDGene)=(t[0],t[2],t[5],t[1])
			
			if int(sScore)>=5:
				
				
				
				if sCNV=="del":
					sOutCNV="Del"
				elif sCNV=="dup":
					sOutCNV="Amp"
				else:
					print("Exception Case!!!!!!!!!!!!!")
					sys.exit()
				
				sGene=sIDGene.split("_")[1]
				
				
				
				fout.write("{0}\n".format("\t".join([sID,sGene,sOutCNV])))
				
			
			
			
			
	
	
	
	
	
	

if __name__=="__main__":
	#dResultDict=dict()
	
	fout=open("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/Result/Total_Result_MQ20_DeviCNV.txt","w")
	
	
	
	
	ParseMQ20Result(fout,"/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/DeviCNV/Target001/08.CNV/Target001.CNVWithScore.MQ20.dupdelTh1.3_0.7.txt")
	ParseMQ20Result(fout,"/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/DeviCNV/Target002/08.CNV/Target002.CNVWithScore.MQ20.dupdelTh1.3_0.7.txt")
	ParseMQ20Result(fout,"/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/DeviCNV/Target003/08.CNV/Target003.CNVWithScore.MQ20.dupdelTh1.3_0.7.txt")
	ParseMQ20Result(fout,"/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/DeviCNV/Target004/08.CNV/Target004.CNVWithScore.MQ20.dupdelTh1.3_0.7.txt")
	ParseMQ20Result(fout,"/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/DeviCNV/Target005/08.CNV/Target005.CNVWithScore.MQ20.dupdelTh1.3_0.7.txt")
	ParseMQ20Result(fout,"/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/CNV/DeviCNV/Target010/08.CNV/Target010.CNVWithScore.MQ20.dupdelTh1.3_0.7.txt")
	
	
	