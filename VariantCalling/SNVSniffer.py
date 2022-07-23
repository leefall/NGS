#!/usr/bin/env python
import gzip
import os, subprocess, glob, time
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
dSampleDict=dict()
sCurrent=os.getcwd()
lCurrent=sCurrent.split("/")
sWorDir="/".join(lCurrent[0:-2])
sCancerType=lCurrent[-3]
#sInputDir="/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/"
sOutputDir="/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/SNVsniffer/"
lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/postrecalBam/Preprocess/postrecal_TCGA*.bam")
#lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
#lPostrecal=lpostrecalbam1+lpostrecalbam2
lPostrecal=lpostrecalbam1

#/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-UF-A7JT-N..bam
for sPathFile in lPostrecal:
	sFile=sPathFile.split("/")[-1]
	sID=sFile.split("_")[1]
	sID=sID.split(".")[0]
	if sID in dSampleDict.keys():
		print(sPathFile)
		#os.system("rm "+sPathFile)
	else:
		dSampleDict[sID]=sPathFile

try:
	os.mkdir(sOutputDir)
except:
	pass

logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)



def SamHeader(sTumorID):
	
	sTumorFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalFile=dSampleDict[sID+"-N"]
	
	
	#sNormalFile=sTumorFile.replace("T.","N.")
#	SampleName(sNormalFile)
#	fp=open('''/data/leefall2/TCGA/Preprocess/'''+sNormalFile+'''.samplename''')
#	sNormalID=fp.readline()
#	sNormalID=sNormalID.strip()
	sSample=sTumorFile.split("/")[-1]
	sSample=sSample.split("-T.")[0]
	subprocess.call('''

samtools view -S -H '''+sNormalFile+''' > '''+sOutputDir+'''/Normal_'''+sSample +'''
samtools view -S -H '''+sTumorFile+''' > '''+sOutputDir+'''/Tumor_'''+sSample+'''




#Remove temporal data
''',shell=True)






def SNVSnifferCall(sTumorID):

	sTumorFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalFile=dSampleDict[sID+"-N"]
#	fp=open('''/data/leefall2/TCGA/Preprocess/'''+sNormalFile+'''.samplename''')
#	sNormalID=fp.readline()
#	sNormalID=sNormalID.strip()
	sSample=sTumorFile.split("/")[-1]
	sSample=sSample.split("-T.")[0]
	subprocess.call('''

SNVSniffer somatic -g /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta '''+sOutputDir+'''/Normal_'''+sSample +''' '''+sOutputDir+'''/Tumor_'''+sSample+''' '''+sNormalFile+''' '''+sTumorFile+''' -o '''+sOutputDir+'''/SNVSniffer_'''+sSample+'''


#Remove temporal data
''',shell=True)



#(dGitTNBCFileID,dFileID)=FileDict()
def ParseAlreadylist():
	dAlldict=dict()
	lAlist=glob.glob("/mnt/Hercules/leefall2/Gustave/Somatic/Paired/*.tbi")
	
	for sFile in lAlist:
		sFile=sFile.split("/")[-1]
		sFile=sFile.split(".")[0]
		sSampleID=sFile.split("_")[-1]
		dAlldict[sSampleID]=0
	
	return dAlldict

def producer_task(q, cosmic_dict):

	sFilelist=[]
	
	for sID in dSampleDict.keys():
		if sID[13]=="T":
			sFilelist.append(sID)
	fTarget=open(sWorDir+"/mc3/ABSOULTE_ID.txt")
	lTarget=[]
	for sLine in fTarget.readlines():
		sLine=sLine.strip()
		lTarget.append(sLine)
	#print(dGitTNBCFileID)
	#lBamlist.remove("postrecal_TCGA-A2-A0EQ-01.bam")
	#lBamlist.remove("TCGA-AR-A5QQ-T.bam")
	#for (sDirpath,sDirname,sFilename) in os.walk("."):
	#for sLine in fFile.readlines():
	
#	sFilelist=["TCGA-AR-A5QQ-T.bam"]
	for i in sFilelist:
		#sLine=sLine.strip()
		#t=sLine.split("\t")
		
		#sSampleID="|".join(t)
		if i[0:12] in lTarget:
		
		#value=sSampleID
			value=i
			q.put(value)
			cosmic_dict[value]=None
	#print(dFileDict)
#		else:
			#pass
	#print(dSpecimentoFileDict)

def consumer_task(q, cosmic_dict):
#	print(cosmic_dict)
	while not q.empty():
		value=q.get(True, 0.05)
		#a,b=0,1
#		for item in range(value):
#			a,b=b,a+b


		SamHeader(value)
		SNVSnifferCall(value)
		
		




		cosmic_dict[value]="complete"
		#logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	#os.chdir("/mnt/towel/leefall2/Gustaveraw/EGAF00002385295")
	os.chdir("/mnt/towel/TCGA/Finalpostrecalbam")
	number_of_cpus=cpu_count()-2
	number_of_cpus=10
	manager=Manager()
	fibo_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, fibo_dict))
	producer.start()
	producer.join()
#	
	consumer_list=[]
	for i in range(number_of_cpus):
		consumer=Process(target=consumer_task, args=(data_queue,fibo_dict))
		consumer.start()
		consumer_list.append(consumer)

	[consumer.join() for consumer in consumer_list]

	logger.info(fibo_dict)





#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")








	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
