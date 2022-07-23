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
sOutputDir="/selene/leefall2/TCGA_Non_TNBC_Sequenza/"
lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/postrecalBam/Preprocess/postrecal_TCGA*.bam")
#lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
#lPostrecal=lpostrecalbam1+lpostrecalbam2
lPostrecal=lpostrecalbam1

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




def Preprocess(sTumorID):
	#sNormalFile=sTumorFile.replace("T.","N.")
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]
	#SampleName(sNormalFile)
	#fp=open('''/data/leefall2/TCGA/Preprocess/'''+sNormalFile+'''.samplename''')
	#sNormalID=fp.readline()
	#sNormalID=sNormalID.strip()
	#sSample=sTumorFile.split("-T.")[0]
	subprocess.call('''



sequenza-utils bam2seqz -n '''+sNormalBam+''' -t '''+sFile+''' --fasta /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta -gc /mnt/QNAP/leefall2/hg19.gc50Base.wig.gz -o '''+sOutputDir+'''/out.'''+sID+'''.gz

sequenza-utils seqz_binning --seqz /mnt/towel/TCGA/CNV/Sequenza/Preprocessed/out.'''+sID+'''.gz -w 50 -o '''+sOutputDir+'''/out_small.'''+sID+'''.gz


''',shell=True)



def RemoveTemporal(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]
	
	#SampleName(sNormalFile)
	#fp=open('''/data/leefall2/TCGA/Preprocess/'''+sNormalFile+'''.samplename''')
	#sNormalID=fp.readline()
	#sNormalID=sNormalID.strip()
	
	subprocess.call('''

rm '''+sOutputDir+'''/out.'''+sID+'''.gz


''',shell=True)




def NonChr(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]
	
	
	#fout=open("/mnt/towel/TCGA/CNV/Sequenza/Rscript/Script_"+sSample+".R","w")
	fp=gzip.open(sOutputDir+"/out_small."+sID+".gz","rt")
	fout=open(sOutputDir+"/nonChr_out_small."+sID,"w")
	fout.write(fp.readline())
	
	lChr=[]
	for i in range(1,23):
		lChr.append(str(i))
	lChr.append("X")
	
	
	for sLine in fp.readlines():
		t=sLine.split("\t")
		t[0]=t[0].replace("chr","")
		if t[0] in lChr:
			fout.write("{0}".format("\t".join(t)))
		
	
	fout.close()
	
	subprocess.call('''
	

bgzip '''+sOutputDir+'''/nonChr_out_small.'''+sID+'''

''',shell=True)




def MakeRscript(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]
	
	fout=open(sOutputDir+"/Script_"+sID+".R","w")
	
	
	
	#subprocess.call('''
	fout.write('''

library(sequenza)
bb<-sequenza.extract("'''+sOutputDir+'''/nonChr_out_small.'''+sID+'''.gz",verbose=TRUE)
CP<-sequenza.fit(bb)
sequenza.results(sequenza.extract = bb, cp.table = CP, sample.id = "'''+sID+'''",out.dir="'''+sOutputDir+'''/'''+sID+'''")
''')
	
	fout.close()
	subprocess.call('''
Rscript '''+sOutputDir+'''/Script_'''+sID+'''.R
''',shell=True)






def producer_task(q, cosmic_dict):

	
	sFilelist=[]
	
	for sID in dSampleDict.keys():
		if sID[13]=="T":
			#print(sID)
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

		
		Preprocess(value)
		RemoveTemporal(value)
		NonChr(value)
		MakeRscript(value)




		cosmic_dict[value]="complete"







if __name__=="__main__":
	
	
	StartTime=(time.ctime())
	data_queue=Queue()
	#os.chdir("/mnt/towel/leefall2/Gustaveraw/EGAF00002385295")
	#os.chdir("/mnt/towel/TCGA/Finalpostrecalbam")
	#number_of_cpus=cpu_count()-2
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
	
	
	
	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
	





