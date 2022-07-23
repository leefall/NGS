#!/usr/bin/env python
import gzip
import os, subprocess, glob, time
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager

dSampleDict=dict()
lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_HNSC/postrecal/Preprocess/postrecal_TCGA*.bam")
lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
lPostrecal=lpostrecalbam1+lpostrecalbam2
for sPathFile in lPostrecal:
	sFile=sPathFile.split("/")[-1]
	sID=sFile.split("_")[1]
	sID=sID.split(".")[0]
	if sID in dSampleDict.keys():
		print(sPathFile)
		#os.system("rm "+sPathFile)
	else:
		dSampleDict[sID]=sPathFile

sInputDir="/mnt/alpha/leefall2/TCGA_HNSC/postrecal/Preprocess"
sOutputDir="/scratch/leefall2/TCGA_HNSC_VCF"
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





def GATKHaplotypecaller(sTCGANormalID):
	sFile=dSampleDict[sTCGANormalID]
	sID=sTCGANormalID[0:12]

	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" HaplotypeCaller \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I '''+str(sFile)+''' \
-L /mnt/QNAP/leefall2/Sureselect6_Regions.bed \
-O '''+sOutputDir+'''/'''+str(sID)+'''.vcf.gz \




#Remove temporal data


''',shell=True)


def producer_task(q, cosmic_dict):


	
	sFilelist=[]
	
	for sID in dSampleDict.keys():
		if sID[13]=="N":
			sFilelist.append(sID)
	fTarget=open("/mnt/alpha/leefall2/TCGA_HNSC/mc3/allmc3Barcode.txt")
	lTarget=[]
	for sLine in fTarget.readlines():
		sLine=sLine.strip()
		lTarget.append(sLine)
	

	sFilelist=sFilelist
	for i in sFilelist:

		value=i
		cosmic_dict[value]=None

		#logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)
#		else:
			#pass


def consumer_task(q, cosmic_dict):
	while not q.empty():
		value=q.get(True, 0.05)

		GATKHaplotypecaller(value)



		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir(sInputDir)
#	number_of_cpus=cpu_count()-2
	number_of_cpus=3
	manager=Manager()
	fibo_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, fibo_dict))
	producer.start()
	producer.join()
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
