#!/usr/bin/env python
sInputBam=""
sHSmetcisPostBam=""

import gzip
import os, subprocess, glob, time
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager


logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)




def producer_task(q, cosmic_dict):

	lBamlist=glob.glob(sInputBam+"/*.bam")
	lPostrecalBam1=glob.glob(dHSmetcisPostBam+"/AlignmentStatus_*.hs_metrics")
	lPostbam=lPostrecalBam1
	lPostID=[]
	
	
	for sPathFile in lPostbam:
		sFile=sPathFile.split("/")[-1]
		sIDBam=sFile.split("_")[1]
		sID=sIDBam[0:14]
		lPostID.append(sID)
	dTargetTCGA=dict()
	fFinal=open("/mnt/alpha/leefall2/TCGA_LUAD/mc3/allmc3Barcode.txt")
	fFinal.readline()
	
	
	
	for sLine in fFinal.readlines():
		sLine=sLine.strip()
		sTCGAID=sLine.split("\t")[0]
		#sFinal.add(t[0])
		dTargetTCGA[sTCGAID+"-T"]=0
		dTargetTCGA[sTCGAID+"-N"]=0
	print("All_Bam")
	print(len(lBamlist))
	
	sFilelist=[]

	for i in lBamlist:
		sFile=i.split("_")[1]
		sID=sFile[0:14]
		if sID in dTargetTCGA:
			if sID in lPostID:
				pass
			else:
				sFilelist.append(i)

	print("Non process files")

	print(len(sFilelist))
	sFinalFilelist=sFilelist

	for i in sFinalFilelist:

		value=i
		cosmic_dict[value]=None

		q.put(value)


def consumer_task(q, cosmic_dict):
	while not q.empty():
		value=q.get(True, 0.05)
		cosmic_dict[value]="complete"




if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir("/mnt/Beta/leefall2/TCGA_LUAD/postrecalBam/")
#	number_of_cpus=cpu_count()-2
	number_of_cpus=4
	manager=Manager()
	fibo_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, fibo_dict))
	producer.start()
	producer.join()
	consumer_list=[]
	print("Number_of_Files_for_process")
	print(len(fibo_dict))
	for i in range(number_of_cpus):
		consumer=Process(target=consumer_task, args=(data_queue,fibo_dict))
		consumer.start()
		consumer_list.append(consumer)

	[consumer.join() for consumer in consumer_list]








	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
