#!/usr/bin/env python
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



def Samtoolsort(sFile):

	
	#print(sFile)
	
	sID=sFile[0:12]
	subprocess.call('''

samtools sort -n -o ../sorted/'''+sID+'''.bam '''+sFile+'''


''',shell=True)

def Bamtofastq(sFile):

	
	sID=sFile[0:12]
	subprocess.call('''

bedtools bamtofastq -i ../sorted/'''+sID+'''.bam -fq ../FASTQ/'''+sID+'''.end1.fq -fq2 ../FASTQ/'''+sID+'''.end2.fq

''',shell=True)



def producer_task(q, cosmic_dict):
	
	
	sFilelist=glob.glob("*.bam")
	#sFilelist=["IonXpress_001_R_2015_11_05_14_47_48_user_PRO-61-20151104_Cancer_Core_Target_001_Auto_user_PRO-61-20151104_Cancer_Core_Target_001_95.txt"]
	for i in sFilelist:
		#value=random.randint(1,len(sFilelist))
		#sCosmicFile=glob.glob("Cosmic_"+i)
		#if not sCosmicFile:
		value=i
		cosmic_dict[value]=None
	
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)
#		else:
			#pass

def consumer_task(q, cosmic_dict):
	#dBarcode_Dict=Parse_Barcode_Dictionary()
	while not q.empty():
		value=q.get(True, 0.05)

		#CalculateDepth("IonXpress_001_R_2015_11_05_14_47_48_user_PRO-61-20151104_Cancer_Core_Target_001_Auto_user_PRO-61-20151104_Cancer_Core_Target_001_95.txt")
		#STARFusion(value)
		Samtoolsort(value)
		Bamtofastq(value)
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))



if __name__=="__main__":

	#CalculateDepth("IonXpress_001_R_2015_11_05_14_47_48_user_PRO-61-20151104_Cancer_Core_Target_001_Auto_user_PRO-61-20151104_Cancer_Core_Target_001_95.txt")
	os.chdir("/mnt/backup/TNBC_RNA/Symbol")
	StartTime=(time.ctime())
	data_queue=Queue()
##	number_of_cpus=cpu_count()-2
	number_of_cpus=8
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
	