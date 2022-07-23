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



def CnvCoverage(sFile):

	sID=sFile.split("_")[0]
	
	subprocess.call('''


cnvkit.py coverage '''+sPath+'''/'''+str(sFile)+''' /mnt/towel/Ophthalmology/CNV/sorted_unique_Intersected_sorted_refFlat.bed -o /mnt/towel/Ophthalmology/CNV/reference/'''+str(sID)+'''.targetcoverage.cnn


#Remove temporal data


''',shell=True)


def cnrtocns(sFile):
	#
	#sID=sFile.split("_")[-1]
	sID=sFile.split(".")[0]
	subprocess.call('''


cnvkit.py segment '''+str(sFile)+''' -o cns/'''+str(sID)+'''.cns


#Remove temporal data


''',shell=True)




def producer_task(q, cosmic_dict):
	
	bamfiles=glob.glob("*.cnr")

	
		
	
	
	sFilelist=bamfiles
	
	#sFilelist=bamfiles[235:276]
	#sFilelist=["IonXpress_001_R_2017_07_05_12_17_40_user_S5XL-0059-39-20170704_KPCDx_Test_Auto_user_S5XL-0059-39-20170704_KPCDx_Test_103.bam"]
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
	while not q.empty():
		value=q.get(True, 0.05)
		#a,b=0,1
#		for item in range(value):
#			a,b=b,a+b
		
		
		#bwa(value)
		#MarkDuplicate(value)
		#BaseRecalibrate(value)
#		IntersectBed(value)
		#GATKHaplotypecaller(value)
		#CnvCoverage(value)
		#
		cnrtocns(value)
		
		
		
		
		
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	#sPath="/mnt/towel/Ophthalmology/CNV/results/flasso/"
	#os.chdir(sPath)
#	number_of_cpus=cpu_count()-2
	number_of_cpus=6
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
	
	
	
	
	
#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")
	
	
	
	
	
	
	
	
	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())





