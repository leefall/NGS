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





def isaac(sFile):
	
	sID=sFile.split("_")[0]
	
	sCurrentDirectory=os.getcwd()
	
	try:
		os.mkdir("Result_"+sID)
	except:
		pass
	
	subprocess.call('''

isaac-align \
-r /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/sorted-reference.xml \
-b '''+sCurrentDirectory+''' \
-m 100 \
--base-calls-format fastq-gz \
-o ./Result_'''+sID+'''


''',shell=True)
	
	
def strelka(sFile):
	#sID=sFile.split("_")[0]
	sID=sFile.split(".")[0]
	try:
		os.mkdir("strelka_"+sID)
	except:
		pass
	
	subprocess.call('''


/storage/home/leefall2/tools/strelka/strelka-2.8.3.centos5_x86_64/bin/configureStrelkaGermlineWorkflow.py \
--bam ./'''+sID+'''.bam \
--referenceFasta /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta \
--callRegions=/storage/home/SNUH/MiSeq_Result/sorted_Targeted_Region.bed.gz \
--runDir ./strelka_'''+sID+'''


./strelka_'''+sID+'''/runWorkflow.py -m local


''',shell=True)



def producer_task(q, cosmic_dict):
	
	#sFilelist=s084+s085+s086+s087
	
	#sFilelist=glob.glob("*_read1.fastq.gz") 
	sFilelist=glob.glob("S150058928.bam") 
	##################### Example of Input File Name####################
	# S150047277_lane1_read1.fastq.gz  S150047277_lane1_read2.fastq.gz #
	####################################################################
	for i in sFilelist:
		value=i
		cosmic_dict[value]=None
	
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)



def consumer_task(q, cosmic_dict):
	while not q.empty():
		value=q.get(True, 0.05)
		isaac(value)
		strelka(value)
		
		#ANNOVAR_Annotation(value)
		
		
		
		
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
#	number_of_cpus=cpu_count()-2
	number_of_cpus=1
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
	
	
	
	
	
	
	
	
	print "Start Time"
	print StartTime
	print "End Time"
	print(time.ctime())





