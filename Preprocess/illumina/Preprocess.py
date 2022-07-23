#!/usr/bin/env python
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




def bwa(sFile):
	try:
		os.mkdir("Preprocess")
	except:
		pass

	sID=sFile[0:15]

	subprocess.call('''





samtools fastq '''+sFile+'''|bwa mem -t 8 -R  "@RG\\tID:'''+sID+'''\\tPL:illumina\\tLB:Exome\\tSM:'''+sID+'''" /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta - > /home/leefall2/TCGA_Non_TNBC_BRCA/Realign_'''+sFile+'''

samtools sort /home/leefall2/TCGA_Non_TNBC_BRCA/Realign_'''+sFile+''' -o /home/leefall2/TCGA_Non_TNBC_BRCA/sorted_Realign_'''+sFile+'''


''',shell=True)

def MarkDuplicate(sFile):

	sID=sFile[0:15]



	subprocess.call('''

java -Xmx8G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/picard_tools_1.108/MarkDuplicates.jar \
AS=TRUE \
I=/home/leefall2/TCGA_Non_TNBC_BRCA/sorted_Realign_'''+sFile+''' \
O=/home/leefall2/TCGA_Non_TNBC_BRCA/dedup_sorted_Realign_'''+sFile+''' \
METRICS_FILE=./Preprocess/'''+sID+'''_duplicates \
REMOVE_DUPLICATES=true \
CREATE_INDEX=True


rm /home/leefall2/TCGA_Non_TNBC_BRCA/Realign_'''+sFile+'''


''',shell=True)

def BaseRecalibrate(sFile):

	sID=sFile[0:15]

	subprocess.call('''


### Base Quality Score Recalibration
##  1) Base Recalibrator


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" BaseRecalibrator \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I /home/leefall2/TCGA_Non_TNBC_BRCA/dedup_sorted_Realign_'''+sFile+''' \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
--known-sites /mnt/QNAP/reference/GATK/hg19/dbsnp_138.hg19.vcf \
-O /home/leefall2/TCGA_Non_TNBC_BRCA/recal_'''+str(sID)+'''.table

## 2) Apply Recalibrate

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" ApplyBQSR \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I /home/leefall2/TCGA_Non_TNBC_BRCA/dedup_sorted_Realign_'''+sFile+''' \
--bqsr-recal-file /home/leefall2/TCGA_Non_TNBC_BRCA/recal_'''+str(sID)+'''.table \
-O /mnt/Beta/leefall2/TCGA_LUAD/postrecalBam/postrecal_'''+str(sID)+'''.bam

rm /home/leefall2/TCGA_Non_TNBC_BRCA/sorted_Realign_'''+sFile+'''

''',shell=True)




def RemoveTemporal(sFile):


	sID=sFile[0:15]
	subprocess.call('''

#Remove temporal data

rm /home/leefall2/TCGA_Non_TNBC_BRCA/dedup_sorted_Realign_'''+sFile+''' 

''',shell=True)




def QC(sFile):
	sID=sFile[0:15]

	subprocess.call('''
java -jar /storage/home/leefall2/tools/picard/picard.jar CollectHsMetrics I=/mnt/Beta/leefall2/TCGA_LUAD/postrecalBam/postrecal_'''+str(sID)+'''.bam O=/mnt/Beta/leefall2/TCGA_LUAD/QC/hsmetrics/AlignmentStatus_'''+sID+'''.hs_metrics R=/mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta TARGET_INTERVALS=/mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals BAIT_INTERVALS=/mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals
#samtools flagstat '''+sFile+''' >/mnt/towel/leefall2/QC/Stats/'''+sID+'''_Stats.txt
#samtools depth '''+sFile+'''  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > /mnt/towel/leefall2/QC/Average/'''+sID+'''_Average.txt
	''',shell=True)

def producer_task(q, cosmic_dict):

	lBamlist=glob.glob("TCGA*.bam")
	lPostbam=glob.glob("/mnt/Beta/leefall2/TCGA_LUAD/postrecalBam/postrecal_TCGA-*bam")
	lPostID=[]
	
	
	for sPathFile in lPostbam:
		sFile=sPathFile.split("/")[-1]
		sIDBam=sFile.split("_")[1]
		sID=sIDBam[0:14]
		lPostID.append(sID)
	print(lPostID)
	
	sFilelist=lBamlist
	print("All_Bam")
	print(len(sFilelist))
	

	for i in sFilelist:
		value=i
		
		sID=i[0:14]
		print(sID)
		
		if sID in lPostID:
			pass
		else:
			cosmic_dict[value]=None
			q.put(value)


def consumer_task(q, cosmic_dict):
	while not q.empty():
		value=q.get(True, 0.05)

		bwa(value)
		MarkDuplicate(value)
		BaseRecalibrate(value)
		RemoveTemporal(value)
		QC(value)


		cosmic_dict[value]="complete"





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir("/mnt/alpha/leefall2/TCGA_LUAD/rawbam/Paired/Symbol")
	number_of_cpus=4
	manager=Manager()
	file_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, file_dict))
	producer.start()
	producer.join()
	consumer_list=[]
	print("Number_of_Files")
	print(len(file_dict))
	print(file_dict.keys())
	for i in range(number_of_cpus):
		consumer=Process(target=consumer_task, args=(data_queue,file_dict))
		consumer.start()
		consumer_list.append(consumer)

	[consumer.join() for consumer in consumer_list]

	logger.info(file_dict)





#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")








	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
