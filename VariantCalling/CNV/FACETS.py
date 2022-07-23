#!/usr/bin/env python
import os, subprocess, glob, time
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
dSampleDict=dict()
lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_HNSC/postrecal/Preprocess/postrecal_TCGA*.bam")
lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
lPostrecal=lpostrecalbam1+lpostrecalbam2
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

#print(dSampleDict["TCGA-UF-A7JT-N"])
#print(dSampleDict["TCGA-UF-A7JT-T"])





sInputDir="/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/"
sOutputDir="/home/leefall2/TCGA_HNSC_CNV/"
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





	
def bwa(sFile):
#	fout=open("Code_"+sFile+".txt","w")
#	try:
#		os.mkdir("Preprocess")
#	except:
#		pass
#	try:
#		os.mkdir("Final_Bam")
#	except:
#		pass
	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''
	
	

bwa mem -t 8 -R "@RG\\tID:'''+sID+'''\\tPL:illumina\\tLB:Exome\\tSM:'''+sID+'''" /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta '''+sFile+''' '''+sPartner+''' > /mnt/towel/BRCA/Preprocess/'''+sID+'''.sam

samtools view -bhS /mnt/towel/BRCA/Preprocess/'''+sID+'''.sam > /mnt/towel/BRCA/Preprocess/'''+sID+'''.bam
samtools sort /mnt/towel/BRCA/Preprocess/'''+sID+'''.bam -o /mnt/towel/BRCA/Preprocess/sorted'''+sID+'''.bam


''',shell=True)

def MarkDuplicate(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	
	subprocess.call('''

java -Xmx8G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/picard_tools_1.108/MarkDuplicates.jar \
AS=TRUE \
I=/mnt/towel/BRCA/Preprocess/sorted'''+sID+'''.bam \
O=/mnt/towel/BRCA/Preprocess/dedup_sorted'''+sID+'''.bam \
METRICS_FILE=/mnt/towel/BRCA/Preprocess/'''+sID+'''_duplicates \
REMOVE_DUPLICATES=true \
CREATE_INDEX=True


''',shell=True)

def BaseRecalibrate(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''


### Base Quality Score Recalibration
##  1) Base Recalibrator


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" BaseRecalibrator \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I /mnt/towel/BRCA/Preprocess/dedup_sorted'''+sID+'''.bam \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
--known-sites /mnt/QNAP/reference/GATK/hg19/dbsnp_138.hg19.vcf \
-O /mnt/towel/BRCA/Preprocess/recal_'''+str(sID)+'''.table

## 2) Apply Recalibrate

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" ApplyBQSR \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I /mnt/towel/BRCA/Preprocess/dedup_sorted'''+sID+'''.bam \
--bqsr-recal-file /mnt/towel/BRCA/Preprocess/recal_'''+str(sID)+'''.table \
-O /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sID)+'''.bam


rm /mnt/towel/BRCA/Preprocess/'''+sID+'''.sam


''',shell=True)



def IntersectBed(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''
	


bedtools intersect -a /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sID)+'''.bam -b /mnt/towel/BRCA/SureSelect_kinome.bed > /mnt/towel/BRCA/Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam

samtools index /mnt/towel/BRCA/Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam


#Remove temporal data



#rm ./Preprocess/'''+sID+'''.bam
#rm ./Preprocess/'''+sID+'''.sam
#rm ./Preprocess/sorted'''+sID+'''.bam
#rm ./Preprocess/dedup_sorted'''+sID+'''.bam 


''',shell=True)



def FACETS(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID.split("-T")[0]
	#sID=sID.split("
	sNormalBam=dSampleDict[sID+"-N"]
	#fout=open("/mnt/towel/TCGA/Code/CNV/Exome_rest_"+sID+".sh","w")
	#print('''/mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam''')
	#print('''/mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam''')
	subprocess.call('''
	
cnv_facets.R \
-t '''+str(sFile)+''' \
-n '''+str(sNormalBam)+''' \
-vcf /mnt/QNAP/reference/GATK/hg19/dbsnp_138.hg19.vcf.gz \
--targets /mnt/QNAP/leefall2/noheader_Sureselect6_Regions.bed \
-g hg19 \
-o '''+sOutputDir+'''/'''+str(sID)+'''



#Remove temporal data


''',shell=True)
#--targets /mnt/QNAP/leefall2/Sureselect6_Regions.bed \

def producer_task(q, cosmic_dict):
	
	
	
	#sFilelist=glob.glob("postrecal_TCGA-*-T..bam")
	
	#dSampleDict[sID]=sPathFile
	
	sFilelist=[]
	
	for sID in dSampleDict.keys():
		if sID[13]=="T":
			sFilelist.append(sID)
	fTarget=open("/mnt/alpha/leefall2/TCGA_HNSC/mc3/allmc3Barcode.txt")
	lTarget=[]
	for sLine in fTarget.readlines():
		sLine=sLine.strip()
		lTarget.append(sLine)
	
	
	lAlready=[]
	lOutputFiles=glob.glob(sOutputDir+"/*vcf.gz")
	for sPathFile in lOutputFiles:
		sFile=sPathFile.split("/")[-1]
		sTCGAID=sFile.split(".")[0]
		lAlready.append(sTCGAID)
		
		
	#print("Total Bam")
	#print(len(sFilelist))
	#sFilelist=sFilelist[0:1]
	
	#sFilelist=["IonXpress_001_R_2017_07_05_12_17_40_user_S5XL-0059-39-20170704_KPCDx_Test_Auto_user_S5XL-0059-39-20170704_KPCDx_Test_103.bam"]
	for i in sFilelist:
		
		#sID=i.split("_")[1]
		#sID=sID.split("-T.")[0]
		sID=i
		#value=random.randint(1,len(sFilelist))
		#sCosmicFile=glob.glob("Cosmic_"+i)
		
		if not sID[0:12] in lTarget:
			#pass
			print(i)
		
		
		if sID in lAlready:
			pass
		else:
			value=i
			cosmic_dict[value]=None
		
			#logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
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
		#IntersectBed(value)
		#GATKHaplotypecaller(value)
		#Mutect2(value)
		FACETS(value)
		
		
		
		cosmic_dict[value]="complete"
		#logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir(sInputDir)
#	number_of_cpus=cpu_count()-2
	number_of_cpus=4
	manager=Manager()
	fibo_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, fibo_dict))
	producer.start()
	producer.join()
	consumer_list=[]
	print("Number of Files to Process")
	print(len(fibo_dict))
	for i in range(number_of_cpus):
		consumer=Process(target=consumer_task, args=(data_queue,fibo_dict))
		consumer.start()
		consumer_list.append(consumer)

	[consumer.join() for consumer in consumer_list]

	#logger.info(fibo_dict)
	
	
	
	
	
#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")
	
	
	
	
	
	
	
	
	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())



#
#
