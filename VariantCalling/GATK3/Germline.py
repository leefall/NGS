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



def GATKHaplotypecaller(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/BRCA/Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O /mnt/towel/BRCA/VCF/Sureselect6_vcf_'''+str(sID)+'''.vcf.gz 



#Remove temporal data


''',shell=True)
#java -Xmx2G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/GenomeAnalysisTK_2.8_1_g932cd3a/GenomeAnalysisTK.jar \
#-R /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta \
#-T UnifiedGenotyper \
#-I ./Preprocess/recal_realign_dedup_sorted_Reordered_'''+str(sID)+'''.bam \
#-o ./GATK_'''+str(sID)+'''.vcf \
#--dbsnp /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/dbsnp_137.hg19.vcf \
#-stand_call_conf 30 \
#-stand_emit_conf 10 \
#-glm BOTH

def CollectOxogArtifact(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail


	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" CollectSequencingArtifactMetrics \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sID)+'''.bam \
-O /mnt/towel/BRCA/Preprocess/Preprocess/oxog_'''+str(sID)+''' \
--FILE_EXTENSION .txt



''',shell=True)


def Mutect2(sFile):

	sID=sFile.split("-")[0]
	sID=sID.split("_")[1]
	
	
	#sID=sID+"-"+sTN
	#sTail=sTail.replace("1","2")
	
	#sPartner=sFile.split("_")[0]+"_"+sTail
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	
	#print('''/mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam''')
	#print('''/mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam''')
	

#/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" Mutect2 \
#--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
#-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
#-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam \
#-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam \
#--normal-sample '''+str(sNormalFile)+''' \
#--germline-resource /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
#-L /mnt/towel/BRCA/SureSelect_kinome.intervals \
#-pon /mnt/QNAP/reference/GATK/hg19/lifted_over_gatk4_mutect2_4136_pon.vcf.gz \
#-O /mnt/towel/BRCA/VCF/rawVCF/rawMutect2_'''+str(sID)+'''.vcf.gz

	subprocess.call('''
#Remove temporal data



java -Xmx8G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T MuTect2 \
-R /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta \
-I:tumor /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam \
-I:normal /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam \
--dbsnp /mnt/QNAP/reference/GATK/hg19//dbsnp_138.hg19.vcf \
--cosmic /mnt/QNAP/reference/Cosmic/Cosmic.hg19.vcf \
-L /mnt/towel/BRCA/SureSelect_kinome.intervals \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o /mnt/towel/BRCA/VCF/GATK3/Somatic_'''+str(sID)+'''.vcf.gz



''',shell=True)


def PileupSummary(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	#sTail=sTail.replace("1","2")
	
	#sPartner=sFile.split("_")[0]+"_"+sTail
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" GetPileupSummaries \
-I ../Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O ../Preprocess/'''+str(sID)+'''.targeted_sequencing.table \
-V /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/QNAP/leefall2/S07604514_Regions.bed \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta


''',shell=True)



def producer_task(q, cosmic_dict):
	
	
	#s089=glob.glob("PMSNH089*-T_1.fastq.gz") #10
	#s090=glob.glob("PMSNH090*-T_1.fastq.gz") #10
	#s089=glob.glob("PMSNH089*-N_1.fastq.gz") #4
	#s090=glob.glob("PMSNH090*-N_1.fastq.gz") #6
	
	#s091=glob.glob("PMSNH091*-T_1.fastq.gz")
	#s092=glob.glob("PMSNH092*_1.fastq.gz")#13
#	s077=glob.glob("postrecal_PMSNH097*-T.bam")#11
#	s078=glob.glob("postrecal_PMSNH100*-T.bam")#11
	s097=glob.glob("postrecal_PMSNH086*-T.bam")#11
	s098=glob.glob("postrecal_PMSNH088*-T.bam")#11
	s099=glob.glob("postrecal_PMSNH010*-T.bam")#11
#	s085=glob.glob("PMSNH085*_1.fastq.gz")#11
	#s099=glob.glob("PMSNH099*-T_1.fastq.gz")
	#s096N=glob.glob("PMSNH096*-N_1.fastq.gz")
	#s097N=glob.glob("PMSNH097*-N_1.fastq.gz")
	
	
#	sPiano=s077+s078
	sPiano=s097+s098+s099
	sFilelist=sPiano
	
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
		#IntersectBed(value)
		#GATKHaplotypecaller(value)
		Mutect2(value)
		#CollectOxogArtifact(value)
		#PileupSummary(value)
		#CalculateContamination(value)
		#FindTumorSampleName(value)
		#Mutect2OnlyTumor(value)
		#FilterVCFcall(value)
		#FilterfromVCFcall(value)
		#FilterBias(value)
		
		
		
		
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir("/mnt/towel/BRCA/Preprocess")
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
	
	
	
	
	
#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")
	
	
	
	
	
	
	
	
	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())





