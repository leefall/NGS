#!/usr/bin/env python
import gzip
import os, subprocess, glob, time
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager

dSampleDict=dict()
#lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_HNSC/postrecal/Preprocess/postrecal_TCGA*.bam")
#lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
#lPostrecal=lpostrecalbam1+lpostrecalbam2
##/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-UF-A7JT-N..bam
#for sPathFile in lPostrecal:
#	sFile=sPathFile.split("/")[-1]
#	sID=sFile.split("_")[1]
#	sID=sID.split(".")[0]
#	if sID in dSampleDict.keys():
#		print(sPathFile)
#		#os.system("rm "+sPathFile)
#	else:
#		dSampleDict[sID]=sPathFile

sInputDir="/mnt/towel/leefall2/GustaveNonTNBC_BRCABam"
sOutputDir="/mnt/towel/leefall2/Gustave_VCF/Non_TNBC_Mutect2"
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




#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $PICARD/ReorderSam.jar \
#INPUT=$BAM/$input".bam" \
#OUTPUT=$Process/$input"_sorted.bam" \
#REFERENCE=$REF/ucsc.hg19.fasta \
#ALLOW_CONTIG_LENGTH_DISCORDANCE=TRUE


def SampleName(sFile):
	subprocess.call('''
	/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" GetSampleName -I '''+str(sFile)+''' -O '''+str(sFile)+'''.samplename''',shell=True)
	



def Mutect2PairedCall(sFile):

#	sID=sFile.split(".")[0]
#	sID=sID.replace("_","")
	#sID=sID.replace("N","")
	sSample=sFile.split("T.")[0]
	sNormalFile=sFile.replace("T.","N.")
#	sSample=sSample.replace("_","")
	#sID=sFile[0:15]
	#sID=sFile.split("_")[1]
	#sID=sFile.split(".")[0]
	#print(sFile+"    "+sSampleID)
	#./EGAF00002385248/s2_B01-007-T1_realigned_recal_sorted.bam|./EGAF00002385247/s2_B01-007-N_realigned_recal_sorted.bam
	#(sSample,sTumorFile,sNormalFile)=sSamplesetID.split("|")
		
	#sNormalFile=dSpecimentoFileDict[sSampleID]["Normal"]
	#sTumorFile=dSpecimentoFileDict[sSampleID]["Tumor"]
	
	
#	IndexBam(sNormalFile)
#	IndexBam(sTumorFile)
	SampleName(sNormalFile)
	fp=open(sNormalFile+'''.samplename''')
	sNormalID=fp.readline()
	sNormalID=sNormalID.strip()
	
	subprocess.call('''
	



/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" Mutect2 \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I '''+str(sFile)+''' \
-I '''+str(sNormalFile)+''' \
--normal '''+str(sNormalID)+''' \
--germline-resource /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
--read-filter NotSecondaryAlignmentReadFilter \
--read-filter NotSupplementaryAlignmentReadFilter \
-L /mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals \
-pon /mnt/QNAP/reference/GATK/hg19/lifted_over_gatk4_mutect2_4136_pon.vcf.gz \
--f1r2-tar-gz /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/f1r2_'''+str(sSample)+'''.tar.gz \
-O /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/rawMutect2_'''+str(sSample)+'''.vcf.gz


#Remove temporal data
''',shell=True)



def PileupSummary(sTumorFile):
#	(sSample,sTumorFile,sNormalFile)=sSamplesetID.split("|")
	
	sSample=sTumorFile.split("T.")[0]
	sNormalFile=sTumorFile.replace("T.","N.")


	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" GetPileupSummaries \
-I '''+str(sTumorFile)+''' \
-V /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals \
-O /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''-T.targeted_sequencing.table 

/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" GetPileupSummaries \
-I '''+str(sNormalFile)+''' \
-V /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals \
-O /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''-N.targeted_sequencing.table


''',shell=True)


def LearnReadOrientationModel(sTumorFile):
	sSample=sTumorFile.split("T.")[0]
	sNormalFile=sTumorFile.replace("T.","N.")

#	(sSample,sTumorFile,sNormalFile)=sSamplesetID.split("|")

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" LearnReadOrientationModel \
-I /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/f1r2_'''+str(sSample)+'''.tar.gz \
-O /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''_read-orientation-model.tar.gz


''',shell=True)


def CalculateContamination(sTumorFile):
	sSample=sTumorFile.split("T.")[0]
	sNormalFile=sTumorFile.replace("T.","N.")

#	(sSample,sTumorFile,sNormalFile)=sSamplesetID.split("|")
	


	subprocess.call('''
	



/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" CalculateContamination \
-I /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''-T.targeted_sequencing.table \
-matched /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''-N.targeted_sequencing.table \
-tumor-segmentation /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''.segments.table \
-O /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''.contamination.table




''',shell=True)
def FilterMutactCall(sTumorFile):
	sSample=sTumorFile.split("T.")[0]
	sNormalFile=sTumorFile.replace("T.","N.")

#	(sSample,sTumorFile,sNormalFile)=sSamplesetID.split("|")

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" FilterMutectCalls \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-V /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/rawMutect2_'''+str(sSample)+'''.vcf.gz \
--contamination-table /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''.contamination.table \
--tumor-segmentation /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''.segments.table \
--orientation-bias-artifact-priors /mnt/Beta/leefall2/Gustave_Non_TNBC_Somatic/'''+str(sSample)+'''_read-orientation-model.tar.gz \
-O /mnt/towel/leefall2/Gustave_VCF/Mutect2_Paired_Non_TNBC/OriFiltered_'''+str(sSample.replace("_",""))+'''.vcf.gz


''',shell=True)




def GATKHaplotypecaller(sFile):
#	sFile=dSampleDict[sTCGANormalID]
	sID=sFile.split(".")[0]
	sID=sID.replace("_","")
	sID=sID.replace("N","")
#	sID=sFile.split("")[1]
#	t=
#	sTemporID=sFile.split(".")[0]
#	sTemporID=sTemporID.split("_")[1]
#	sTemporlist=sTemporID.split("-")
#	sID="-".join(sTemporlist[0:3])
	
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" HaplotypeCaller \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I '''+str(sFile)+''' \
-L /mnt/QNAP/leefall2/Sureselect6_Regions.bed \
-O '''+sOutputDir+'''/'''+str(sID)+'''.vcf.gz \




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
#


def producer_task(q, cosmic_dict):

	#lBamlist=glob.glob("postrecal_TCGA-A2-A0EQ*01.bam")
	lBamlist=glob.glob("*T.bam")
	
	#lBamlist=glob.glob("postrecal_TCGA-*N..bam")
	#lBamlist2=glob.glob("postrecal_TCGA-*11.bam")
	
#	sFilelist=[]
	
#	for sID in dSampleDict.keys():
#		if sID[13]=="N":
#			sFilelist.append(sID)
#	fTarget=open("/mnt/alpha/leefall2/TCGA_HNSC/mc3/allmc3Barcode.txt")
#	lTarget=[]
#	for sLine in fTarget.readlines():
#		sLine=sLine.strip()
#		lTarget.append(sLine)
#	
	
	#lBamlist.remove("postrecal_TCGA-A2-A0EQ-01.bam")

#	lBamlist.remove("postrecal_TCGA-AN-A0XU-01A.bam.bam")
	#sFriday=lBamlist[0:9]
	#sJulia=lBamlist[9:18]
	#sAvatar=lBamlist[18:27]
	#sMeme=lBamlist[27:36]
	#sPiano=lBamlist[36:45]
	#sNeo=lBamlist[45:54]
	#sCipher=lBamlist[54:63]
	#sQueen=lBamlist
#	sFilelist=sFilelist
	#sFilelist=["IonXpress_001_R_2017_07_05_12_17_40_user_S5XL-0059-39-20170704_KPCDx_Test_Auto_user_S5XL-0059-39-20170704_KPCDx_Test_103.bam"]
	for i in lBamlist:
		#value=random.randint(1,len(sFilelist))
		#sCosmicFile=glob.glob("Cosmic_"+i)
		#if not sCosmicFile:
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
		#Mutect2PONcalling(value)
		#CollectOxogArtifact(value)
		#PileupSummary(value)
		#CalculateContamination(value)
		#Mutect2PairedCall(value)
		#FindTumorSampleName(value)
		#Mutect2OnlyTumor(value)
		#FilterVCFcall(value)
		#FilterfromVCFcall(value)
		#FilterBias(value)
		#PileupSummary(value)
		#CalculateContamination(value)
		#LearnReadOrientationModel(value)
		FilterMutactCall(value)




		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir(sInputDir)
#	number_of_cpus=cpu_count()-2
	number_of_cpus=7
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

#	logger.info(fibo_dict)
	print(len(fibo_dict))




#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")








	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
