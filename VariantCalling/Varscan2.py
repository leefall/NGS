#!/usr/bin/env python
import gzip
import os, subprocess, glob, time
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager

sCurrent=os.getcwd()
lCurrent=sCurrent.split("/")
sWorDir="/".join(lCurrent[0:-2])
sCancerType=lCurrent[-3]
#sInputDir="/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/"

dSampleDict=dict()
lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_HNSC/postrecal/Preprocess/postrecal_TCGA*.bam")
lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
lPostrecal=lpostrecalbam1+lpostrecalbam2

sOutputDir="/data/leefall2//TCGA_HNSC_Varscan2/"



#dSampleDict=dict()
#sCurrent=os.getcwd()
#lCurrent=sCurrent.split("/")
#sWorDir="/".join(lCurrent[0:-2])
#sCancerType=lCurrent[-3]
#sInputDir="/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/"
#sOutputDir="/scratch/leefall2/TCGA_Non_TNBC_BRCA_Varscan2/"
#lpostrecalbam1=glob.glob("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/postrecalBam/Preprocess/postrecal_TCGA*.bam")
#lpostrecalbam2=glob.glob("/mnt/Beta/leefall2/TCGA_HNSC/postrecalbam/postrecal_TCGA-*.bam")
#lPostrecal=lpostrecalbam1+lpostrecalbam2
#lPostrecal=lpostrecalbam1

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
	sID=sFile[0:15]
	#sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	#sTail=sTail.replace("1","2")

	#sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''





samtools fastq '''+sFile+'''|bwa mem -t 8 -R  "@RG\\tID:'''+sID+'''\\tPL:illumina\\tLB:Exome\\tSM:'''+sID+'''" /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta - > ./Preprocess/Realign_'''+sFile+'''

samtools sort ./Preprocess/Realign_'''+sFile+''' -o ./Preprocess/sorted_Realign_'''+sFile+'''

rm ./Preprocess/Realign_'''+sFile+'''

''',shell=True)

def MarkDuplicate(sFile):

	sID=sFile[0:15]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN



	subprocess.call('''

java -Xmx8G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/picard_tools_1.108/MarkDuplicates.jar \
AS=TRUE \
I=./Preprocess/sorted_Realign_'''+sFile+''' \
O=./Preprocess/dedup_sorted_Realign_'''+sFile+''' \
METRICS_FILE=./Preprocess/'''+sID+'''_duplicates \
REMOVE_DUPLICATES=true \
CREATE_INDEX=True


''',shell=True)

def BaseRecalibrate(sFile):

	sID=sFile[0:15]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN

	subprocess.call('''


### Base Quality Score Recalibration
##  1) Base Recalibrator


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" BaseRecalibrator \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I ./Preprocess/dedup_sorted_Realign_'''+sFile+''' \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
--known-sites /mnt/QNAP/reference/GATK/hg19/dbsnp_138.hg19.vcf \
-O ./Preprocess/recal_'''+str(sID)+'''.table

## 2) Apply Recalibrate

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" ApplyBQSR \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I ./Preprocess/dedup_sorted_Realign_'''+sFile+''' \
--bqsr-recal-file ./Preprocess/recal_'''+str(sID)+'''.table \
-O ./Preprocess/postrecal_'''+str(sID)+'''.bam

''',shell=True)



def IntersectBed(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''



bedtools intersect -a ../Preprocess/postrecal_'''+str(sID)+'''.bam -b /mnt/QNAP/leefall2/S31285117_Regions.bed > ../Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam

samtools index ../Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam


#Remove temporal data



rm ./Preprocess/'''+sID+'''.bam
#rm ./Preprocess/'''+sID+'''.sam
rm ./Preprocess/sorted'''+sID+'''.bam
rm ./Preprocess/dedup_sorted'''+sID+'''.bam


''',shell=True)

def IndexBam(sFile):
	subprocess.call('''
	samtools index '''+str(sFile),shell=True)


def SampleName(sFile):
	subprocess.call('''
	/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" GetSampleName -I '''+str(sFile)+''' -O '''+str(sFile)+'''.samplename''',shell=True)
	




def GATKHaplotypecaller(sSampleID,dSpecimentoFileDict):

	#sID=sFile[0:15]
	#sID=sFile.split("_")[1]
	#sID=sFile.split(".")[0]
	#print(sFile+"    "+sSampleID)
	
	
	sNormalFile=dSpecimentoFileDict[sSampleID]["Normal"]
	sTumorFile=dSpecimentoFileDict[sSampleID]["Tumor"]
	
	
	#IndexBam(sNormalFile)
	#IndexBam(sTumorFile)
	#SampleName(sNormalFile)
	#subprocess.call('''
	print('''



/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" Mutect2 \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I '''+str(sTumorFile)+''' \
-I '''+str(sNormalFile)+''' \
--normal '''+str(sNormalFile)+'''.samplename \
--germline-resource /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
--read-filter NotSecondaryAlignmentReadFilter \
--read-filter NotSupplementaryAlignmentReadFilter \
-L /mnt/towel/BRCA/SureSelect_kinome.bed \
-pon /mnt/QNAP/reference/GATK/hg19/lifted_over_gatk4_mutect2_4136_pon.vcf.gz \
--f1r2-tar-gz /mnt/Hercules/leefall2/Gustave/Preprocess/f1r2_'''+str(sSampleID)+'''.tar.gz \
-O /mnt/Hercules/leefall2/Gustave/Somatic/Paired/rawMutect2_'''+str(sSampleID)+'''.vcf.gz


#Remove temporal data

''')
#''',shell=True)
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



#''',shell=True)


def Mutect2PairedCall(sSamplesetID):

	#sID=sFile[0:15]
	#sID=sFile.split("_")[1]
	#sID=sFile.split(".")[0]
	#print(sFile+"    "+sSampleID)
	#./EGAF00002385248/s2_B01-007-T1_realigned_recal_sorted.bam|./EGAF00002385247/s2_B01-007-N_realigned_recal_sorted.bam
	(sSample,sTumorFile,sNormalFile)=sSamplesetID.split("|")
	
	#sNormalFile=dSpecimentoFileDict[sSampleID]["Normal"]
	#sTumorFile=dSpecimentoFileDict[sSampleID]["Tumor"]
	
	
	IndexBam(sNormalFile)
	IndexBam(sTumorFile)
	SampleName(sNormalFile)
	fp=open(sNormalFile+'''.samplename''')
	sNormalID=fp.readline()
	sNormalID=sNormalID.strip()
	
	subprocess.call('''
	



/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" Mutect2 \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I '''+str(sTumorFile)+''' \
-I '''+str(sNormalFile)+''' \
--normal '''+str(sNormalID)+''' \
--germline-resource /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
--read-filter NotSecondaryAlignmentReadFilter \
--read-filter NotSupplementaryAlignmentReadFilter \
-L /mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals \
-pon /mnt/QNAP/reference/GATK/hg19/lifted_over_gatk4_mutect2_4136_pon.vcf.gz \
--f1r2-tar-gz /mnt/Hercules/leefall2/Gustave/Preprocess/f1r2_'''+str(sSample)+'''.tar.gz \
-O /mnt/Hercules/leefall2/Gustave/Somatic/Paired/rawMutect2_'''+str(sSample)+'''.vcf.gz


#Remove temporal data
''',shell=True)





def Mutect2PONcalling(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" Mutect2 \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/Ophthalmology/Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
--panel-of-normals /mnt/towel/Ophthalmology/Preprocess/OPLL_pon.vcf.gz \
-L /mnt/QNAP/leefall2/Sureselect7_target_interval_file.intervals \
-O /mnt/towel/Ophthalmology/VCF/Mutect2/Mutect2_'''+str(sID)+'''.vcf.gz



#Remove temporal data


''',shell=True)



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
-I ../Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O ../Preprocess/oxog_'''+str(sID)+''' \
--FILE_EXTENSION .txt



''',shell=True)

def PileupSummary(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail


	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" GetPileupSummaries \
-I ../Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O ../Preprocess/'''+str(sID)+'''.targeted_sequencing.table \
-V /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/QNAP/leefall2/S31285117_Regions.bed \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta


''',shell=True)

def CalculateContamination(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail


	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" CalculateContamination \
-I ../Preprocess/'''+str(sID)+'''.targeted_sequencing.table \
-O ../Preprocess/'''+str(sID)+'''.targeted_sequencing.contamination.table


''',shell=True)


def FindTumorSampleName(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail


	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" GetSampleName \
-I ../Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O ../Preprocess/'''+str(sID)+'''.targeted_sequencing.sample_name


''',shell=True)

def Mutect2OnlyTumor(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" Mutect2 \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/Ophthalmology/Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O /mnt/towel/Ophthalmology/VCF/Mutect2/Mutect2TumorOnly_'''+str(sID)+'''.vcf.gz \
-tumor ../Preprocess/'''+str(sID)+'''.targeted_sequencing.sample_name \
--af-of-alleles-not-in-resource 2.5e-06 \
--germline-resource /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/QNAP/leefall2/Sureselect7_target_interval_file.intervals \
-pon /mnt/QNAP/reference/GATK/hg19/lifted_over_gatk4_mutect2_4136_pon.vcf.gz


#Remove temporal data


''',shell=True)



def FilterVCFcall(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" FilterMutectCalls \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-O /mnt/towel/Ophthalmology/VCF/Mutect2/Mutect2TumorOnly_'''+str(sID)+'''.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \
-V /mnt/towel/Ophthalmology/VCF/Mutect2/Mutect2TumorOnly_'''+str(sID)+'''.vcf.gz \
--contamination-table ../Preprocess/'''+str(sID)+'''.targeted_sequencing.contamination.table \
-L /mnt/QNAP/leefall2/Sureselect7_target_interval_file.intervals



#Remove temporal data


''',shell=True)

def FilterfromVCFcall(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail
	
	fp=gzip.open("/mnt/towel/Ophthalmology/VCF/Mutect2/Mutect2TumorOnly_"+str(sID)+".targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz","rb")
	fout=open("/mnt/towel/Ophthalmology/VCF/Mutect2/Filtered_Mutect2TumorOnly_"+str(sID)+".vcf","w")
	
	
	for sLine in fp.readlines():
		if sLine[0]=="#":
			fout.write(sLine)
		else:
			t=sLine.split("\t")
			
			
			sFilter=t[6]
			
			if sFilter=="PASS":
				fout.write(sLine)
			else:
				pass
	
	
	fout.close()
	subprocess.call('''

bgzip '''+"/mnt/towel/Ophthalmology/VCF/Mutect2/Filtered_Mutect2TumorOnly_"+str(sID)+".vcf"+'''
tabix -p vcf '''+"/mnt/towel/Ophthalmology/VCF/Mutect2/Filtered_Mutect2TumorOnly_"+str(sID)+".vcf.gz"+'''


''',shell=True)


def FilterBias(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" FilterByOrientationBias \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-O /mnt/towel/Ophthalmology/VCF/Mutect2/FinalMutect2TumorOnly_'''+str(sID)+'''.vcf.gz \
-P ../Preprocess/oxog_'''+str(sID)+'''.pre_adapter_detail_metrics.txt \
-V /mnt/towel/Ophthalmology/VCF/Mutect2/Filtered_Mutect2TumorOnly_'''+str(sID)+'''.vcf.gz \
-L /mnt/QNAP/leefall2/Sureselect7_target_interval_file.intervals \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-AM G/T \
-AM C/T



#Remove temporal data


''',shell=True)




def mPileup(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]
	subprocess.call('''


samtools mpileup -q 1 -f /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta -B '''+sNormalBam+''' '''+sFile+''' > '''+sOutputDir+'''/Two_'''+sID+'''pileup

#Remove temporal data


''',shell=True)


def Varscan2(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]

#	sID=sFile.split("-T.")[0]
#	sNormalBam=sFile.replace("T.","N.")

	subprocess.call('''

java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/tools/VarScan.v2.3.9.jar somatic '''+sOutputDir+'''/Two_'''+sID+'''pileup '''+sOutputDir+'''/Varscan2_'''+sID+''' --mpileup 1
#Remove temporal data


''',shell=True)

def Removetempor(sTumorID):
	sFile=dSampleDict[sTumorID]
	sID=sTumorID[0:12]
	sNormalBam=dSampleDict[sID+"-N"]

#	sID=sFile.split("-T.")[0]
#	sNormalBam=sFile.replace("T.","N.")

	subprocess.call('''


rm '''+sOutputDir+'''/Two_'''+sID+'''pileup

#Remove temporal data


''',shell=True)
	



def FileDict():
#dCliTNBC=dict()
	sCliTNBC=set()
	sGitTNBC=set()
	dGitTNBCFileID=dict()
	dFileID=dict()
	fGit=open("/storage/home/leefall2/mypro/BRCA_EGA/Gustave/mBC_WES_Fabrice_Andre_2019_private-master/Data/Corresponding_ID_EGA_mBC.txt","r")
	fClinic=open("/storage/home/leefall2/mypro/BRCA_EGA/Gustave/mBC_WES_Fabrice_Andre_2019_private-master/Data/mbc_sample_info_withoutWhitespaceatLastline.txt")
	
	
	
	
	
	
	
	for sLine in fClinic.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sID,sType)=(t[2],t[3])
		if sType=="TNBC":
			sCliTNBC.add(sID[:-1])
			
		
	#print(sCliTNBC)
	#print(len(sCliTNBC))
	
	fGit.readline()
	
	for sLine in fGit.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#sGitTNBC.add(t[-1].replace("_",""))
		
		sSpecimenID=t[-1].replace("_","")
		
		sSampleID=sSpecimenID[:-1]
		
		
		if sSampleID in sCliTNBC:
			if sSpecimenID in dGitTNBCFileID.keys():
				pass
			else:
				dGitTNBCFileID[sSpecimenID]=[]
			
			
			if "N" in sSpecimenID:
				#dGitTNBCFileID[t[1]]=t[-1]
				dFileID[t[1]]=sSpecimenID
				dGitTNBCFileID[sSpecimenID].append("Normal")
				dGitTNBCFileID[sSpecimenID].append(t[1])
				dGitTNBCFileID[sSpecimenID].append(sSampleID)
				#dGitTNBCFileID[sSpecimenID]["Normal"].append(t[1])
				#dGitTNBCFileID[sSampleID]["Normal"].append(sSpecimenID)
			elif "T" in sSpecimenID:
				dFileID[t[1]]=sSpecimenID
				#dGitTNBCFileID[t[1]]=t[-1]
				#dGitTNBCFileID[sSpecimenID]["Tumor"].append(t[1])
				dGitTNBCFileID[sSpecimenID].append("Tumor")
				dGitTNBCFileID[sSpecimenID].append(t[1])
				dGitTNBCFileID[sSpecimenID].append(sSampleID)
				#dGitTNBCFileID[sSampleID]["Tumor"].append(sSpecimenID)
				
		
		
		
	#print(dGitTNBCFileID)
	#print(len(dGitTNBCFileID))
	#print (len(sCliTNBC&sGitTNBC))
	
	
	
	return dGitTNBCFileID,dFileID


#(dGitTNBCFileID,dFileID)=FileDict()


def producer_task(q, cosmic_dict):

	#lBamlist=glob.glob("postrecal_TCGA-A2-A0EQ*01.bam")
	
	sFilelist=[]
	
	for sID in dSampleDict.keys():
		if sID[13]=="T":
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
	#print(cosmic_dict)
	while not q.empty():
		value=q.get(True, 0.05)
		#a,b=0,1
#		for item in range(value):
#			a,b=b,a+b


#		
		mPileup(value)
		Varscan2(value)
		Removetempor(value)
#		
		#Mutect2PONcalling(value)
		#CollectOxogArtifact(value)
		#PileupSummary(value)
		#CalculateContamination(value)
		#FindTumorSampleName(value)
		#Mutect2OnlyTumor(value)
		#FilterVCFcall(value)
		#FilterfromVCFcall(value)
		#FilterBias(value)



		cosmic_dict[value]="complete"
		#logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	#os.chdir("/mnt/towel/leefall2/Gustaveraw/EGAF00002385295")
	#os.chdir("/mnt/towel/TCGA/Finalpostrecalbam")
	number_of_cpus=cpu_count()-2
	number_of_cpus=6
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

	#logger.info(fibo_dict)
	print(len(fibo_dict))




#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")








	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
