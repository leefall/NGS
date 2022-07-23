#!/usr/bin/env python
import gzip
import numpy as np
import os, subprocess, glob, time, tabix
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager

dWorkDir="/mnt/Beta/leefall2/TCGA_BLCA_Work"
logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)




class SNV:

	def __init__(self):
		self.sHGVSGenomicChange=""
		self.sHGVSCodingChange=""
		self.sHGVSProteinChange=""
		self.sHGNCGeneSymbol=""
		self.nEntrezGeneID=""
		self.sEnsemblGeneID=""
		self.sChromosome=""
		self.nPosition=""
		self.sReferenceAllele=""
		self.sAlternativeAllele=""
		self.sCytogeneticLocation=""
		self.sStrandOrientation=""
		self.nCodon=""
		self.nExon=""
		self.sMolecularEffects=""
		self.sVariantType=""
		self.sGenotype=""
		self.sdbSNP=""
		self.sClinvar=""
		self.sCosmic=""
		self.sFunctionalDomain=""


	def Parse_ScoreResult(self,sLine):
		sLine=sLine.strip()
		t=sLine.split("\t")

		(self.nSIFT, self.nCADD, self.nPolyPhen)=(t[5],t[22],t[7])




	def Parse_Zygosity(self,sZygosity):
		self.sZygosity=sZygosity


	def Parse_Scorelist(self,lList):
		self.nPolyPhen=lList[7]
		self.nCADD=lList[23]

	def Parse_Cosmic(self,lList):
		sCosmic=lList[1]
		sCosmic=sCosmic.split(";")[0]
		sCosmic=sCosmic.split("=")[-1]
		if self.sCosmic=="":
			self.sCosmic=sCosmic
		else:
			self.sCosmic=self.sCosmic+"|"+sCosmic

	def Parse_dbSNP(self,lList):
		sRSID=lList[1]
		self.sRsID=sRSID




	def Parse_Frequency(self,sLine):
		(sChr, nStartPosition, nEndPosition, sRef, sAlt, nPopFreqMax, n1000G_ALL,n1000G_AFR,n1000G_AMR,n1000G_EAS,n1000G_EUR,n1000G_SAS,nExAC_ALL,nExAC_AFR,nExAC_AMR,nExAC_EAS,nExAC_FIN,nExAC_NFE,nExAC_OTH,nExAC_SAS,nESP6500siv2_ALL,nESP6500siv2_AA,nESP6500siv2_EA,nCG46)=\
		sLine.split("\t")

		self.n1000GALL=n1000G_ALL
		self.n1000GEAS=n1000G_EAS
		#self.nExACALL=nExAC_ALL
		#self.nExACEAS=nExAC_EAS
		#self.nExACSAS=nExAC_SAS
		self.ESP6500ALL=nESP6500siv2_ALL
		self.ESP6500EA=nESP6500siv2_EA

	def Parse_ExAC(self,sExACs):
		(nExACALL, nExACEAS)=(sExACs.split(",")[0],sExACs.split(",")[3])


		self.nExACALL=nExACALL
		self.nExACEAS=nExACEAS
		#self.nExACSAS=nExAC_SAS

	def ExonicParse(self,t):
		#t=[lNMIDs,nExons,sCodons,Allele,sGene]

		lNMIDs=t[0]
		nExons=t[1]
		sCodons=t[2]
		Allele=t[3]
		sGene=t[4]
		lCodonNumber=[]


		if "|" in self.sHGVSProteinChange:
			lRefseqPID=self.sHGVSProteinChange.split("|")
		else:
			lRefseqPID=[self.sHGVSProteinChange]


		for n in xrange(len(sCodons)):
			try:
				sCodons[n]=lRefseqPID[n]+":"+sCodons[n]
			except IndexError:
				sCodons[n]=lRefseqPID[0]+":"+sCodons[n]

		if not lNMIDs==[]:
			for sNMID in lNMIDs:
				sCodon=sNMID.split(":")[-1]
				sCodon=sCodon.replace("c.","")
				if "ins" in sCodon:
					sCodon=sCodon.split("ins")[0]
				elif "del" in sCodon:
					sCodon=sCodon.split("del")[0]
				elif "dup" in sCodon:
					sCodon=sCodon.split("dup")[0]
				else:
					sCodon=sCodon[1:-1]


				lCodonNumber.append(sCodon)


		self.nCodon="|".join(lCodonNumber)
		self.nExon="|".join(nExons)
		self.sHGVSCodingChange="|".join(lNMIDs)
		self.sHGVSProteinChange="|".join(sCodons)




	def parsesCytoband(self,sCytoband):
		self.sCytogeneticLocation=sCytoband




	def Parse_Deleterious(self,lList):
		self.sGene=lList[5]
		self.sMolecularEffect=lList[6]
		self.nSIFT=lList[9]
		self.nCADD=lList[11]
		self.nPolyPhen=lList[10]
		self.n1000GALL=lList[12]
		self.n1000GEAS=lList[13]
		self.nExACALL=lList[14]
		self.nExACEAS=lList[15]
		self.nExACSAS="."
		self.ESP6500ALL="."
		self.ESP6500EA="."
		self.sZygosity=lList[4]
		self.sCosmic=lList[8]
		self.sRsID=lList[7]




def GATKVCFFilter(sFile):
	
	try:
		os.mkdir(dWorkDir+"/Temporaly")
	except:
		pass
	
	try:
		os.mkdir(dWorkDir+"/Filter")
	except:
		pass
	
	
	
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" SelectVariants \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-V '''+sFile+''' \
-select-type SNP \
-O '''+dWorkDir+'''/Temporaly/SNV_'''+sFile+''' 




/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" VariantFiltration \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-V '''+dWorkDir+'''/Temporaly/SNV_'''+sFile+''' \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "DP < 10.0" --filter-name "LowDP10" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--cluster-size 3 \
--cluster-window-size 50 \
-O '''+dWorkDir+'''/Temporaly/Filtered_SNV_'''+sFile+'''






#Remove temporal data


''',shell=True)
	fout=open(dWorkDir+"/Temporaly/Unsorted_"+sFile.replace(".gz",""),"w")
	fSNV=gzip.open(dWorkDir+'''/Temporaly/Filtered_SNV_'''+sFile,"rt")
	
	for sLine in fSNV.readlines():
		if sLine[0]=="#":
			fout.write(sLine)
		else:
			t=sLine.split("\t")
			if t[6]=="PASS":
				fout.write(sLine)
	fSNV.close()
	
#	fINDEL=gzip.open('''./Temporaly/Filtered_INDEL_'''+sFile,"rt")
	
#	for sLine in fINDEL.readlines():
#		if sLine[0]=="#":
#			pass
#		else:
#			t=sLine.split("\t")
#			if t[6]=="PASS":
#				fout.write(sLine)
			
	
	fout.close()
	os.system("vcf-sort "+dWorkDir+"/Temporaly/Unsorted_"+sFile.replace(".gz","")+" >"+dWorkDir+"/Filter/GATK_PASS_"+sFile.replace(".gz",""))
	os.system("bgzip "+dWorkDir+"/Filter/GATK_PASS_"+sFile.replace(".gz",""))
	


def VCF_to_Variant(sFile):
	sTarget="/mnt/towel/BRCA/SureSelect_kinome.bed"
#	sTarget="ureSelect_kinome.bed"

	GATKVCFFilter(sFile)
	#Filter_VCF(sFile,sTarget)

	#ANNOVAR_Annotation(sFile)




def consumer_task(q, cosmic_dict):
	dDeleteriousDict=dict()
	while not q.empty():
		value=q.get(True, 0.05)
		VCF_to_Variant(value)





		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))


def producer_task(q, cosmic_dict):

	sFilelist=glob.glob("*.vcf.gz")
	#sFilelist=glob.glob("TCGA-YZ-A985-10*.vcf.gz")
	sFilelist=sFilelist
	for i in sFilelist:
		value=i
		cosmic_dict[value]=None
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)
#		else:
			#pass





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()

	number_of_cpus=6


	#os.chdir("/second_storage/leefall2/CancerPanel/vcf/Ion/Target_"+sDir)
	os.chdir("/mnt/alpha/leefall2/TCGA_COAD/Germline")
#	try:
#		os.mkdir("Filtering")
#	except:
#		pass


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
