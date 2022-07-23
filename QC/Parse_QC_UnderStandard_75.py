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


def CheckUnderStandard(sFile,fout):
	sID=sFile.split("_")[1]
	sID=sID[0:15]
	
	sInput=sFile
	HM=open(sInput)
	
	for i in range(0,7):
		
		HM.readline()				
	HM_line=HM.readline()
	HM_line_list=HM_line.split("\t")

	(nTotalReads, nTotalAlignedRead, ratioAlignedRead ,nMappedbases ,nOntargetBases, nMeanTargetCoverage,nOntargetCoverage,nMedianCoverage)=\
	(HM_line_list[5], HM_line_list[10], HM_line_list[11], HM_line_list[12], HM_line_list[17], HM_line_list[21],HM_line_list[35],HM_line_list[23])
	
	if ((round(float(nOntargetCoverage),2)>=0.75) and (float(nMeanTargetCoverage)>=20)):
		pass
		fout.write("{0}\t{1}\t{2}\n".format(sID,nOntargetCoverage,nMeanTargetCoverage))
	else:
		#print(sFile)
		print("Low Quality Sample\t{0}\t{1}\t{2}".format(sID,nOntargetCoverage,nMeanTargetCoverage))
	
	
	
	
	






if __name__=="__main__":
	StartTime=(time.ctime())
	
	os.chdir("/mnt/alpha/leefall2/TCGA_OV/QC/hsmetrics/")
#	number_of_cpus=cpu_count()-2
	number_of_cpus=4
	p = Pool(number_of_cpus)
	manager=Manager()
	lBamlist=glob.glob("AlignmentStatus_TCGA*.hs_metrics")

	sFilelist=lBamlist

	fout=open("/mnt/alpha/leefall2/TCGA_OV/QC/RemoveFile_from_QC.txt","w")

	for sFile in sFilelist:
		CheckUnderStandard(sFile,fout)
	




	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())
