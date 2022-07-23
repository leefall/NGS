#!/usr/bin/env python
import os
import gzip
import subprocess
import sys, time, random, re ,requests, logging, glob, tabix
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
import multiprocessing

def Cosmic_Filter(value, queue):
	print(value)
	os.system("Rscript "+value)
	return
	


def QCmetrics(sFile):
	sID=sFile.split("/")[-1]
	sID=sID.split(".")[0]
	
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk CollectHsMetrics -I '''+sFile+''' -O /mnt/towel/leefall2/QC/HSmetrics/AlignmentStatus_'''+sID+'''.hs_metrics -BI /mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals -TI /mnt/QNAP/leefall2/Sureselect6_target_interval_file.intervals

#samtools flagstat '''+sFile+''' >/mnt/towel/leefall2/QC/Stats/'''+sID+'''_Stats.txt
#samtools depth '''+sFile+'''  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > /mnt/towel/leefall2/QC/Average/'''+sID+'''_Average.txt
	''',shell=True)
	

if __name__=="__main__":
	StartTime=(time.ctime())

	number_of_cpus=4
	manager=Manager()

	p = Pool(number_of_cpus)
	#sFilelist=glob.glob("/mnt/towel/UKBiobank/FE_vcfs/*.gvcf.gz")#4011848
	
	sFinalFilelist=glob.glob("/mnt/towel/leefall2/GustaveTNBCBam/*.bam")
	
	sFinalFilelist=sFinalFilelist

	Pool_map = p.map(QCmetrics, sFinalFilelist)
	
	
	
	
	
	
	

	print(StartTime)
	print((time.ctime()))
	

