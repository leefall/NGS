#!/usr/bin/bash


### Make target intervals from target bed using PICARD ###
java -jar /storage/home/leefall2/tools/picard/picard.jar BedToIntervalList I=/mnt/QNAP/leefall2/noChr_S07604514_Regions.bed O=/mnt/QNAP/leefall2/noChrSureselect6_target_interval_file.intervals SD=/mnt/QNAP/reference/GATK/hg19/b37/TCGA/Homo_sapiens_assembly19.fasta



