#!/usr/bin/env python
import sys, os, subprocess, time, commands, glob
from optparse import OptionParser
from subprocess import PIPE,Popen
print(time.ctime())

#Get input file as argument sInputFile
parser=OptionParser(usage = "usage: %prog -i Input_File(without '.bam') -b Input Directory -g GATK_Directory name -r Reference_Directory -p Picard_Directory ([Option]-o Output Directory -s Samtools_Directroy -v Vcftools_Directory)")
parser.add_option("-i","--Input_File", dest="Input_file",
		help="Your Input File (Only File name without '.bam', and Neglect this parameter if you want to run batch)", default='')
parser.add_option("-b","--Input_Dir", dest="Input_Dir",
		help="Your Input_[REQUIRED]", default='Input')
parser.add_option("-g","--GATK_Dir", dest="GATK_Dir",
		help="Your GATK Directory [REQUIRED]", default='')		
parser.add_option("-p","--Picard_Dir", dest="Picard_Dir",
		help="Your Picard_Directory [REQUIRED]", default='')		
parser.add_option("-r","--ref_Dir", dest="Ref_Dir",
		help="Your Reference_Directory for align [REQUIRED]", default='')		
parser.add_option("-s","--Samtools_Dir", dest="Samtools_Dir",
		help="Your Samtools_Directory[Option]", default='')		
parser.add_option("-v","--Vcf_Dir", dest="Vcf_Dir",
		help="Your vcftools_Directory[Option]", default='')		
parser.add_option("-o","--Output_Dir", dest="Output_Dir",
		help="Your Output_Directory[Option] 'current Directory/Result' is defalut", default='')				
parser.add_option("-f","--Prefix", dest="Prefix",
		help="Your Output file name[Option]", default='')

(options, args)=parser.parse_args()
sInputFile=options.Input_file
sInputDir=options.Input_Dir
sGATKDir=options.GATK_Dir
sPicardDir=options.Picard_Dir
sRefDir=options.Ref_Dir
sSamtoolsDir=options.Samtools_Dir
sVcftoolsDir=options.Vcf_Dir
sOutputDir=options.Output_Dir
sPrefix=options.Prefix
CurrentLocation=os.getcwd()
PackageLocation=os.path.dirname( os.path.abspath( __file__ ) )
if not options.Output_Dir:
	os.system('mkdir ./Result')
	sOutputDir=CurrentLocation+'/Result'
else:
	pass
	#print options.Output_Dir
	#sys.exit()

	

if not options.GATK_Dir:
	print "Missing param: GATK Directory"
	parser.print_usage()
	exit()
if not options.Picard_Dir:
	print "Missing param: Picard Directory"
	parser.print_usage()
	exit()
if not options.Ref_Dir:
	print "Missing param: Reference Directory"
	parser.print_usage()
	exit()
if options.Samtools_Dir:
	sSamtoolsDir=sSamtoolsDir+'/'
if options.Vcf_Dir:
	sVcftoolsDir=sVcftoolsDir+'/'
	
if not options.Input_file:
	if sInputDir=="Input":
		print "Missing param: Input file name or Input Directory for batch"
		parser.print_usage()
		exit()
	else:
		os.chdir(sInputDir)
		lFile_list=glob.glob("*.bam")
		os.chdir(CurrentLocation)
		for sFile in lFile_list:
			sFile=sFile.replace(".bam","")
			comm=PackageLocation+'/F_code.py -i '+sFile+' -b '+sInputDir+' -g '+sGATKDir+' -p '+sPicardDir+' -r '+sRefDir+" -o "+sOutputDir+' &'
			os.system(comm)
else:
	print(time.ctime())
	sInputFile=sInputFile.replace(".bam","")
#	os.chdir(PackageLocation+'/code')
#	os.chdir(PackageLocation+'/code')
	command='python Curation_Bronj_one.py '+sInputFile+' '+PackageLocation+' '+sOutputDir+' '+sRefDir+' '+sPrefix+' '+sVcftoolsDir
	os.system(command)
	os.chdir(PackageLocation+'/Process')
	print "Finish"
	print(time.ctime())

