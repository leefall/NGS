#!/usr/bin/env python
import os
import sys
import glob
import gzip
import commands
import subprocess

sSample=sys.argv[1]
sPackageLocation=sys.argv[2]
try:
	sVcftoolsDir=sys.argv[6]
except IndexError:
	sVcftoolsDir=''
try:
	sPrefix='/'+sys.argv[5]
except IndexError:
	sPrefix=''
sCurrentLocation=os.path.dirname( os.path.abspath( __file__ ) )
sRefDir=sys.argv[4]
sOutputDir=sys.argv[3]

def Make_Interval(sSample):
	pos_list = []
	infile = open(sPackageLocation+"/Process/"+sSample+"_all.vcf")
	lines = infile.readlines()
	
	for line in lines :
		if not line.startswith('#') :
			items = line.strip().split('\t')
			chr = items[0]
			pos = items[1]
			merge = chr + '_' + pos
	
			pos_list.append(merge)
	
	pos_set = set(pos_list)
	print str(len(pos_set))
	
	outfile = open(sPackageLocation+'/Process/'+sSample+'_targets.interval_list','w')
	for pos in pos_set :
		itms = pos.strip().split('_')
		chrom = itms[0]
		posit = itms[1]
		start = int(posit) - 40
		end = int(posit) + 40
	
		outfile.write(chrom + ':' + str(start) + '-' + str(end) + '\n')
	
	outfile.close()

def Make_Code(sSample):
	print "Start Making Code"
	fCount=1
	nCount=1
	
	sCmd=['wc','-l',sPackageLocation+'/Process/Bronj/'+'Target_Amplicon_set.txt']
	fd_popen=subprocess.Popen(sCmd,stdout=subprocess.PIPE).stdout
	data=fd_popen.read().strip()
	nLine=int(data.split(" ")[0])
	nCycle=nLine/11
	sample=sSample
	outfile = open(sCurrentLocation+'/code/' + sample+"_"+str(fCount) + '.sh','w')
	infile1 = open(sPackageLocation+'/Process/Bronj/'+'Target_Amplicon_set.txt')
	lines1 = infile1.readlines()
	for line1 in lines1 :
		outfile.write("java -Xmx10g -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary/ -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/GenomeAnalysisTK_2.8_1_g932cd3a/GenomeAnalysisTK.jar \\" + '\n')
		outfile.write("-R "+sRefDir+"/ucsc.hg19.fasta \\" + '\n')
		outfile.write("-T HaplotypeCaller \\" + '\n')
		outfile.write("-I "+sOutputDir+'/' + sample + ".bam \\" + '\n')
		outfile.write("-o "+sCurrentLocation+"/vcf/" + sample+"_"+str(fCount) + "_" + line1.strip().replace(":","_").replace("-","_") + ".vcf \\" + '\n')
		outfile.write("--dbsnp "+sRefDir+"/dbsnp_137.hg19.vcf \\" + '\n')
		outfile.write("-stand_call_conf 30 \\" + '\n')
		outfile.write("-stand_emit_conf 10 \\" + '\n')
		outfile.write("-bamout "+sCurrentLocation+"/bam/" + sample+"_"+str(fCount)+ "_" + line1.strip().replace(":","_").replace("-","_") + ".bam \\" + '\n')
		outfile.write("--bamWriterType ALL_POSSIBLE_HAPLOTYPES \\" + "\n")
		outfile.write("-L " + line1.strip() + ' \\' + '\n')
		outfile.write("--disable_auto_index_creation_and_locking_when_reading_rods" + '\n\n')
		nCount+=1
		if nCount%nCycle==0:
			outfile.close()
			fCount+=1
			outfile = open(sCurrentLocation+'/code/' + sample+"_"+str(fCount) + '.sh','w')
		else:
			pass
		

def Execute_code(sSample):
	#x="sh "+sCurrentLocation+"/code/"+sSample+"*.sh"
	lFiles_list=glob.glob(sCurrentLocation+"/code/"+sSample+"*.sh")
	for sFile in lFiles_list:
		#print sFile
		x="sh "+sFile+" &"
		#print x
		os.system(x)
	


def Compare_diff_loci(sSample):
	uni_all = []
	for file_uni in sorted(glob.glob(sPackageLocation+'/Process/'+sSample+'*.vcf')) :
		infile_uni = open(file_uni)
		uni_lines = infile_uni.readlines()
		for uni_line in uni_lines :
			if not uni_line.startswith('#') :
				uni_items = uni_line.strip().split('\t')
				uni_chr = uni_items[0]
				uni_pos = uni_items[1]
				uni_merge = uni_chr + '_' + uni_pos
				uni_all.append(uni_merge)
	uni_set = set(uni_all)
	
	recall_all = []
	recall_list_full = []
	
	
	for file in sorted(glob.glob(sCurrentLocation+'/vcf/'+sSample+'*.vcf')) :
		infile = open(file)
		lines = infile.readlines()
		for line in lines :
			if not line.startswith('#') :
				items = line.strip().split('\t')
				chr = items[0]
				pos = items[1]
				ref = items[3]
				alt = items[4]
				recall_list_full.append(chr + '_' + pos + '_' + ref + '_' + alt)
				recall_all.append(chr + '_' + pos)
			elif not line :
				continue
	all_set = set(recall_all)
	cnt = len(uni_set)
	
	ing_vanish = uni_set.difference(all_set)
	ing_intersection = uni_set.intersection(all_set)
	ing_trueNegative = all_set.difference(uni_set)
	
	outfile = open(sCurrentLocation+'/TXT/'+sSample+'_uni_vs_Ldenovo_stat.txt','w')
	outfile.write('Type	FalsePositive	TruePositive	FalseNegative 	FP(%) \n')
	outfile.write("TOTAL" + '\t' + str(len(ing_vanish)) + '\t' + str(len(ing_intersection)) + '\t' + str(len(ing_trueNegative)) + '\t' + str(round(float(len(ing_vanish))/float(cnt)*100,2)) + '\n')
	
	
	outfile1 = open(sCurrentLocation+'/TXT/'+sSample+'_FalsePositive.txt','w')
	for i in ing_vanish :
		outfile1.write(i + '\n')
	outfile2 = open(sCurrentLocation+'/TXT/'+sSample+'_TruePositive.txt','w')
	for i in ing_intersection :
		outfile2.write(i + '\n')
	outfile3 = open(sCurrentLocation+'/TXT/'+sSample+'_FalseNegative.txt','w')
	for i in ing_trueNegative :
		outfile3.write(i + '\n')
	
	
	all_set_full = set(recall_list_full)
	
	indel_cnt = 0
	snv_cnt = 0
	out_indel = open(sCurrentLocation+'/TXT/'+sSample+'_FalseNegative_indel.txt','w')
	out_snv = open(sCurrentLocation+'/TXT/'+sSample+'_FalseNegative_snv.txt','w')
	for full_item in all_set_full :
		items = full_item.strip().split('_')
		chr_pos = items[0] + '_' + items[1]
		if chr_pos in ing_trueNegative :
			ref = items[2]
			alt = items[3]
	
			if len(list(ref)) > 1 or len(list(alt)) > 1 :
				indel_cnt = indel_cnt + 1
				out_indel.write(full_item.strip() + '\n')
			else :
				snv_cnt = snv_cnt + 1
				out_snv.write(full_item.strip() + '\n')
	
	outfile.write("FN INDEL cnt : " + str(indel_cnt) + '\n')
	outfile.write("FN SNV cnt : " + str(snv_cnt) + '\n')
	
	outfile.close()


def Make_UniqVCF(sSample):
	os.system("cat "+sCurrentLocation+"/vcf/" + sSample + '*.vcf | grep -v "#" | sort | uniq > '+sCurrentLocation+'/vcf_after_true/' + sSample + '.vcf')
	
	

def Sort_VCFOne(sSample):
	for file in sorted(glob.glob(sCurrentLocation+'/vcf_after_true/*'+sSample+'.vcf')):
		outname = file.strip().split('/')[-1]
		os.system("cat " + file + " | "+sVcftoolsDir+"vcf-sort -c > "+sCurrentLocation+"/vcf_after_true/sorted_" + outname)
def Make_TrueVCF(sSample):
	
	
	true_pos = open(sCurrentLocation+'/TXT/'+sSample+'_TruePositive.txt')
	true_lines = true_pos.readlines()
	
	infile1 = open(sCurrentLocation+'/TXT/'+sSample+'_FalseNegative_snv.txt')
	lines1 = infile1.readlines()
	FN_snv_list = []
	for line1 in lines1 :
		FN_snv_list.append(line1.strip().split('_')[0] + '_' + line1.strip().split('_')[1])
	
	
	true_list = []
	for true_line in true_lines :
		true_list.append(true_line.strip())
	
	outfile = open(sCurrentLocation+'/F_vcf/' + sSample + '_1.vcf','w')

	ld_pos_list = []
	vcf_list = []
	infile_ld = open(sCurrentLocation+'/vcf_after_true/sorted_' + sSample + '.vcf')
	ld_lines = infile_ld.readlines()
	for ld_line in ld_lines :
		chr = ld_line.strip().split('\t')[0]
		pos = ld_line.strip().split('\t')[1]
		ref = ld_line.strip().split('\t')[3]
		alt = ld_line.strip().split('\t')[4]
		merge = chr + '_' + pos
		if merge in true_list :
			if merge not in ld_pos_list :
				ld_pos_list.append(chr + '_' + pos)
				vcf_list.append(ld_line.strip())
		elif merge in FN_snv_list :
			if merge not in ld_pos_list :
				if len(list(ref)) < 2 and len(list(alt)) < 2 :
					ld_pos_list.append(chr + '_' + pos)
					vcf_list.append(ld_line.strip())

	for i in vcf_list :
		outfile.write(i + '\n')
	outfile.close()
	
	
def Sort_VCFTwo(sSample):
	file=sSample+'_1.vcf'
	os.system("cat "+sCurrentLocation+"/F_vcf/" + file + " | vcf-sort -c > "+sCurrentLocation+"/F_vcf/" + sSample + ".vcf")
	print "Sort VCFTwo finished"

def Bgzip(sSample):
	file=sSample+".vcf"
	os.chdir(sCurrentLocation+"/F_vcf/")
	os.system("bgzip " + file)
	os.system("tabix -p vcf " + sSample + '.vcf.gz')
	os.system("mv "+file+".gz* "+sOutputDir+'/'+sPrefix)
	print "mv "+file+".gz* "+sOutputDir+'/'+sPrefix
	print "Bgzip finished"


def Main():
#	sSample=sys.argv[1]
	Make_Interval(sSample)
	#print "Done Before Make_Code"
	Make_Code(sSample)
	#print "Done Before Execute"
	Execute_code(sSample)
	Compare_diff_loci(sSample)
	Make_UniqVCF(sSample)
	Sort_VCFOne(sSample)
	Make_TrueVCF(sSample)
	Sort_VCFTwo(sSample)
	Bgzip(sSample)

	
Main()
