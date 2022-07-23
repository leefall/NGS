#!/bin/sh
#export GATK=$(cd ..; pwd)/Tools/GenomeAnalysisTK_2.8_1_g932cd3a
#export PICARD=$(cd ..; pwd)/Tools/picard_tools_1.108
#export REF=$(cd ..; pwd)/data/hg19
#export SAMTOOLS=$(cd ..; pwd)/Tools/samtools_0.1.19
#export TMP=$(cd ..; pwd)/temporary
#export BAM=$(cd ..; pwd)/Input
#export Process=$(cd ..; pwd)/Process
export GATK=$3
export PICARD=$4
export REF=$5
export SAMTOOLS=$9
export TMP=$8
export BAM=$2
export Process=$6






export GATK, PICARD, REF, TMP, BAM, SAMTOOLS, Process

input=$1
sam=samtools
InputDir=$2
OutputDir=$7
#date
#export BAM=$(cd ..; pwd)/$InputDir



#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $PICARD/ReorderSam.jar \
#INPUT=$BAM/$input".bam" \
#OUTPUT=$Process/$input"_sorted.bam" \
#REFERENCE=$REF/ucsc.hg19.fasta \
#ALLOW_CONTIG_LENGTH_DISCORDANCE=TRUE

#$SAMTOOLS$sam index $Process/$input"_sorted.bam"

#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-R $REF/ucsc.hg19.fasta \
#-T RealignerTargetCreator \
#-known $REF/Mills_and_1000G_gold_standard.indels.hg19.vcf \
#-I $Process/$input"_sorted.bam" \
#-o $Process/$input"_realigner.intervals"
#
#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-T IndelRealigner \
#-R $REF/ucsc.hg19.fasta \
#-I $Process/$input"_sorted.bam" \
#-targetIntervals $Process/$input"_realigner.intervals" \
#-o $Process/$input"_sorted_realigned.bam"

#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-T BaseRecalibrator \
#-R $REF/ucsc.hg19.fasta \
#-I $Process/$input"_sorted_realigned.bam"  \
#-knownSites $REF/dbsnp_137.hg19.vcf \
#-o $Process/$input"_sorted_realigned_recal.grp"

#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-T BaseRecalibrator \
#-R $REF/ucsc.hg19.fasta \
#-I $Process/$input"_sorted_realigned.bam" \
#-knownSites $REF/dbsnp_137.hg19.vcf \
#-BQSR $Process/$input"_sorted_realigned_recal.grp" \
#-o $Process/$input"_sorted_realigned_post.grp"

#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-T PrintReads \
#-R $REF/ucsc.hg19.fasta \
#-I $Process/$input"_sorted_realigned.bam" \
#-BQSR $Process/$input"_sorted_realigned_post.grp" \
#--out $Process/$input"_sorted_realigned_recaled.bam"

#intersectBed -abam $Process/$input"_sorted_realigned_recaled.bam" -b sequencing_target.bed > $OutputDir/$input".bam"

$SAMTOOLS$sam index $OutputDir/$input".bam"

java -Xmx10G -Djava.io.tmpdir=$TMP \
-jar $GATK/GenomeAnalysisTK.jar \
-R $REF/ucsc.hg19.fasta \
-T UnifiedGenotyper \
-I $OutputDir/$input".bam" \
-o $OutputDir/$input"_all.vcf" \
--dbsnp $REF/dbsnp_137.hg19.vcf \
-stand_call_conf 30 \
-stand_emit_conf 10 \
-glm BOTH

#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-R $REF/ucsc.hg19.fasta \
#-T UnifiedGenotyper \
#-I $OutputDir/$input".bam" \
#-o $OutputDir/$input"_5.vcf" \
#--dbsnp $REF/dbsnp_137.hg19.vcf \
#-stand_call_conf 30 \
#-stand_emit_conf 5 \
#-glm BOTH


#java -Xmx10G -Djava.io.tmpdir=$TMP \
#-jar $GATK/GenomeAnalysisTK.jar \
#-R $REF/ucsc.hg19.fasta \
#-T UnifiedGenotyper \
#-I $OutputDir/$input".bam" \
#-o $OutputDir/$input"_1.vcf" \
#--dbsnp $REF/dbsnp_137.hg19.vcf \
#-stand_call_conf 30 \
#-stand_emit_conf 1 \
#-glm BOTH
