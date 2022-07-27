#!/usr/bin/bash
position="$(pwd)"
echo "$position"
for i in KRAS_IonXpress_*.bam;
do
echo $i
#samtools mpileup -uf /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta $position/$i | bcftools view -bvcg - > $i".raw.bcf"
#bcftools view $i".raw.bcf"|vcfutils.pl varFilter -D100000 > $i".vcf"
freebayes -f /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta --region chr12:25362705-25398385 $i  > $i."freebayes.vcf"
done


