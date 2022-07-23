#!/usr/bin/bash
for sFile in *.bam;
do
samtools flagstat $sFile > $sFile.stats
done
