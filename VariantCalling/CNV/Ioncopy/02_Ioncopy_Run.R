#Ioncopy
install.packages("ioncopy")
library(ioncopy)
BiocManager::install("multtest")
#Load data

Amplicon.data<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/Header_GC_Contects_per_merged_Bed.txt",header=FALSE,stringsAsFactors = FALSE))
head(Amplicon.data)
Target_001_Depth<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_001.txt",header=FALSE,stringsAsFactors = FALSE))
Target_002_Depth<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_002.txt",header=FALSE,stringsAsFactors = FALSE))
Target_003_Depth<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_003.txt",header=FALSE,stringsAsFactors = FALSE))
Target_004_Depth<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_004.txt",header=FALSE,stringsAsFactors = FALSE))
Target_005_Depth<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_005.txt",header=FALSE,stringsAsFactors = FALSE))
Target_010_Depth<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_010.txt",header=FALSE,stringsAsFactors = FALSE))


Target_001_ID<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_001_ID.txt",header=FALSE,stringsAsFactors = FALSE))
Target_002_ID<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_002_ID.txt",header=FALSE,stringsAsFactors = FALSE))
Target_003_ID<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_003_ID.txt",header=FALSE,stringsAsFactors = FALSE))
Target_004_ID<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_004_ID.txt",header=FALSE,stringsAsFactors = FALSE))
Target_005_ID<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_005_ID.txt",header=FALSE,stringsAsFactors = FALSE))
Target_010_ID<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Average_Original/Target_010_ID.txt",header=FALSE,stringsAsFactors = FALSE))

#Add Column
colnames(Amplicon.data)<-c("Start_Position","End_Position","Amplicon","Gene","GC_Contents")



Ioncopy.header<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/Code/Ioncopy_Header.txt",header = FALSE,stringsAsFactors = FALSE))
head(Ioncopy.header)
head(Target_001_Depth)

#Batch 1
Ioncopy_Batch.1<-Target_001_Depth

row.names(Ioncopy_Batch.1)<-Ioncopy.header[,1]

colnames(Ioncopy_Batch.1)<-Target_001_ID[,1]

head(Ioncopy_Batch.1)

#Significant test
CNV.1<-assess.CNA(Ioncopy_Batch.1)
CNV.1$CN
#bonferroni
vanilla.1<-call.CNA(CNV.1,method.mt="none")
call.vanilla.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0)
call.vanilla.1<-call.CNA(CNV.1)
call.vanilla.1$tab

write.table(call.vanilla.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch1_Bonferroni.txt",quote=FALSE,sep="\t")

#FDR
call.fdr.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="fdr")

call.fdr.1$tab
write.table(call.fdr.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch1_FDR.txt",quote=FALSE,sep="\t")


#Target002
Ioncopy_Batch.1<-Target_002_Depth

row.names(Ioncopy_Batch.1)<-Ioncopy.header[,1]
head(Target_002_ID)
colnames(Ioncopy_Batch.1)<-Target_002_ID[,1]

head(Ioncopy_Batch.1)

#Significant test
CNV.1<-assess.CNA(Ioncopy_Batch.1)
CNV.1$CN
#bonferroni
vanilla.1<-call.CNA(CNV.1,method.mt="none")
call.vanilla.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0)
call.vanilla.1<-call.CNA(CNV.1)
#call.vanilla.1$tab

write.table(call.vanilla.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch2_Bonferroni.txt",quote=FALSE,sep="\t")

#FDR
call.fdr.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="fdr")

#call.fdr.1$tab
write.table(call.fdr.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch2_FDR.txt",quote=FALSE,sep="\t")


#Target_003
Ioncopy_Batch.1<-Target_003_Depth

row.names(Ioncopy_Batch.1)<-Ioncopy.header[,1]
head(Target_003_ID)
colnames(Ioncopy_Batch.1)<-Target_003_ID[,1]

head(Ioncopy_Batch.1)

#Significant test
CNV.1<-assess.CNA(Ioncopy_Batch.1)
CNV.1$CN
#bonferroni
vanilla.1<-call.CNA(CNV.1,method.mt="none")
call.vanilla.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0)
call.vanilla.1<-call.CNA(CNV.1)
#call.vanilla.1$tab

write.table(call.vanilla.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch3_Bonferroni.txt",quote=FALSE,sep="\t")

#FDR
call.fdr.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="fdr")

#call.fdr.1$tab
write.table(call.fdr.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch3_FDR.txt",quote=FALSE,sep="\t")


###Target 004

Ioncopy_Batch.1<-Target_004_Depth

row.names(Ioncopy_Batch.1)<-Ioncopy.header[,1]
head(Target_004_ID)
colnames(Ioncopy_Batch.1)<-Target_004_ID[,1]

head(Ioncopy_Batch.1)

#Significant test
CNV.1<-assess.CNA(Ioncopy_Batch.1)
CNV.1$CN
#bonferroni
vanilla.1<-call.CNA(CNV.1,method.mt="none")
call.vanilla.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0)
call.vanilla.1<-call.CNA(CNV.1)
#call.vanilla.1$tab

write.table(call.vanilla.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch4_Bonferroni.txt",quote=FALSE,sep="\t")

#FDR
call.fdr.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="fdr")

#call.fdr.1$tab
write.table(call.fdr.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch4_FDR.txt",quote=FALSE,sep="\t")

###Target 005

Ioncopy_Batch.1<-Target_005_Depth

row.names(Ioncopy_Batch.1)<-Ioncopy.header[,1]
head(Target_005_ID)
colnames(Ioncopy_Batch.1)<-Target_005_ID[,1]

head(Ioncopy_Batch.1)

#Significant test
CNV.1<-assess.CNA(Ioncopy_Batch.1)
CNV.1$CN
#bonferroni
vanilla.1<-call.CNA(CNV.1,method.mt="none")
call.vanilla.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0)
call.vanilla.1<-call.CNA(CNV.1)
#call.vanilla.1$tab

write.table(call.vanilla.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch5_Bonferroni.txt",quote=FALSE,sep="\t")

#FDR
call.fdr.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="fdr")

#call.fdr.1$tab
write.table(call.fdr.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch5_FDR.txt",quote=FALSE,sep="\t")




###Target 004

Ioncopy_Batch.1<-Target_010_Depth

row.names(Ioncopy_Batch.1)<-Ioncopy.header[,1]
head(Target_010_ID)
colnames(Ioncopy_Batch.1)<-Target_010_ID[,1]

head(Ioncopy_Batch.1)

#Significant test
CNV.1<-assess.CNA(Ioncopy_Batch.1)
CNV.1$CN
#bonferroni
vanilla.1<-call.CNA(CNV.1,method.mt="none")
call.vanilla.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0)
call.vanilla.1<-call.CNA(CNV.1)
#call.vanilla.1$tab

write.table(call.vanilla.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch10_Bonferroni.txt",quote=FALSE,sep="\t")

#FDR
call.fdr.1<-call.CNA(CNV.1, analysis.mode="gene-wise", method.p="samples_genes/amplicons",method.mt="fdr")

#call.fdr.1$tab
write.table(call.fdr.1$tab,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/Ioncopy/Batch10_FDR.txt",quote=FALSE,sep="\t")


#Shniy

runIoncopy()
?assess.CNA



#With Control.

.libPaths()


