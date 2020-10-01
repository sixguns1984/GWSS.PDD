library("plyr")

args<-commandArgs(TRUE)
genotypeinput=args[1]
HRinput=args[2]
PHSoutput=args[3]

#read genotyping data of subject, subject ID as columan name, SNPID as row name
MEGAPGenotype<-read.delim(file=genotypeinput, header=F, check.names=F,row.names=1)

#read HR dataset with SNPID and HR two columns
MEGADementia.HR<-read.delim(file=HRinput, header=T, check.names=F)

###########
#calculate PHS function
########### 
PHSscore<- function(MEGAPGenotype, MEGADementia.HR)  
{
  IndivdiualPHS<-NULL
  for (i in 1:ncol(MEGAPGenotype))
  {
  
  Tempdataframe<-data.frame(SNPs=MEGAPGenotype[,i],MEGADementia.HR[,c("HR","SNPID")])
  
  NASNP<-sum(is.na(Tempdataframe[,"SNPs"]))
  
  Comb.PHS<-sum(with(Tempdataframe,SNPs*log(HR)),na.rm = T)
  
  IndivdiualPHS<-rbind(IndivdiualPHS,data.frame(ID=colnames(MEGAPGenotype)[i],Comb.PHS,SNPNumber=nrow(MEGADementia.HR)-NASNP))
  }
  return(IndivdiualPHS)
}


Subject.PHS<-PHSscore(MEGAPGenotype, MEGADementia.HR)

write.table(Subject.PHS,file=PHSoutput,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")



