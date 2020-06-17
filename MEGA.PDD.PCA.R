library("car")
library("survival")
library("MASS")

args<-commandArgs(TRUE)
subjectinput=args[1]
Longitudinalinput=args[2]
sampleinput=args[3]
genotypeinput=args[4]
resultoutput=args[5]

# read Subject data
MEGASubject<-read.delim(subjectinput, header=T, check.names=F)

# read longitude data
MEGALongitude<-read.delim(Longitudinalinput, header=T, check.names=F)

### merge longitude and subject dataframe
MEGALongSub<-merge(MEGALongitude,MEGASubject, by=c("ID","STUDY_NAME"))
MEGALongSub[,"AGEVISIT"]<-with(MEGALongSub,AGE+YEARS)
MEGALongSub[,"DURVISIT"]<-with(MEGALongSub,YEARS+DURBASE)

# read sample name with the same order as genotyping data
MEGAPlate<-read.delim(sampleinput, header=T, check.names=F)

# read genotyping data
MEGAPGenotype<-read.delim(file=genotypeinput, header=F, check.names=F,row.names=1)
colnames(MEGAPGenotype)<-MEGAPlate$ID

###### exclude LS1 dataset, as it use different criteria to determinate "PD dementia"
LongSubSNP<-subset(subset(MEGALongSub,STUDY_NAME!="LS1"),DURVISIT<=12)
DementiaLSub<-subset(LongSubSNP, DEMENTIA==1)
LSDementia<-merge(aggregate(DURVISIT~ID, data=DementiaLSub, min),LongSubSNP,by=c("ID","DURVISIT"))

DementiaLargeLSub<-subset(LongSubSNP, !ID %in% as.character(LSDementia$ID))
DementiaLargeLSub<-subset(DementiaLargeLSub,DEMENTIA==0)
LSDementiaLarge<-merge(aggregate(DURVISIT~ID, data=DementiaLargeLSub, max),LongSubSNP,by=c("ID","DURVISIT"))

LSDementia<-subset(LSDementia,YEARS>0)
SubDementia.raw<-rbind(LSDementia,LSDementiaLarge)

###### For LS1 dataset, as it use different criteria to determinate "PD dementia"
LongSubSNPLS1<-subset(subset(MEGALongSub,STUDY_NAME=="LS1"),DURVISIT<=12)
DementiaLSubLS1<-subset(LongSubSNPLS1, SCOPACOG<=22)
LSDementiaLS1<-merge(aggregate(DURVISIT~ID, data=DementiaLSubLS1, min),LongSubSNPLS1,by=c("ID","DURVISIT"))
LSDementiaLS1$DEMENTIA<-1

DementiaLargeLSubLS1<-subset(LongSubSNPLS1, !ID %in% as.character(LSDementiaLS1$ID))
DementiaLargeLSubLS1<-subset(DementiaLargeLSubLS1, SCOPACOG>22)
LSDementiaLargeLS1<-merge(aggregate(DURVISIT~ID, data=DementiaLargeLSubLS1, max),DementiaLargeLSubLS1,by=c("ID","DURVISIT"))

LSDementiaLargeLS1$DEMENTIA<-0
LSDementiaLS1<-subset(LSDementiaLS1,YEARS>0)
SubDementia.raw.LS1<-rbind(LSDementiaLS1,LSDementiaLargeLS1)

####### all avialable PD subjects for further Cox analysis
SubDementia.raw<-rbind(SubDementia.raw,SubDementia.raw.LS1)


##### Carrier or non-carrier of a allele for each subject
SetCarrier<-function(DF) {
  
  Temp<-DF
  colnames(Temp)<-c("ID","SNP")
  Code2.count=nrow(subset(Temp,SNP == 2))
  Code0.count=nrow(subset(Temp,SNP == 0))
  
  TempCarrier<-subset(Temp,SNP > 0)
  TempCarrier[,"CARRIER"]<-"Y"
  TempNoCarrier<-subset(Temp,SNP == 0)
  TempNoCarrier[,"CARRIER"]<-"N"
  
  Temp<-rbind(TempCarrier,TempNoCarrier)
  Temp$CARRIER<-as.factor(Temp$CARRIER)
  return(Temp)
}

###### Genome-wide association study using Cox regression anlaysis for PD dementia

SNPCoxP<-NULL

for (i in 1:nrow(MEGAPGenotype))  
{
  SNPID<-rownames(MEGAPGenotype)[i]
  MEGATemp<-t(MEGAPGenotype[i,]) 
  MEGATemp<-data.frame(ID=rownames(MEGATemp), SNP=MEGATemp,row.names=NULL,check.names=FALSE)
	  
    tryCatch({
      MEGATemps<-SetCarrier(MEGATemp)
      SubDementia<-merge(SubDementia.raw,MEGATemps, by=c("ID"))

# using survdiff function to check allele effects on PD dementia only	    
      Dementiadif<-survdiff(Surv(DURVISIT,DEMENTIA)~SNP, data=SubDementia)
      difp.val <- 1 - pchisq(Dementiadif$chisq, length(Dementiadif$n) - 1)

# using coxph function to check allele effects on PD dementiaby adjusting covariates
      Dementiafmcox<-coxph(Surv(DURVISIT,DEMENTIA)~AAO+SNP+YEARSEDUC+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+frailty(STUDY_NAME),data=SubDementia, method="breslow")
      Dementiafmcoxp.val<-summary(Dementiafmcox)$coefficients[2,6]
      Dementiafmcox.Conf<-summary(Dementiafmcox)$conf.int[2,]

      Temp<-data.frame(cbind(difp.val,Dementiafmcoxp.val),t(Dementiafmcox.Conf),SNP=SNPID,TN=nrow(SubDementia),CarrierN=nrow(subset(SubDementia,CARRIER=="Y")),NonCarrier=nrow(subset(SubDementia,CARRIER=="N")))
      SNPCoxP<-rbind(SNPCoxP,Temp)},error=function(e){})

}  

colnames(SNPCoxP)<-c("DifPvalue","CoxPvalue","exp(coef)","exp(-coef)","L95","U95","SNP","TotalN","Carrier","NonCarrier")


####### Final Output
if(file.exists(resultoutput))
{
  write.table(SNPCoxP,file=resultoutput,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}else
{
  write.table(SNPCoxP,file=resultoutput,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)  
}


