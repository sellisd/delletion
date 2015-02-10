#read required libraries
library(stringr)
#read required files
conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")
dubious<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)

minBarNo <- 1 # >minBarNo should be available for each gene 
eoL <- numeric(0) # number of extreme overdominant in each condition
hetBenL <- numeric(0)
totalL <- numeric(0)
counter <- 1

# remove from dubious not trusted
notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)
neutral <- dubious[which(!dubious[,1] %in% therio),1] # neutral set without not trusted

for(i in conditions$V1){
  ua <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ua.dat",sep=""),header=F)
  da <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".da.dat",sep=""),header=F)
  us <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".us.dat",sep=""),header=F)
  ds <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ds.dat",sep=""),header=F)
  ua$V1 <- str_extract(ua$V1,"[^:]*")
  us$V1 <- str_extract(us$V1,"[^:]*")
  da$V1 <- str_extract(da$V1,"[^:]*")
  ds$V1 <- str_extract(ds$V1,"[^:]*") 
  #neutral set
  ua.ind<-which(ua$V1%in%dubious$dubiousORFS.Whitlock)
  da.ind<-which(da$V1%in%dubious$dubiousORFS.Whitlock)
  us.ind<-which(us$V1%in%dubious$dubiousORFS.Whitlock)
  ds.ind<-which(ds$V1%in%dubious$dubiousORFS.Whitlock)
  
  ua.hetD<-quantile(ua$V3[ua.ind],probs=c(.05,.95),na.rm=T,names=F)
  ua.homD<-quantile(ua$V2[ua.ind],probs=c(.05,.95),na.rm=T,names=F)
  da.hetD<-quantile(da$V3[da.ind],probs=c(.05,.95),na.rm=T,names=F)
  da.homD<-quantile(da$V2[da.ind],probs=c(.05,.95),na.rm=T,names=F)
  us.hetD<-quantile(us$V3[us.ind],probs=c(.05,.95),na.rm=T,names=F)
  us.homD<-quantile(us$V2[us.ind],probs=c(.05,.95),na.rm=T,names=F)
  ds.hetD<-quantile(ds$V3[ds.ind],probs=c(.05,.95),na.rm=T,names=F)
  ds.homD<-quantile(ds$V2[ds.ind],probs=c(.05,.95),na.rm=T,names=F)
  #heterozygously beneficial
  ua.hetBen <- which(ua$V3>ua.hetD[2])
  #heterozygously deleterious
  ua.hetDel <- which(ua$V3<ua.hetD[1])
  #homozygously beneficial
  ua.homBen <- which(ua$V2>ua.homD[2])
  #homozygously deleterious
  ua.homDel <- which(ua$V2<ua.homD[1])
  
  da.hetBen <- which(da$V3>da.hetD[2])
  da.hetDel <- which(da$V3<da.hetD[1])
  da.homBen <- which(da$V2>da.homD[2])
  da.homDel <- which(da$V2<da.homD[1])
  
  us.hetBen <- which(us$V3>us.hetD[2])
  us.hetDel <- which(us$V3<us.hetD[1])
  us.homBen <- which(us$V2>us.homD[2])
  us.homDel <- which(us$V2<us.homD[1])
  ds.hetBen <- which(ds$V3>ds.hetD[2])
  ds.hetDel <- which(ds$V3<ds.hetD[1])
  ds.homBen <- which(ds$V2>ds.homD[2])
  ds.homDel <- which(ds$V2<ds.homD[1])
  
  ua.eo <- intersect(ua.hetBen,ua.homDel)
  da.eo <- intersect(da.hetBen,da.homDel)
  us.eo <- intersect(us.hetBen,us.homDel)
  ds.eo <- intersect(ds.hetBen,ds.homDel)
  
  # for each het Ben find what evidence we have
  # how many barcodes and how many independent deletions
  # if it has >2 barcodes and all of them agree consider it hetBen
  
  #make list of all ORFs across all barcodes and for each find in how many barcodes it is
  allORFs <- unique(c(ua$V1,da$V1,us$V1,ds$V1))
  uaBool <- allORFs %in% ua$V1
  daBool <- allORFs %in% da$V1
  usBool <- allORFs %in% us$V1
  dsBool <- allORFs %in% ds$V1
  barcodeRepr <- uaBool+daBool+usBool+dsBool #how many barcodes we have for each deletion
  #for each ORF find if hetBen for each barcode, keep only those that are for all barcodes > 2
  #for each ORF in allORFs:
  barcodehetBen <- allORFs %in% ua[ua.hetBen,1] + allORFs %in% us[us.hetBen,1] + allORFs %in% da[da.hetBen,1] + allORFs %in% ds[ds.hetBen,1]
  barcodehomBen <- allORFs %in% ua[ua.homBen,1] + allORFs %in% us[us.homBen,1] + allORFs %in% da[da.homBen,1] + allORFs %in% ds[ds.homBen,1]
  barcodeeo <- allORFs %in% ua[ua.eo,1] + allORFs %in% da[da.eo,1] + allORFs %in% us[us.eo,1] + allORFs %in% ds[ds.eo,1]
  hetBen <- allORFs[barcodehetBen>minBarNo & barcodehetBen/barcodeRepr ==1]
  homBen <- allORFs[barcodehomBen>minBarNo & barcodehomBen/barcodeRepr ==1]
  eo <- allORFs[barcodeeo>minBarNo & barcodeeo/barcodeRepr == 1]
  YPDdel.hetBen <- which(hetBen %in% YPDdel[,1])
  YPDben.hetBen <- which(hetBen %in% YPDben[,1])
  YPDdel.eo <- which(eo %in% YPDdel[,1])
  YPDben.eo <- which(eo %in% YPDben[,1])
  #filter out MDR ORFs
  hetBen <- setdiff(hetBen,therio)
  eo <- setdiff(eo,therio)
  #filter out YPD beneficial or deleterious mutations
  hetBen <- setdiff(hetBen,YPDdel.hetBen)
  hetBen <- setdiff(hetBen,YPDben.hetBen)
  eo <- setdiff(eo,YPDdel.eo)
  eo <- setdiff(eo,YPDben.eo)
  
  eoL <- append(eoL,length(eo))
  hetBenL <- append(hetBenL,length(hetBen))
  totalL <- append(totalL,length(allORFs))
  
  write.table(hetBen,file=paste("~/projects/yeastOverdominance/collection/analysis/odORFs/",i,"hetBen.dat",sep=""),row.names=F,col.names=F)
  write.table(eo,file=paste("~/projects/yeastOverdominance/collection/analysis/odORFs/",i,"eo.dat",sep=""),row.names=F,col.names=F)
  print(i)
  counter <- counter + 1
}

save(eoL,file="~/projects/yeastOverdominance/collection/data/eoL.dat")
save(hetBenL,file="~/projects/yeastOverdominance/collection/data/hetBenL.dat")
#mean extreme overdominance
mean(eoL/hetBenL)*100
# [1] 5.670341

