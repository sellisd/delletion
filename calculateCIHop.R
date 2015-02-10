hip<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HIP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5912)))
hop<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HOP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5846)))
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs
d<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
dI<-which(hop[,1]%in%d[,1]) # 120 neutrals
n<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/compoundNames.csv") 
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)
m<-read.table("~/projects/yeastOverdominance/collection/Hoepfner/groups.dat")
unq <- unique(paste(m$V1,m$V2)) #unique combinations of compound and concentration
load("~/projects/yeastOverdominance/collection/data/nss.dat")
load("~/projects/yeastOverdominance/collection/data/dss.dat")
load("~/projects/yeastOverdominance/collection/data/rss.dat")

notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)


totalL<-length(unq)
nsum <- matrix(nrow=totalL,ncol=5)
dsum <- matrix(nrow=totalL,ncol=5)
rsum <- matrix(nrow=totalL,ncol=5)

counter <- 1
condNames <- character()
condID <- character()
concs <- numeric(0)

# which gene deletions are not to be trusted
filterOut <- which(hop[,1] %in% therio)

for(i in unq){ #for each combination.
  cmbid <- strsplit(as.character(i)," ")[[1]][1]
  condID <- append(condID,cmbid)
  concentration <- strsplit(as.character(i)," ")[[1]][2]
  name<-n[which(cmbid==n[,1]),2]
  if(length(name)==0){
    titleString<-paste("CMB",cmbid,sep="")
  }else{
    titleString<-name
  }
  print(paste(counter,"/",length(unq),titleString))
  row.index <- which(paste(m$V1,m$V2)==i)
  hipi  <- m[row.index,4] == "HIP" & m[row.index,3] == 0
  hipZi <- m[row.index,4] == "HIP" & m[row.index,3] == 1
  hopi  <- m[row.index,4] == "HOP" & m[row.index,3] == 0
  hopZi <- m[row.index,4] == "HOP" & m[row.index,3] == 1
  
  hipIndex <- m[row.index,][hipi,5]
  hopIndex <- m[row.index,][hopi,5]
  hipIndexZ <- m[row.index,][hipZi,5]
  hopIndexZ <- m[row.index,][hopZi,5]
  
  if(length(hipIndex)>1){     # calculate averages for each column if more than one
    hipColumn <- rowMeans(hip[,hipIndex],na.rm=T)
  }else{
    hipColumn <- hip[,hipIndex]
  }
  if(length(hopIndex)>1){
    hopColumn <- rowMeans(hop[,hopIndex],na.rm=T)
  }else{
    hopColumn <- hop[,hopIndex]
  }
  if(length(hipIndexZ)>1){
    hipColumnZ <- rowMeans(hip[,hipIndexZ],na.rm=T)
  }else{
    hipColumnZ <- hip[,hipIndexZ]
  }
  if(length(hopIndexZ)>1){
    hopColumnZ <- rowMeans(hop[,hopIndexZ],na.rm=T)
  }else{
    hopColumnZ <- hop[,hopIndexZ]
  }
  #nss[hip,hipZ, hop,hopZ]
  ## 1 hetD lower
  ## 2 hetD upper
  ## 3 hetDZ lower
  ## 4 hetDZ upper
  ## 5 homD lower
  ## 6 homD upper
  ## 7 homDZ lower
  ## 8 homDZ upper
  
  #calculate proportion of HetBen, HetDel, HomBen and HetDel
  HetBen <- which(hipColumn>nss[counter,2])
  HetBenZ <- which(hipColumnZ>nss[counter,4])
  HetBeneficial <- intersect(HetBen,HetBenZ)
  HetBeneficial <- setdiff(HetBeneficial,filterOut)
  
  HetDel <- which(hipColumn<nss[counter,1])
  HetDelZ <- which(hipColumnZ<nss[counter,3])
  HetDeleterious <- intersect(HetDel,HetDelZ)
  HetDeleterious <- setdiff(HetDeleterious,filterOut)
  
  HomBen <- which(hopColumn>nss[counter,6])
  HomBenZ <- which(hopColumnZ>nss[counter,8])
  HomBeneficial <- intersect(HomBen,HomBenZ)
  HomBeneficial <- setdiff(HomBeneficial,filterOut)
  
  HomDel <- which(hopColumn<nss[counter,5])
  HomDelZ <- which(hopColumnZ<nss[counter,7])
  HomDeleterious <- intersect(HomDel,HomDelZ)
  HomDeleterious <- setdiff(HomDeleterious,filterOut)
  
  #remove from YPD del and YPD ben
  YPDdel.hetBen <- which(HetBeneficial %in% YPDdel[,1])
  YPDben.hetBen <- which(HetBeneficial %in% YPDben[,1])
  HetBeneficial <- setdiff(HetBeneficial,YPDdel.hetBen)
  HetBeneficial <- setdiff(HetBeneficial,YPDben.hetBen)
  
  #extreme overdominant
  eo <- intersect(HetBeneficial,HomDeleterious)
  
  #save ORFs in file for HetBeneficial and eo
  
  fileName <- paste(titleString,".",concentration,sep="")
  fileName <- gsub("/","_",fileName)
  condNames <- append(condNames,as.character(titleString))
  concs <- append(concs,concentration)
  
  write.table(hip[HetBeneficial,1],file=paste("~/projects/yeastOverdominance/collection/Hoepfner/odORFs/",fileName,"hetBen.dat",sep=""),row.names=F,col.names=F)
  write.table(hip[eo,1],file=paste("~/projects/yeastOverdominance/collection/Hoepfner/odORFs/",fileName,"eo.dat",sep=""),row.names=F,col.names=F)
  
  #save
  nsum[counter,] <- c(length(HetBeneficial),length(HetDeleterious),length(HomBeneficial),length(HomDeleterious),length(eo))
  
  #repeat for dubious and random
  dHetBen <- which(hipColumn>dss[counter,2])
  dHetBenZ <- which(hipColumnZ>dss[counter,4])
  dHetBeneficial <- intersect(dHetBen,dHetBenZ)
  dHetBeneficial <- setdiff(dHetBeneficial,filterOut)   
  dHetDel <- which(hipColumn<dss[counter,1])
  dHetDelZ <- which(hipColumnZ<dss[counter,3])
  dHetDeleterious <- intersect(dHetDel,dHetDelZ)
  dHetDeleterious <- setdiff(dHetDeleterious,filterOut)
  dHomBen <- which(hopColumn>dss[counter,6])
  dHomBenZ <- which(hopColumnZ>dss[counter,8])
  dHomBeneficial <- intersect(dHomBen,dHomBenZ)
  dHomBeneficial <- setdiff(dHomBeneficial,filterOut)
  dHomDel <- which(hopColumn<dss[counter,5])
  dHomDelZ <- which(hopColumnZ<dss[counter,7])
  dHomDeleterious <- intersect(dHomDel,dHomDelZ)
  dHomDeleterious <- setdiff(dHomDeleterious,filterOut)
  deo <- intersect(dHetBeneficial,dHomDeleterious)
  dsum[counter,] <- c(length(dHetBeneficial),length(dHetDeleterious),length(dHomBeneficial),length(dHomDeleterious),length(deo))
  
  rHetBen <- which(hipColumn>rss[counter,2])
  rHetBenZ <- which(hipColumnZ>rss[counter,4])
  rHetBeneficial <- intersect(rHetBen,rHetBenZ)
  rHetBeneficial <- setdiff(rHetBeneficial,filterOut)   
  rHetDel <- which(hipColumn<rss[counter,1])
  rHetDelZ <- which(hipColumnZ<rss[counter,3])
  rHetDeleterious <- intersect(rHetDel,rHetDelZ)
  rHetDeleterious <- setdiff(rHetDeleterious,filterOut)
  rHomBen <- which(hopColumn>rss[counter,6])
  rHomBenZ <- which(hopColumnZ>rss[counter,8])
  rHomBeneficial <- intersect(rHomBen,rHomBenZ)
  rHomBeneficial <- setdiff(rHomBeneficial,filterOut)
  rHomDel <- which(hopColumn<rss[counter,5])
  rHomDelZ <- which(hopColumnZ<rss[counter,7])
  rHomDeleterious <- intersect(rHomDel,rHomDelZ)
  rHomDeleterious <- setdiff(rHomDeleterious,filterOut)
  reo <- intersect(rHetBeneficial,rHomDeleterious)
  rsum[counter,] <- c(length(rHetBeneficial),length(rHetDeleterious),length(rHomBeneficial),length(rHomDeleterious),length(reo))
  counter <- counter + 1
}


save(nsum,file="~/projects/yeastOverdominance/collection/data/nsum.dat")
save(dsum,file="~/projects/yeastOverdominance/collection/data/dsum.dat")
save(rsum,file="~/projects/yeastOverdominance/collection/data/rsum.dat")
save(condNames,file="~/projects/yeastOverdominance/collection/data/condNames.dat")
save(concs,file="~/projects/yeastOverdominance/collection/data/concs.dat")
save(condID,file="~/projects/yeastOverdominance/collection/data/condID.dat")



