# jackknife resampling of Dataset II

maxJ <- 1000 # number of jackknife resampling

# read required libraries

# read required data files
load("~/projects/yeastOverdominance/collection/data/nss.dat")   #nss
hip<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HIP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5912)))
hop<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HOP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5846)))
d<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
n<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/compoundNames.csv")
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)
m<-read.table("~/projects/yeastOverdominance/collection/Hoepfner/groups.dat")
unq <- unique(paste(m$V1,m$V2)) #unique combinations of compound and concentrat

# load not trusted genes
notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)
neutral <- d[which(!d[,1]%in%therio),1] # neutral set without not trusted
dI<-which(hop[,1]%in%neutral) # neutral set without not trusted

# which gene deletions are not to be trusted
filterOut <- which(hop[,1] %in% therio)
allORF <- hop[,1]
allORF <- which(!allORF %in% therio) # exclude
counter <- 1
Weo <- numeric(0)
Whb <- numeric(0)
for(i in unq){ #for each combination of compound and condition
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
  eoLN <- numeric(0)
  hetBenLN <- numeric(0)
  eoLR <- numeric(0)
  hetBenLR <- numeric(0)
  for (j in c(1:maxJ)){ #   for maxJ
    half <- sample(dI,60)  # subsample 55 from neutral set
    rest <- setdiff(dI,half) # the rest 60
    rrest <- sample(setdiff(allORF,half),60) #random other (but not the excluded ones)
    # calculate hetBen, hetBenZ, homDel, homDelZ
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
    # how many are in half, rest and rrest
    eoLN <- append(eoLN,sum(rest%in%eo))
    eoLR <- append(eoLR,sum(rrest%in%eo))
    hetBenLN <- append(hetBenLN,sum(rest%in%HetBeneficial))
    hetBenLR <- append(hetBenLR,sum(rrest%in%HetBeneficial))
  }
  Weo <- append(Weo,wilcox.test(eoLN,eoLR,alternative="less")$p.value)
  Whb <- append(Whb,wilcox.test(hetBenLN,hetBenLR,alternative="less")$p.value)
  print(paste(counter, i))
  counter <- counter + 1
}

if(0){
    plot(hipColumn,hipColumnZ,pch=19,cex=0.3,col="#00000050",xlim=c(-6,2))
    points(hipColumn[half],hipColumnZ[half],pch=19,cex=0.5,col="#ff000050")
    points(hipColumn[rest],hipColumnZ[rest],pch=19,cex=0.5,col="#0000ff50")
    points(hipColumn[rrest],hipColumnZ[rrest],pch=19,cex=0.5,col="#00ff0050")
    abline(v=nss[counter,1:2])
    abline(h=nss[counter,3:4])
    points(hipColumn[eo],hipColumnZ[eo],pch=1,col="#aaaa00")
}

# save
write.table(data.frame(Weo,Whb),file="~/projects/yeastOverdominance/collection/data/jk.Hop.dat",quote=FALSE,row.names=FALSE, col.names=TRUE)
