#MDRfilter.Hop.R
# calculate MDR distribution for dataset II

#read necessary libraries
library(stringr)
#read data files
hip<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HIP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5912)))
hop<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HOP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5846)))
m<-read.table("~/projects/yeastOverdominance/collection/Hoepfner/groups.dat")
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs
d<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
dI<-which(hop[,1]%in%d[,1]) # 120 neutrals
n<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/compoundNames.csv") 
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)
#m$V1 compound
#m$V2 concentration
#m$V3 z-score?
#m$V4 HIP/HOP
#m$V5 column index
unq <- unique(paste(m$V1,m$V2)) #unique combinations of compound and concentration
filteredORFs<-rep(0,length(hip[,1])) # list of counters, one for each gene, to count how many times is heterozygously beneficail
counter <- 0
for(i in unq){ #for each combination.
    print(paste(counter,"/",length(unq)))
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

    #calculate percentiles
    hipD<-quantile(hipColumn[dI],probs=c(.05,.95),na.rm=T,names=F)
    hopD<-quantile(hopColumn[dI],probs=c(.05,.95),na.rm=T,names=F)
    hipDZ<-quantile(hipColumnZ[dI],probs=c(.05,.95),na.rm=T,names=F)
    hopDZ<-quantile(hopColumnZ[dI],probs=c(.05,.95),na.rm=T,names=F)
    #find outliers for HIP experiments
    HetBen<-which(hipColumn>hipD[2]) # Row Index of strains with HIP larger than 95% 
    HetBenZ<-which(hipColumnZ>hipDZ[2]) # Row Index of strains with HIP larger than 95% for z-score
    #intersection
    HetResistant <- intersect(HetBen,HetBenZ)
    #find multidrug resistant ORFs
    filteredORFs[HetResistant] <- filteredORFs[HetResistant] + 1
    counter <- counter + 1
}


#calculate MDR distribution
h.Hop <- hist(filteredORFs,breaks=c(0:600),plot=F)
#estimate Poisson distribution with the same mean
bootPoisson.Hop <- matrix(nrow=1000,ncol=600)
for (r in c(1:1000)){
    hb <- hist(rpois(length(filteredORFs[filteredORFs!=0]),mean(filteredORFs[filteredORFs!=0])),breaks=c(0:600),plot=F)
    bootPoisson.Hop[r,] <- hb$counts
}

# calculate not Trusted
notTrusted.Hop <- hop[which(filteredORFs>200),1]

#save gene names of not Trusted
write.table(notTrusted.Hop,file="~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat",quote=FALSE, row.names=FALSE, col.names=FALSE)

# make MDR plot
setEPS()
postscript("~/projects/yeastOverdominance/collection/report/figures/Hop.MDR.eps")
par(las=1,bty="l")
heights <- numeric(0)
ul <- numeric(0)
ll <- numeric(0)
for(i in c(1:600)){
  qboot <- quantile(bootPoisson.Hop[,i],probs=c(0.25,0.5,0.75),names=F)
  heights <- append(heights,qboot[2])
  ul <- append(ul,qboot[3])
  ll <- append(ll,qboot[1])
}
logHeights <- log10(heights)
logHeights[which(!is.finite(logHeights))]<-0
plot(h.Hop$mids+0.5,log10(h.Hop$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,450)) #set axes etc
rect(h.Hop$mids,0,h.Hop$mids+1,logHeights,col="grey",border="grey")
arrows(h.Hop$mids+.5,log10(ll),h.Hop$mids+.5,log10(ul),code=0,angle=90)
points(h.Hop$mids+0.5,log10(h.Hop$counts),pch=19)
dev.off()
