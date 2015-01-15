#MDRfilter.Hil.R
# calculate MDR distribution for dataset I

#load the required libraries
library(stringr)
library(binom)

#and the datasets:
# - conditions: list of all conditions
# - dubious   : list of 120 dubious ORFS
# - allDubious: all dubious ORFs
# - YPDDel    : deleteriousin rich conditions
# - YPDBen    : beneficial in rich conditions

conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")
dubious<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)

allhomBen <- numeric(0) #homozygote Beneficial across all conditions
allhetBen <- numeric(0)
conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")
#loop through conditions, first pass
for(i in conditions$V1){
    ua <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ua.dat",sep=""),header=F)
    da <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".da.dat",sep=""),header=F)
    us <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".us.dat",sep=""),header=F)
    ds <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ds.dat",sep=""),header=F)
    ua$V1 <- str_extract(ua$V1,"[^:]*")
    us$V1 <- str_extract(us$V1,"[^:]*")
    da$V1 <- str_extract(da$V1,"[^:]*")
    ds$V1 <- str_extract(ds$V1,"[^:]*") 
    ua.ind<-which(ua$V1%in%dubious$dubiousORFS.Whitlock)
    da.ind<-which(da$V1%in%dubious$dubiousORFS.Whitlock)
    us.ind<-which(us$V1%in%dubious$dubiousORFS.Whitlock)
    ds.ind<-which(ds$V1%in%dubious$dubiousORFS.Whitlock)

    ua.hetD<-quantile(ua$V3[ua.ind],probs=c(.05,.95),na.rm=T,names=F)
    da.hetD<-quantile(da$V3[da.ind],probs=c(.05,.95),na.rm=T,names=F)
    us.hetD<-quantile(us$V3[us.ind],probs=c(.05,.95),na.rm=T,names=F)
    ds.hetD<-quantile(ds$V3[ds.ind],probs=c(.05,.95),na.rm=T,names=F)

    ua.hetBen <- which(ua$V3>ua.hetD[2]) #heterozygously beneficial
    da.hetBen <- which(da$V3>da.hetD[2])
    us.hetBen <- which(us$V3>us.hetD[2])
    ds.hetBen <- which(ds$V3>ds.hetD[2])
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
    
    hetBen <- allORFs[barcodehetBen>1 & barcodehetBen/barcodeRepr ==1]
    
    allhetBen <- append(allhetBen,hetBen)
    print(i)
}


# calculate MDR distribution
h.Hil <- hist(table(allhetBen),breaks=c(0:60),plot=FALSE)

#estimate Poisson distribution with the same mean
mdr.Hil <- table(allhetBen)
bootPoisson.Hil <- data.frame()
for (r in c(1:1000)){
    hb <- hist(rpois(length(mdr.Hil),mean(mdr.Hil)),breaks=c(0:60),plot=FALSE)
    bootPoisson.Hil <- rbind(bootPoisson.Hil,hb$counts)
}

# calculate not Trusted
notTrusted.Hil <- names(mdr.Hil)[mdr.Hil>10]

# save gene names of not Trusted
write.table(notTrusted.Hil,file="~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat",quote=FALSE, row.names=FALSE, col.names=FALSE)

# make MDR plot
setEPS()
postscript("~/projects/yeastOverdominance/collection/report/figures/Hil.MDR.eps")
par(las=1,bty="l")
heights <- numeric(0)
ul <- numeric(0)
ll <- numeric(0)
for(i in c(1:60)){
    qboot <- quantile(bootPoisson.Hil[,i],probs=c(0.25,0.5,0.75),names=F)
    heights <- append(heights,qboot[2])
    ul <- append(ul,qboot[3])
    ll <- append(ll,qboot[1])
}
logHeights <- log10(heights)
logHeights[which(!is.finite(logHeights))]<-0
plot(h.Hil$mids+0.5,log10(h.Hil$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,45)) #set axes etc
rect(h.Hil$mids,0,h.Hil$mids+1,logHeights,col="grey")
arrows(h.Hil$mids+.5,log10(ll),h.Hil$mids+.5,log10(ul),code=0,angle=90)
points(h.Hil$mids+0.5,log10(h.Hil$counts),pch=19)
dev.off()
