## Deletion collection analysis
## -----------------------------
##     analysis of deletion collection data. Separately for Hillenmeyer etal 2008 and Hoepfner et al., 2014 but with similar methodology

 
#read required libraries
library(stringr)
library(binom)
#read required files
conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")
dubious<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)

# first pass: find MDR
#---------------------
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

h <- hist(table(allhetBen),breaks=c(0:60),plot=FALSE)
                                        #estimate Poisson distribution with the same mean
mdr <- table(allhetBen)
bootPoisson <- data.frame()
for (r in c(1:1000)){
    hb <- hist(rpois(length(mdr),mean(mdr)),breaks=c(0:60),plot=FALSE)
    bootPoisson <- rbind(bootPoisson,hb$counts)
}
notTrusted <- names(mdr)[mdr>10] #Hil

hopnt<-read.table("~/projects/yeastOverdominance/collection/Hoepfner/HopNotTrusted.dat")
hilnt<-read.table("~/projects/yeastOverdominance/collection/analysis/HilNotTrusted.dat")
therio <- union(hopnt$V1,hilnt$V1)
## why therio?
## > length(therio)
## [1] 666

notTrusted<-therio
                                        #each row is a bootstrap. for each integer from min to max get mean and quartiles

pdf("~/projects/dissertation/thesis/img/chapter3/multidrugHillenmeyer.pdf")
par(las=1,bty="l")
heights <- numeric(0)
ul <- numeric(0)
ll <- numeric(0)
for(i in c(1:60)){
    qboot <- quantile(bootPoisson[,i],probs=c(0.25,0.5,0.75),names=F)
    heights <- append(heights,qboot[2])
    ul <- append(ul,qboot[3])
    ll <- append(ll,qboot[1])
}
logHeights <- log10(heights)
logHeights[which(!is.finite(logHeights))]<-0
plot(h$mids+0.5,log10(h$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,45)) #set axes etc
rect(h$mids,0,h$mids+1,logHeights,col="grey")
arrows(h$mids+.5,log10(ll),h$mids+.5,log10(ul),code=0,angle=90)
points(h$mids+0.5,log10(h$counts),pch=19)
dev.off()

#plot for both datasets
#-----------------------
pdf("~/projects/dissertation/thesis/img/chapter3/multidrug.pdf",width=10,height=6)
par(las=1,bty="l",mfrow=c(1,2))
heights <- numeric(0)
ul <- numeric(0)
ll <- numeric(0)
for(i in c(1:60)){
    qboot <- quantile(bootPoisson[,i],probs=c(0.25,0.5,0.75),names=F)
    heights <- append(heights,qboot[2])
    ul <- append(ul,qboot[3])
    ll <- append(ll,qboot[1])
}
logHeights <- log10(heights)
logHeights[which(!is.finite(logHeights))]<-0
plot(h$mids+0.5,log10(h$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,45)) #set axes etc
rect(h$mids,0,h$mids+1,logHeights,col="grey")
arrows(h$mids+.5,log10(ll),h$mids+.5,log10(ul),code=0,angle=90)
points(h$mids+0.5,log10(h$counts),pch=19)
#-----------
load("~/projects/yeastOverdominance/collection/Hoepfner/multidrugHop.h.dat")
load("~/projects/yeastOverdominance/collection/Hoepfner/multidrugHop.filteredORFs.dat")
bootPoisson <- matrix(nrow=1000,ncol=600)
for (r in c(1:1000)){
    hb <- hist(rpois(length(filteredORFs[filteredORFs!=0]),mean(filteredORFs[filteredORFs!=0])),breaks=c(0:600),plot=F)
    bootPoisson[r,] <- hb$counts
}
heights <- numeric(0)
ul <- numeric(0)
ll <- numeric(0)
for(i in c(1:600)){
  qboot <- quantile(bootPoisson[,i],probs=c(0.25,0.5,0.75),names=F)
  heights <- append(heights,qboot[2])
  ul <- append(ul,qboot[3])
  ll <- append(ll,qboot[1])
}
logHeights <- log10(heights)
logHeights[which(!is.finite(logHeights))]<-0
plot(h$mids+0.5,log10(h$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,450)) #set axes etc
rect(h$mids,0,h$mids+1,logHeights,col="grey",border="grey")
arrows(h$mids+.5,log10(ll),h$mids+.5,log10(ul),code=0,angle=90)
points(h$mids+0.5,log10(h$counts),pch=19)
dev.off()


notTrusted <- names(mdr)[mdr>10]
write.table(notTrusted,file="~/projects/yeastOverdominance/collection/analysis/HilNotTrusted.dat",quote=F,row.names=F,col.names=F)

# second pass: dubious and random sets
#-------------------------------------
#explore set of 120 and calculate CI for random and dubious sets

ruahet <-numeric(0)
rdahet <-numeric(0)
rushet <-numeric(0)
rdshet <-numeric(0)
         
ruahom <-numeric(0)
rdahom <-numeric(0)
rushom <-numeric(0)
rdshom <-numeric(0)
         
duahet <-numeric(0)
ddahet <-numeric(0)
dushet <-numeric(0)
ddshet <-numeric(0)
         
duahom <-numeric(0)
ddahom <-numeric(0)
dushom <-numeric(0)
ddshom <-numeric(0)
         
uahet <- numeric(0)
dahet <- numeric(0)
ushet <- numeric(0)
dshet <- numeric(0)
         
uahom <- numeric(0)
dahom <- numeric(0)
ushom <- numeric(0)
dshom <- numeric(0)


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

                                        # confidence intervals from the manually curated set of 120 dubious ORFs
    ua.ind <- which(ua$V1%in%dubious$dubiousORFS.Whitlock)
    da.ind <- which(da$V1%in%dubious$dubiousORFS.Whitlock)
    us.ind <- which(us$V1%in%dubious$dubiousORFS.Whitlock)
    ds.ind <- which(ds$V1%in%dubious$dubiousORFS.Whitlock)

    ua.hetD <- quantile(ua$V3[ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    ua.homD <- quantile(ua$V2[ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    da.hetD <- quantile(da$V3[da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    da.homD <- quantile(da$V2[da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    us.hetD <- quantile(us$V3[us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    us.homD <- quantile(us$V2[us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    ds.hetD <- quantile(ds$V3[ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
    ds.homD <- quantile(ds$V2[ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)

                                        #dubious dubious ORFs
    dds.ua.hetD <- numeric(0)
    dds.ua.homD <- numeric(0)
    dds.da.hetD <- numeric(0)
    dds.da.homD <- numeric(0)
    dds.us.hetD <- numeric(0)
    dds.us.homD <- numeric(0)
    dds.ds.hetD <- numeric(0)
    dds.ds.homD <- numeric(0)

    rds.ua.hetD <- numeric(0)
    rds.ua.homD <- numeric(0)
    rds.da.hetD <- numeric(0)
    rds.da.homD <- numeric(0)
    rds.us.hetD <- numeric(0)
    rds.us.homD <- numeric(0)
    rds.ds.hetD <- numeric(0)
    rds.ds.homD <- numeric(0)

    allORFs<-unique(c(ua$V1,da$V1,us$V1,ds$V1))

    for(boot in c(1:1000)){
        dds <- sample(allDubious$V2,120) #dubious dubious ORFs

        dd.ua.ind <- which(ua$V1%in%dds)
        dd.da.ind <- which(da$V1%in%dds)
        dd.us.ind <- which(us$V1%in%dds)
        dd.ds.ind <- which(ds$V1%in%dds)

        rds <- sample(allORFs,120)  #random dubious ORFs out of all possible

        rdd.ua.ind <- which(ua$V1%in%rds)
        rdd.da.ind <- which(da$V1%in%rds)
        rdd.us.ind <- which(us$V1%in%rds)
        rdd.ds.ind <- which(ds$V1%in%rds)

                                        #calculate quantiles
        dd.ua.hetD <- quantile(ua$V3[dd.ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.ua.homD <- quantile(ua$V2[dd.ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.da.hetD <- quantile(da$V3[dd.da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.da.homD <- quantile(da$V2[dd.da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.us.hetD <- quantile(us$V3[dd.us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.us.homD <- quantile(us$V2[dd.us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.ds.hetD <- quantile(ds$V3[dd.ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        dd.ds.homD <- quantile(ds$V2[dd.ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)

        rd.ua.hetD <- quantile(ua$V3[rdd.ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.ua.homD <- quantile(ua$V2[rdd.ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.da.hetD <- quantile(da$V3[rdd.da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.da.homD <- quantile(da$V2[rdd.da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.us.hetD <- quantile(us$V3[rdd.us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.us.homD <- quantile(us$V2[rdd.us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.ds.hetD <- quantile(ds$V3[rdd.ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        rd.ds.homD <- quantile(ds$V2[rdd.ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)

                                        #save across replicates
        dds.ua.hetD <- rbind(dds.ua.hetD,dd.ua.hetD)
        dds.ua.homD <- rbind(dds.ua.homD,dd.ua.homD)
        dds.da.hetD <- rbind(dds.da.hetD,dd.da.hetD)
        dds.da.homD <- rbind(dds.da.homD,dd.da.homD)
        dds.us.hetD <- rbind(dds.us.hetD,dd.us.hetD)
        dds.us.homD <- rbind(dds.us.homD,dd.us.homD)
        dds.ds.hetD <- rbind(dds.ds.hetD,dd.ds.hetD)
        dds.ds.homD <- rbind(dds.ds.homD,dd.ds.homD)

        rds.ua.hetD <- rbind(rds.ua.hetD,rd.ua.hetD)
        rds.ua.homD <- rbind(rds.ua.homD,rd.ua.homD)
        rds.da.hetD <- rbind(rds.da.hetD,rd.da.hetD)
        rds.da.homD <- rbind(rds.da.homD,rd.da.homD)
        rds.us.hetD <- rbind(rds.us.hetD,rd.us.hetD)
        rds.us.homD <- rbind(rds.us.homD,rd.us.homD)
        rds.ds.hetD <- rbind(rds.ds.hetD,rd.ds.hetD)
        rds.ds.homD <- rbind(rds.ds.homD,rd.ds.homD)
    }

                                        #save list of CI for random and dubious subsets
    ruahet <- rbind(ruahet, colMeans(rds.ua.hetD))
    rdahet <- rbind(rdahet, colMeans(rds.da.hetD))
    rushet <- rbind(rushet, colMeans(rds.us.hetD))
    rdshet <- rbind(rdshet, colMeans(rds.ds.hetD))
    
    ruahom <- rbind(ruahom, colMeans(rds.ua.homD))
    rdahom <- rbind(rdahom, colMeans(rds.da.homD))
    rushom <- rbind(rushom, colMeans(rds.us.homD))
    rdshom <- rbind(rdshom, colMeans(rds.ds.homD))
    
    duahet <- rbind(duahet, colMeans(dds.ua.hetD))
    ddahet <- rbind(ddahet, colMeans(dds.da.hetD))
    dushet <- rbind(dushet, colMeans(dds.us.hetD))
    ddshet <- rbind(ddshet, colMeans(dds.ds.hetD))
    
    duahom <- rbind(duahom, colMeans(dds.ua.homD))
    ddahom <- rbind(ddahom, colMeans(dds.da.homD))
    dushom <- rbind(dushom, colMeans(dds.us.homD))
    ddshom <- rbind(ddshom, colMeans(dds.ds.homD))

    uahet <- rbind(uahet, ua.hetD)
    dahet <- rbind(dahet, da.hetD)
    ushet <- rbind(ushet, us.hetD)
    dshet <- rbind(dshet, ds.hetD)
    
    uahom <- rbind(uahom, ua.homD)
    dahom <- rbind(dahom, da.homD)
    ushom <- rbind(ushom, us.homD)
    dshom <- rbind(dshom, ds.homD)
    print(i)
}

                                        #compare upper CI
          #   neutral        dubious      random
muahom <- data.frame(uahom[,3], duahom[,3],ruahom[,3])
mdahom <- data.frame(dahom[,3], ddahom[,3],rdahom[,3])
mushom <- data.frame(ushom[,3], dushom[,3],rushom[,3])
mdshom <- data.frame(dshom[,3], ddshom[,3],rdshom[,3])

muahet <- data.frame(uahet[,3], duahet[,3],ruahet[,3])
mdahet <- data.frame(dahet[,3], ddahet[,3],rdahet[,3])
mushet <- data.frame(ushet[,3], dushet[,3],rushet[,3])
mdshet <- data.frame(dshet[,3], ddshet[,3],rdshet[,3])

#upper CI
#neutral > dubious
sum(muahom[,1]-muahom[,2]>0)
sum(muahom[,1]-muahom[,2]>0)
sum(muahom[,1]-muahom[,2]>0)
sum(muahom[,1]-muahom[,2]>0)
#neutral > random
sum(muahom[,1]-muahom[,3]>0)
sum(muahom[,1]-muahom[,3]>0)
sum(muahom[,1]-muahom[,3]>0)
sum(muahom[,1]-muahom[,3]>0)
#dubious > random
sum(muahom[,2]-muahom[,3]>0)
sum(muahom[,2]-muahom[,3]>0)
sum(muahom[,2]-muahom[,3]>0)
sum(muahom[,2]-muahom[,3]>0)


#compare ranges for homozygotes and heterozygotes for each barcode
uadfhet <- data.frame(uahet[,3] - uahet[,1], duahet[,3]-duahet[,1],ruahet[,3]-ruahet[,1])
dadfhet <- data.frame(dahet[,3] - dahet[,1], ddahet[,3]-ddahet[,1],rdahet[,3]-rdahet[,1])
usdfhet <- data.frame(ushet[,3] - ushet[,1], dushet[,3]-dushet[,1],rushet[,3]-rushet[,1])
dsdfhet <- data.frame(dshet[,3] - dshet[,1], ddshet[,3]-ddshet[,1],rdshet[,3]-rdshet[,1])

uadfhom <- data.frame(uahom[,3] - uahom[,1], duahom[,3]-duahom[,1],ruahom[,3]-ruahom[,1])
dadfhom <- data.frame(dahom[,3] - dahom[,1], ddahom[,3]-ddahom[,1],rdahom[,3]-rdahom[,1])
usdfhom <- data.frame(ushom[,3] - ushom[,1], dushom[,3]-dushom[,1],rushom[,3]-rushom[,1])
dsdfhom <- data.frame(dshom[,3] - dshom[,1], ddshom[,3]-ddshom[,1],rdshom[,3]-rdshom[,1])

(sum(uadfhom[,1]-uadfhom[,2]<0)/111*100 +
sum(dadfhom[,1]-dadfhom[,2]<0)/111*100 +
sum(usdfhom[,1]-usdfhom[,2]<0)/111*100 +
sum(dsdfhom[,1]-dsdfhom[,2]<0)/111*100)/4 # Homozygotes: neutral ORFs set range - dubious ORFs set range
#[1] 98.64865

(sum(uadfhom[,1]-uadfhom[,3]<0)/111*100 +
sum(dadfhom[,1]-dadfhom[,3]<0)/111*100 +
sum(usdfhom[,1]-usdfhom[,3]<0)/111*100 +
sum(dsdfhom[,1]-dsdfhom[,3]<0)/111*100)/4 #Homozygotes: range - range Random
#[1] 97.74775

(sum(uadfhom[,2]-uadfhom[,3]<0)/111*100 +
sum(dadfhom[,2]-dadfhom[,3]<0)/111*100 +
sum(usdfhom[,2]-usdfhom[,3]<0)/111*100 +
sum(dsdfhom[,2]-dsdfhom[,3]<0)/111*100)/4 #Homozygotes: range dubious - range random
#[1] 75


(sum(uadfhet[,1]-uadfhet[,2]<0)/111*100 +
sum(dadfhet[,1]-dadfhet[,2]<0)/111*100 +
sum(usdfhet[,1]-usdfhet[,2]<0)/111*100 +
sum(dsdfhet[,1]-dsdfhet[,2]<0)/111*100)/4 #Heterozygotes
#[1] 59.23423

(sum(uadfhet[,1]-uadfhet[,3]<0)/111*100 +
sum(dadfhet[,1]-dadfhet[,3]<0)/111*100 +
sum(usdfhet[,1]-usdfhet[,3]<0)/111*100 +
sum(dsdfhet[,1]-dsdfhet[,3]<0)/111*100)/4
#[1] 56.75676

(sum(uadfhet[,2]-uadfhet[,3]<0)/111*100 +
sum(dadfhet[,2]-dadfhet[,3]<0)/111*100 +
sum(usdfhet[,2]-usdfhet[,3]<0)/111*100 +
sum(dsdfhet[,2]-dsdfhet[,3]<0)/111*100)/4
#[1] 49.54955

# third pass 
#-------------
minBarNo <- 1 # >minBarNo should be available for each gene 
eoL <- numeric(0) # number of extreme overdominant in each condition
hetBenL <- numeric(0)
eoD<- numeric(0) # number of extreme overdominant in each condition
hetBenD<- numeric(0)
eoR <- numeric(0) # number of extreme overdominant in each condition
hetBenR <- numeric(0)
totalL <- numeric(0)
counter <- 1

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
    hetBen <- setdiff(hetBen,notTrusted)
    eo <- setdiff(eo,notTrusted)

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

#repeat for dubious ORFs CI

                                        #heterozygously beneficial
    ua.hetBen <- which(ua$V3>duahet[counter,3])
                                        #heterozygously deleterious
    ua.hetDel <- which(ua$V3<duahet[counter,1])
                                        #homozygously beneficial
    ua.homBen <- which(ua$V2>duahom[counter,3])
                                        #homozygously deleterious
    ua.homDel <- which(ua$V2<duahom[counter,1])

    da.hetBen <- which(da$V3>ddahet[counter,3])
    da.hetDel <- which(da$V3<ddahet[counter,1])
    da.homBen <- which(da$V2>ddahom[counter,3])
    da.homDel <- which(da$V2<ddahom[counter,1])

    us.hetBen <- which(us$V3>dushet[counter,3])
    us.hetDel <- which(us$V3<dushet[counter,1])
    us.homBen <- which(us$V2>dushom[counter,3])
    us.homDel <- which(us$V2<dushom[counter,1])
    ds.hetBen <- which(ds$V3>ddshet[counter,3])
    ds.hetDel <- which(ds$V3<ddshet[counter,1])
    ds.homBen <- which(ds$V2>ddshom[counter,3])
    ds.homDel <- which(ds$V2<ddshom[counter,1])

    ua.eo <- intersect(ua.hetBen,ua.homDel)
    da.eo <- intersect(da.hetBen,da.homDel)
    us.eo <- intersect(us.hetBen,us.homDel)
    ds.eo <- intersect(ds.hetBen,ds.homDel)
    
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
    hetBen <- setdiff(hetBen,notTrusted)
    eo <- setdiff(eo,notTrusted)

                                        #filter out YPD beneficial or deleterious mutations
    hetBen <- setdiff(hetBen,YPDdel.hetBen)
    hetBen <- setdiff(hetBen,YPDben.hetBen)
    eo <- setdiff(eo,YPDdel.eo)
    eo <- setdiff(eo,YPDben.eo)

    eoD <- append(eoD,length(eo))
    hetBenD <- append(hetBenD,length(hetBen))

#repeat for random ORFs CI

                                        #heterozygously beneficial
    ua.hetBen <- which(ua$V3>ruahet[counter,3])
                                        #heterozygously deleterious
    ua.hetDel <- which(ua$V3<ruahet[counter,1])
                                        #homozygously beneficial
    ua.homBen <- which(ua$V2>ruahom[counter,3])
                                        #homozygously deleterious
    ua.homDel <- which(ua$V2<ruahom[counter,1])
    
    da.hetBen <- which(da$V3>rdahet[counter,3])
    da.hetDel <- which(da$V3<rdahet[counter,1])
    da.homBen <- which(da$V2>rdahom[counter,3])
    da.homDel <- which(da$V2<rdahom[counter,1])
    
    us.hetBen <- which(us$V3>rushet[counter,3])
    us.hetDel <- which(us$V3<rushet[counter,1])
    us.homBen <- which(us$V2>rushom[counter,3])
    us.homDel <- which(us$V2<rushom[counter,1])
    ds.hetBen <- which(ds$V3>rdshet[counter,3])
    ds.hetDel <- which(ds$V3<rdshet[counter,1])
    ds.homBen <- which(ds$V2>rdshom[counter,3])
    ds.homDel <- which(ds$V2<rdshom[counter,1])
    
    ua.eo <- intersect(ua.hetBen,ua.homDel)
    da.eo <- intersect(da.hetBen,da.homDel)
    us.eo <- intersect(us.hetBen,us.homDel)
    ds.eo <- intersect(ds.hetBen,ds.homDel)

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
    hetBen <- setdiff(hetBen,notTrusted)
    eo <- setdiff(eo,notTrusted)
    
                                        #filter out YPD beneficial or deleterious mutations
    hetBen <- setdiff(hetBen,YPDdel.hetBen)
    hetBen <- setdiff(hetBen,YPDben.hetBen)
    eo <- setdiff(eo,YPDdel.eo)
    eo <- setdiff(eo,YPDben.eo)

    eoR <- append(eoR,length(eo))
    hetBenR <- append(hetBenR,length(hetBen))

    print(i)
    counter <- counter + 1
}

#mean extreme overdominance
mean(eoL/hetBenL)*100
# [1] 5.670341
mean(eoD/hetBenD)*100
# [1] 3.415634
mean(eoR/hetBenR)*100
#[1] 2.628345

png("eoRatioHillenmeyer.png",height=800,width=640)
par(las=1,bty="l",mai=c(.5,3,0,0.2),mgp=c(1.5,0.5,0),cex=0.8)
hb <- data.frame(conditions$V1,eoL/hetBenL,eoL,hetBenL)
orderedhB <- hb[order(hb[,2]),]
orderedhB <- orderedhB[which(orderedhB[,2]>0),]
ci <- binom.confint(orderedhB[,3],orderedhB[,4],method="exact")
bar <- barplot(orderedhB[,2],names.arg=orderedhB[,1],horiz=TRUE,xlim=c(0,0.55),xlab="extreme overdominance")
arrows(ci$lower,bar,ci$upper,bar,code=0)
dev.off()
write.csv(hb,file="~/projects/yeastOverdominance/collection/analysis/odRatio.csv",row.names=F)

l <- length(hb[,1])
pdf("eoRatioHillenmeyer.pdf")
par(las=1,bty="l")
orderedHi <- hb[order(hb[,2]),]
ci <- binom.confint(orderedHi[,3],orderedHi[,4],method="exact")
plot(100*orderedHi[,2],ylab="extreme overdominance (%)",xlab="conditions",type="l",ylim=c(0,50))
polygon(c(c(1:l),c(l:1)),c(ci$lower,rev(ci$upper))*100,col="grey60",border=NA)
points(orderedHi[,2]*100,type="l",lwd=2)
dev.off()

#fourth pass
#-----------
## for each condition
## fot 1000 repeats
#   randomly split in half the neutral ORFs set
#   use the first half as a set of neutral ORFs and the other will be tested to see how different they are
#   pick 60 random gene deletion strains excluding the 60 of the first subset
#  find how many of the second set are considered extreme overdominant and how many of the random ones
#  compare the overlap of the 60 subset to the 60 random to extreme overdominance and hetben
#  we expect that the 60 random should more often be extreme overdominant and het ben than the 60 from the second subset
# this is true for XXX80% of eo and YYY29% of het ben
maxJ <- 1000
Weo <- numeric(0)
Whb <- numeric(0)
for(i in conditions$V1){
    ua <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ua.dat",sep=""),header=F)
    da <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".da.dat",sep=""),header=F)
    us <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".us.dat",sep=""),header=F)
    ds <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ds.dat",sep=""),header=F)
    ua$V1 <- str_extract(ua$V1,"[^:]*")
    us$V1 <- str_extract(us$V1,"[^:]*")
    da$V1 <- str_extract(da$V1,"[^:]*")
    ds$V1 <- str_extract(ds$V1,"[^:]*")
    allORF <- unique(intersect(intersect(ua$V1,us$V1),intersect(da$V1,ds$V1)))
    eoLN <- numeric(0)
    hetBenLN <- numeric(0)
    eoLR <- numeric(0)
    hetBenLR <- numeric(0)
    totalLN <- numeric(0)
    for (j in c(1:maxJ)){
        halph <- sample(dubious$dubiousORFS.Whitlock,60)  # subsample 60 from neutral set
        rest <- setdiff(dubious$dubiousORFS.Whitlock,halph) #the rest
        rrest <- sample(setdiff(allORF,halph),60) #random other
        ua.ind<-which(ua$V1%in%halph)
        da.ind<-which(da$V1%in%halph)
        us.ind<-which(us$V1%in%halph)
        ds.ind<-which(ds$V1%in%halph)

        ua.hetD <- quantile(ua$V3[ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        ua.homD <- quantile(ua$V2[ua.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        da.hetD <- quantile(da$V3[da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        da.homD <- quantile(da$V2[da.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        us.hetD <- quantile(us$V3[us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        us.homD <- quantile(us$V2[us.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        ds.hetD <- quantile(ds$V3[ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        ds.homD <- quantile(ds$V2[ds.ind],probs=c(.05,0.5,.95),na.rm=T,names=F)
        ## calculate extreme overdominance
                                        #heterozygously beneficial
        ua.hetBen <- which(ua$V3>ua.hetD[2])
                                        #homozygously deleterious
        ua.homDel <- which(ua$V2<ua.homD[1])

        da.hetBen <- which(da$V3>da.hetD[2])
        da.homDel <- which(da$V2<da.homD[1])
        us.hetBen <- which(us$V3>us.hetD[2])
        us.homDel <- which(us$V2<us.homD[1])
        ds.hetBen <- which(ds$V3>ds.hetD[2])
        ds.homDel <- which(ds$V2<ds.homD[1])

        ua.eo <- intersect(ua.hetBen,ua.homDel)
        da.eo <- intersect(da.hetBen,da.homDel)
        us.eo <- intersect(us.hetBen,us.homDel)
        ds.eo <- intersect(ds.hetBen,ds.homDel)
        allORFs <- unique(c(ua$V1,da$V1,us$V1,ds$V1))
        uaBool <- allORFs %in% ua$V1
        daBool <- allORFs %in% da$V1
        usBool <- allORFs %in% us$V1
        dsBool <- allORFs %in% ds$V1
        barcodeRepr <- uaBool+daBool+usBool+dsBool #how many barcodes we have for each deletion
                                        #for each ORF in allORFs:
        barcodehetBen <- allORFs %in% ua[ua.hetBen,1] + allORFs %in% us[us.hetBen,1] + allORFs %in% da[da.hetBen,1] + allORFs %in% ds[ds.hetBen,1]
        
        barcodehomDel <- allORFs %in% ua[ua.homDel,1] + allORFs %in% us[us.homDel,1] + allORFs %in% da[da.homDel,1] + allORFs %in% ds[ds.homDel,1]
        
        barcodeeo <- allORFs %in% ua[ua.eo,1] + allORFs %in% da[da.eo,1] + allORFs %in% us[us.eo,1] + allORFs %in% ds[ds.eo,1]
        
        hetBen <- allORFs[barcodehetBen>minBarNo & barcodehetBen/barcodeRepr == 1]
        homDel <- allORFs[barcodehomDel>minBarNo & barcodehomDel/barcodeRepr == 1]
        eo <- allORFs[barcodeeo>minBarNo & barcodeeo/barcodeRepr == 1]
        
        YPDdel.hetBen <- which(hetBen %in% YPDdel[,1])
        YPDben.hetBen <- which(hetBen %in% YPDben[,1])
        YPDdel.eo <- which(eo %in% YPDdel[,1])
        YPDben.eo <- which(eo %in% YPDben[,1])
        
                                        #filter out MDR ORFs
        hetBen <- setdiff(hetBen,notTrusted)
        eo <- setdiff(eo,notTrusted)
        
#filter out YPD beneficial or deleterious mutations
        hetBen <- setdiff(hetBen,YPDdel.hetBen)
        hetBen <- setdiff(hetBen,YPDben.hetBen)
        eo <- setdiff(eo,YPDdel.eo)
        eo <- setdiff(eo,YPDben.eo)

        #how many are in the rest 60 and how many are in a random sample of 60
        eoLN <- append(eoLN,sum(rest%in%eo))
        eoLR <- append(eoLR,sum(rrest%in%eo))
        hetBenLN <- append(hetBenLN,sum(rest%in%hetBen))
        hetBenLR <- append(hetBenLR,sum(rrest%in%hetBen))
                                        #        totalLN <- append(totalLN,length(allORFs))
    }
    Weo <- append(Weo,wilcox.test(eoLN,eoLR,alternative="less")$p.value)
    Whb <- append(Whb,wilcox.test(hetBenLN,hetBenLR,alternative="less")$p.value)
    print(paste(i))
}
sum(p.adjust(Weo)<0.05)/length(Weo)*100
sum(p.adjust(Whb)<0.05)/length(Whb)*100
#with old notTrusted 123
#[1] 80.18018
#[1] 29.72973



########################
# Hoepfner et al.,2014 #
########################

#./matchHoepfner.pl > ../Hoepfner/groups.dat

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

#first pass find MDR
#--------------------

# experiment with parallel loops for faster processing
## library(doParallel)
## registerDoParallel(cores=8)
## foreach(i=c(1:300)) %dopar% sqrt(i)
############

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

#make plot of multidrug resistance
h <- hist(filteredORFs,breaks=c(0:600),plot=F)
save(h,file="multidrugHop.h.dat")
save(filteredORFs,file="multidrugHop.filteredORFs.dat")

#estimate Poisson distribution with the same mean
bootPoisson <- matrix(nrow=1000,ncol=600)
for (r in c(1:1000)){
    hb <- hist(rpois(length(filteredORFs[filteredORFs!=0]),mean(filteredORFs[filteredORFs!=0])),breaks=c(0:600),plot=F)
    bootPoisson[r,] <- hb$counts
}


pdf("~/projects/dissertation/thesis/img/chapter3/multidrugHoepfner.pdf")
par(las=1,bty="l")
heights <- numeric(0)
ul <- numeric(0)
ll <- numeric(0)
for(i in c(1:600)){
  qboot <- quantile(bootPoisson[,i],probs=c(0.25,0.5,0.75),names=F)
  heights <- append(heights,qboot[2])
  ul <- append(ul,qboot[3])
  ll <- append(ll,qboot[1])
}
logHeights <- log10(heights)
logHeights[which(!is.finite(logHeights))]<-0
plot(h$mids+0.5,log10(h$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,450)) #set axes etc
rect(h$mids,0,h$mids+1,logHeights,col="grey",border="grey")
arrows(h$mids+.5,log10(ll),h$mids+.5,log10(ul),code=0,angle=90)
points(h$mids+0.5,log10(h$counts),pch=19)
dev.off()
  #remove ORFs with multidrug resistance as possibly carrying secondary mutations
#use therio
#notTrusted<-which(filteredORFs>200)
#write.table(hop[notTrusted,1],file="~/projects/yeastOverdominance/collection/Hoepfner/HopNotTrusted.dat",quote=F,row.names=F,col.names=F)
        
#second pass: dubious and random sets
#-----------------------------------

#calculate ranges nss, dss and rss to be used in the thrird pass
#where we can compare the number of HetBen, HetDel, HomBen and HomDel for each group We expect that under the random and dubious sets we will have less HetDel and HomDel (larger variance) across all conditions

totalL<-length(unq)

rss <- matrix(nrow=totalL,ncol=8)
dss <- matrix(nrow=totalL,ncol=8)
nss <- matrix(nrow=totalL,ncol=8)
counter <- 1

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
    bootReplicates <- 1000
    ds.homD  <- matrix(nrow=bootReplicates,ncol=2)
    ds.homDZ <- matrix(nrow=bootReplicates,ncol=2)
    ds.hetD  <- matrix(nrow=bootReplicates,ncol=2)
    ds.hetDZ <- matrix(nrow=bootReplicates,ncol=2)
    rs.homD  <- matrix(nrow=bootReplicates,ncol=2)
    rs.homDZ <- matrix(nrow=bootReplicates,ncol=2)
    rs.hetD  <- matrix(nrow=bootReplicates,ncol=2)
    rs.hetDZ <- matrix(nrow=bootReplicates,ncol=2)
    
    for(boot in c(1:bootReplicates)){
        dds <- sample(allDubious$V2,120) #dubious dubious ORFs
        dboot <- which(hop[,1] %in% dds)
        d.homD<-quantile(hopColumn[dboot],probs=c(.05,.95),na.rm=T,names=F)
        d.homDZ<-quantile(hopColumnZ[dboot],probs=c(.05,.95),na.rm=T,names=F)
        d.hetD<-quantile(hipColumn[dboot],probs=c(.05,.95),na.rm=T,names=F)
        d.hetDZ<-quantile(hipColumnZ[dboot],probs=c(.05,.95),na.rm=T,names=F)

        rds <- sample(length(hip[,1]),120) #random ORFs
        r.homD  <- quantile(hopColumn[rds],probs=c(.05,.95),na.rm=T,names=F)
        r.homDZ <- quantile(hopColumnZ[rds],probs=c(.05,.95),na.rm=T,names=F)
        r.hetD  <- quantile(hipColumn[rds],probs=c(.05,.95),na.rm=T,names=F)
        r.hetDZ <- quantile(hipColumnZ[rds],probs=c(.05,.95),na.rm=T,names=F)
         
        #save
        ds.homD[boot,] <-  d.homD
        ds.homDZ[boot,] <- d.homDZ
        ds.hetD[boot,]  <- d.hetD
        ds.hetDZ[boot,] <- d.hetDZ

        rs.homD[boot,]  <- r.homD
        rs.homDZ[boot,] <- r.homDZ
        rs.hetD[boot,]  <- r.hetD
        rs.hetDZ[boot,] <- r.hetDZ
    }

    nss[counter,] <- c(hipD,hipDZ,hopD,hopDZ)
    dss[counter,] <- c(colMeans(ds.hetD),colMeans(ds.hetDZ),colMeans(ds.homD),colMeans(ds.homDZ))
    rss[counter,] <- c(colMeans(rs.hetD),colMeans(rs.hetDZ),colMeans(rs.homD),colMeans(rs.homDZ))
    counter <- counter + 1
}

#third pass
#----------
nsum <- matrix(nrow=totalL,ncol=5)
dsum <- matrix(nrow=totalL,ncol=5)
rsum <- matrix(nrow=totalL,ncol=5)
## pdf("temp.pdf")
## par(mfrow=c(2,2))
counter <- 1
condNames <- character()
condID <- character()
concs <- numeric(0)
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

   print(paste(titleString,counter,"/",length(unq)))
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
   HetBeneficial <- setdiff(HetBeneficial,notTrusted)
   
   HetDel <- which(hipColumn<nss[counter,1])
   HetDelZ <- which(hipColumnZ<nss[counter,3])
   HetDeleterious <- intersect(HetDel,HetDelZ)
   HetDeleterious <- setdiff(HetDeleterious,notTrusted)
   
   HomBen <- which(hopColumn>nss[counter,6])
   HomBenZ <- which(hopColumnZ>nss[counter,8])
   HomBeneficial <- intersect(HomBen,HomBenZ)
   HomBeneficial <- setdiff(HomBeneficial,notTrusted)
   
   HomDel <- which(hopColumn<nss[counter,5])
   HomDelZ <- which(hopColumnZ<nss[counter,7])
   HomDeleterious <- intersect(HomDel,HomDelZ)
   HomDeleterious <- setdiff(HomDeleterious,notTrusted)

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
   dHetBeneficial <- setdiff(dHetBeneficial,notTrusted)   
   dHetDel <- which(hipColumn<dss[counter,1])
   dHetDelZ <- which(hipColumnZ<dss[counter,3])
   dHetDeleterious <- intersect(dHetDel,dHetDelZ)
   dHetDeleterious <- setdiff(dHetDeleterious,notTrusted)
   dHomBen <- which(hopColumn>dss[counter,6])
   dHomBenZ <- which(hopColumnZ>dss[counter,8])
   dHomBeneficial <- intersect(dHomBen,dHomBenZ)
   dHomBeneficial <- setdiff(dHomBeneficial,notTrusted)
   dHomDel <- which(hopColumn<dss[counter,5])
   dHomDelZ <- which(hopColumnZ<dss[counter,7])
   dHomDeleterious <- intersect(dHomDel,dHomDelZ)
   dHomDeleterious <- setdiff(dHomDeleterious,notTrusted)
   deo <- intersect(dHetBeneficial,dHomDeleterious)
   dsum[counter,] <- c(length(dHetBeneficial),length(dHetDeleterious),length(dHomBeneficial),length(dHomDeleterious),length(deo))

   rHetBen <- which(hipColumn>rss[counter,2])
   rHetBenZ <- which(hipColumnZ>rss[counter,4])
   rHetBeneficial <- intersect(rHetBen,rHetBenZ)
   rHetBeneficial <- setdiff(rHetBeneficial,notTrusted)   
   rHetDel <- which(hipColumn<rss[counter,1])
   rHetDelZ <- which(hipColumnZ<rss[counter,3])
   rHetDeleterious <- intersect(rHetDel,rHetDelZ)
   rHetDeleterious <- setdiff(rHetDeleterious,notTrusted)
   rHomBen <- which(hopColumn>rss[counter,6])
   rHomBenZ <- which(hopColumnZ>rss[counter,8])
   rHomBeneficial <- intersect(rHomBen,rHomBenZ)
   rHomBeneficial <- setdiff(rHomBeneficial,notTrusted)
   rHomDel <- which(hopColumn<rss[counter,5])
   rHomDelZ <- which(hopColumnZ<rss[counter,7])
   rHomDeleterious <- intersect(rHomDel,rHomDelZ)
   rHomDeleterious <- setdiff(rHomDeleterious,notTrusted)
   reo <- intersect(rHetBeneficial,rHomDeleterious)
   rsum[counter,] <- c(length(rHetBeneficial),length(rHetDeleterious),length(rHomBeneficial),length(rHomDeleterious),length(reo))

   
   if(0){
       pdf("temp.pdf")
       par(mfcol=c(2,2))
                                        #plot
       plot(hipColumnZ,hipColumn,cex=0.5,col="#00000070",main=paste("HIP",titleString),pch=19)
       points(hipColumnZ[HetBeneficial],hipColumn[HetBeneficial],col="#ff000060",pch=19)
       points(hipColumnZ[HomDeleterious],hipColumn[HomDeleterious],col="#0000ff60",pch=19)

       points(hipColumnZ[5573],hipColumn[5573],col="#00ff0060",pch=19)
       
       abline(h=nss[counter,c(1,2)],col="red")
       abline(v=nss[counter,c(3,4)],col="red")
       abline(h=dss[counter,c(1,2)],col="blue")
       abline(v=dss[counter,c(3,4)],col="blue")
       abline(h=rss[counter,c(1,2)],col="black")
       abline(v=rss[counter,c(3,4)],col="black")
       
       plot(hopColumnZ,hopColumn,cex=0.5,col="#00000070",main=paste("HOP",titleString),pch=19)
       points(hopColumnZ[HetBeneficial],hopColumn[HetBeneficial],col="#ff000060",pch=19)
       points(hopColumnZ[HomDeleterious],hopColumn[HomDeleterious],col="#0000ff60",pch=19)

       points(hopColumnZ[5573],hopColumn[5573],col="#00ff0060",pch=19)
       
       abline(h=nss[counter,c(5,6)],col="red")
       abline(v=nss[counter,c(7,8)],col="red")
       abline(h=dss[counter,c(5,6)],col="blue")
       abline(v=dss[counter,c(7,8)],col="blue")
       abline(h=rss[counter,c(5,6)],col="black")
       abline(v=rss[counter,c(7,8)],col="black")
   }
   counter <- counter + 1
}
# dev.off()


#compare ranges by comparing proportion of hetben hetdel homben homdel
colSums((nsum[,c(1:4)]-dsum[,c(1:4)])<0)/totalL*100 #neutral - dubious
colSums((nsum[,c(1:4)]-rsum[,c(1:4)])<0)/totalL*100 #neutral - random
colSums((dsum[,c(1:4)]-rsum[,c(1:4)])<0)/totalL*100 #dubious - random

#[1] 44.47876 46.94981 29.80695 20.69498
#[1] 43.47490 45.55985 26.67954 18.49421
#[1] 51.50579 48.33977 36.79537 33.89961

 ## 56.52510 53.28185 71.11969 79.49807
 ## 56.71815 54.74903 73.62934 81.73745
 ## 46.91120 51.00386 61.93050 64.63320



#######>>>>>>>>>>mexri edo meta to 666 <<<<<<<,,
#mean percentage of het. beneficial across conditions
mean(nsum[,1])/length(hip[,1])*100
#[1] 4.397615

#mean percentage of extreme overdominance
mean(nsum[,5]/nsum[,1])*100
#[1] 3.911677

hetben   eo
eoH<-data.frame(condNames,concs,nsum[,5]/nsum[,1],nsum[,1],nsum[,5])
eoH[order(eoH[,3]),]
write.csv(eoH,file="~/projects/yeastOverdominance/collection/Hoepfner/odRatio.csv",row.names=F)
l <- length(eoH[,1])
pdf("eoRatioHoepfner.pdf")
par(las=1,bty="l")
orderedHo <- eoH[order(eoH[,3]),]
ci <- binom.confint(orderedHo[,5],orderedHo[,4],method="exact")
plot(orderedHo[,3]*100,ylab="extreme overdominance (%)",xlab="conditions",type="l",ylim=c(0,30))
polygon(c(c(1:l),c(l:1)),c(ci$lower,rev(ci$upper))*100,col="grey60",border=NA)
points(orderedHo[,3]*100,type="l",lwd=2)
dev.off()

#fourth pass
#-----------


# Overlap of hetBen across experiments
######################################
#Overlap hetben Hillenmeyer et al
#--------------------------------
#find conditions with multiple timepoints
#calculate overlap
path <- "~/projects/yeastOverdominance/collection/analysis/odORFs/"
files <- character()
files[1] <- paste(path,"23DEGREESC15oldscannerhetBen.dat",sep="")
files[2] <- paste(path,"23DEGREESC5oldscannerhetBen.dat",sep="")
files[3] <- paste(path,"37DEGREESC20oldscannerhetBen.dat",sep="")
files[4] <- paste(path,"37DEGREESC5oldscannerhetBen.dat",sep="")
files[5] <- paste(path,"5FLUOROURACIL15oldscannerhetBen.dat",sep="")
files[6] <- paste(path,"5FLUOROURACIL20oldscannerhetBen.dat",sep="")
files[7] <- paste(path,"5FLUOROURACIL5oldscannerhetBen.dat",sep="")
files[8] <- paste(path,"FLOXURIDINE20oldscannerhetBen.dat",sep="")
files[9] <- paste(path,"FLOXURIDINE5oldscannerhetBen.dat",sep="")
files[10] <- paste(path,"YPGLYCEROL15oldscannerhetBen.dat",sep="")
files[11] <- paste(path,"YPGLYCEROL5oldscannerhetBen.dat",sep="")
generations <- c(15,5,20,5,15,20,5,20,5,15,5)
env <- c("23DEGREESC","23DEGREESC","37DEGREESC","37DEGREESC","5FLUOROURACIL","5FLUOROURACIL","5FLUOROURACIL","FLOXURIDINE","FLOXURIDINE","YPGLYCEROL","YPGLYCEROL")
#pick random 2 with the same condition but different generation and caclulate overlap
bcomSame <- numeric(0)
bootNumber <- 10000
while(bootNumber > 0){
    pair <- sample(length(files),2) #pick a pair of files (without replacements)
    if(env[pair[1]]==env[pair[2]]){ #same condition
        fileA <- read.table(files[pair[1]])
        fileB <- read.table(files[pair[2]])
        bcomSame <- append(bcomSame,length(intersect(fileA[,1],fileB[,1]))/length(union(fileA[,1],fileB[,1])))
        bootNumber <- bootNumber - 1
        print(bootNumber)
    }
}
#pick random 2 with different condition
bcomDiff <- numeric(0)
bootNumber <- 10000
while(bootNumber > 0){
    pair <- sample(length(files),2) #pick a pair of files (without replacements)
    if(env[pair[1]]!=env[pair[2]]){ #diff condition
        fileA <- read.table(files[pair[1]])
        fileB <- read.table(files[pair[2]])
        bcomDiff <- append(bcomDiff,length(intersect(fileA[,1],fileB[,1]))/length(union(fileA[,1],fileB[,1])))
        bootNumber <- bootNumber - 1
        print(bootNumber)
    }
}
hist(bcomSame)
mean(bcomSame)
hist(bcomDiff)
mean(bcomDiff)
wilcox.test(bcomSame,bcomDiff,alternative="greater")


#Overlap hetben Hoepfner et al
#------------------------------
df <- data.frame(condID,condNames,as.numeric(concs),nsum[,1:5]) #hetben, hetdel, homben, homdel
indexD <- duplicated(df[,1]) #compounds with only one concentration
dpl <- unique(df[indexD,1])


#compare overlap of heterozygous beneficial across compounds and across conditions
path <- "~/projects/yeastOverdominance/collection/Hoepfner/odORFs/"
files <- dir(path=path, pattern="*hetBen.dat")
#find conditions for which we have more than one concentration
conditions <- character(0)
for(f in files){
    conditions <- append(conditions,strsplit(f,".",fixed=T)[[1]][1])
}

filesD <- which(duplicated(conditions))
#perform 1000 bootstrap measurements, to measure overlap across measurements with the same compound
bcomSame <- numeric(0)
bootNumber <- 10000
while(bootNumber > 0){
    pair <- sample(filesD,2) #pick a pair of files
    fileA <- read.table(paste(path,files[pair[1]],sep=""))
    fileB <- read.table(paste(path,files[pair[2]],sep=""))
    compoundA <- strsplit(files[pair[1]],".",fixed=T)[[1]][1]
    compoundB <- strsplit(files[pair[2]],".",fixed=T)[[1]][1]
    if(compoundA == compoundB){
        bcomSame <- append(bcomSame,length(intersect(fileA[,1],fileB[,1]))/length(union(fileA[,1],fileB[,1])))
        bootNumber <- bootNumber - 1
        print(bootNumber)
    }
}

#repeat bootstrap for comp
bootNumber <- 10000
bcomDiff <- numeric(0)
while(bootNumber > 0){
    pair <- sample(filesD,2) #pick a pair
    fileA <- read.table(paste(path,files[pair[1]],sep=""))
    fileB <- read.table(paste(path,files[pair[2]],sep=""))
    compoundA <- strsplit(files[pair[1]],".",fixed=T)[[1]][1]
    compoundB <- strsplit(files[pair[2]],".",fixed=T)[[1]][1]
    if(compoundA != compoundB){
        bcomDiff <- append(bcomDiff,length(intersect(fileA[,1],fileB[,1]))/length(union(fileA[,1],fileB[,1])))
        bootNumber <- bootNumber - 1
        print(bootNumber)
    }
}

wilcox.test(bcomSame,bcomDiff,alternative="greater")


#Overlap between datasets
#-----------------------

#common conditions
pathHop <- "~/projects/yeastOverdominance/collection/Hoepfner/odORFs/"
pathHil <- "~/projects/yeastOverdominance/collection/analysis/odORFs/"
fileHil <- numeric()
fileHop <- numeric()

fileHil[1] <- "ACLACINOMYCINA20newscannerhetBen.dat"
fileHil[2] <- "COCL220newscannerhetBen.dat"
fileHil[3] <- "CUSO420newscannerhetBen.dat"
fileHil[4] <- "CYCLOHEXIMIDE20newscannerhetBen.dat"
fileHil[5] <- "LATRUNCULIN20newscannerhetBen.dat"
fileHil[6] <- "MECHLORETHAMINE20oldscannerhetBen.dat"
fileHil[7] <- "HGCL220newscannerhetBen.dat"
fileHil[8] <- "METHOTREXATE20oldscannerhetBen.dat"
fileHil[9] <- "MYCOPHENOLICACID20oldscannerhetBen.dat"
fileHil[10] <- "RAPAMYCIN20newscannerhetBen.dat"
fileHil[11] <- "ZNCL220newscannerhetBen.dat"
compoundHil <- c("AclacinomycinA","CoCl2","CuSO4","Cycloheximide","Latrunculin","Mechlorethamine","HgCl2","Methotrexate","MycophenolicAcid","Rapamycin","ZnCl2")

fileHop[1] <- "Aclacinomycin A, Aclarubicin.4.857hetBen.dat"
fileHop[2] <- "Cobalt(Il)-Chlorid.200hetBen.dat"
fileHop[3] <- "Cobalt(Il)-Chlorid.92.63hetBen.dat"
fileHop[4] <- "Copper(II)Sulfate.120hetBen.dat"
fileHop[5] <- "Copper(II)Sulfate.200hetBen.dat"
fileHop[6] <- "Copper(II)Sulfate.50hetBen.dat"
fileHop[7] <- "Copper(II)Sulfate.75hetBen.dat"
fileHop[8] <- "Cycloheximide.0.03hetBen.dat"
fileHop[9] <- "Cycloheximide.0.05hetBen.dat"
fileHop[10] <- "Latrunculin A.0.7hetBen.dat"
fileHop[11] <- "Latrunculin A.0.9hetBen.dat"
fileHop[12] <- "Mechlorethamine.95hetBen.dat"
fileHop[13] <- "Mercury(II) chloride.35hetBen.dat"
fileHop[14] <- "Mercury(II) chloride.50hetBen.dat"
fileHop[15] <- "Methotrexate.200hetBen.dat"
fileHop[16] <- "Mycophenolic Acid (Myfortic).100hetBen.dat"
fileHop[17] <- "Mycophenolic Acid (Myfortic).120hetBen.dat"
fileHop[18] <- "Mycophenolic Acid (Myfortic).70hetBen.dat"
fileHop[19] <- "Rapamycin.0.001hetBen.dat"
fileHop[20] <- "Rapamycin.2e-04hetBen.dat"
fileHop[21] <- "Rapamycin.5e-04hetBen.dat"
fileHop[22] <- "Zinc Chloride.20000hetBen.dat"
fileHop[23] <- "Zinc Chloride.200hetBen.dat"
fileHop[24] <- "Zinc Chloride.40000hetBen.dat"
compoundHop <- c("AclacinomycinA","CoCl2","CoCl2","CuSO4","CuSO4","CuSO4","CuSO4","Cycloheximide","Cycloheximide","Latrunculin","Latrunculin","Mechlorethamine","HgCl2","HgCl2","Methotrexate","MycophenolicAcid","MycophenolicAcid","MycophenolicAcid","Rapamycin","Rapamycin","Rapamycin","ZnCl2","ZnCl2","ZnCl2")

#pick one file from each dataset and calculate overlap
overlap <- function(fileHil,fileHop){
    ORFHil1 <- read.table(paste(pathHil,fileHil,sep=""))
    ORFHop1 <- read.table(paste(pathHop,fileHop,sep=""))
    o <- length(intersect(ORFHil1[,1],ORFHop1[,1]))
    n <- length(union(ORFHil1[,1],ORFHop1[,1]))
    c(o/n,n)
}

overlap(fileHil[1],fileHop[1])

overlap(fileHil[2],fileHop[2])
overlap(fileHil[2],fileHop[3])

overlap(fileHil[3],fileHop[4])
overlap(fileHil[3],fileHop[5])
overlap(fileHil[3],fileHop[6])
overlap(fileHil[3],fileHop[7])

overlap(fileHil[4],fileHop[8])
overlap(fileHil[4],fileHop[9])

overlap(fileHil[5],fileHop[10])
overlap(fileHil[5],fileHop[11])

overlap(fileHil[6],fileHop[12])

overlap(fileHil[7],fileHop[13])
overlap(fileHil[7],fileHop[14])

overlap(fileHil[8],fileHop[15])

overlap(fileHil[9],fileHop[16])
overlap(fileHil[9],fileHop[17])
overlap(fileHil[9],fileHop[18])

overlap(fileHil[10],fileHop[19])
overlap(fileHil[10],fileHop[20])
overlap(fileHil[10],fileHop[21])

overlap(fileHil[11],fileHop[22])
overlap(fileHil[11],fileHop[23])
overlap(fileHil[11],fileHop[24])


#same compounds
bSame <- numeric(0)
bootNumber <- 10000
while(bootNumber > 0){
    indexHil <- sample(length(compoundHil),1)
    indexHop <- sample(length(compoundHop),1)
    if(compoundHil[indexHil] == compoundHop[indexHop]){
        ORFHil <- read.table(paste(pathHil,fileHil[indexHil],sep=""))
        ORFHop <- read.table(paste(pathHop,fileHop[indexHop],sep=""))
        bSame <- append(bSame,length(intersect(ORFHil[,1],ORFHop[,1]))/length(union(ORFHil[,1],ORFHop[,1])))
        bootNumber <- bootNumber - 1
        print(bootNumber)
    }   
}

#different compounds
bDiff <- numeric(0)
bootNumber <- 10000
while(bootNumber > 0){
    indexHil <- sample(length(compoundHil),1)
    indexHop <- sample(length(compoundHop),1)
    if(compoundHil[indexHil] != compoundHop[indexHop]){
        ORFHil <- read.table(paste(pathHil,fileHil[indexHil],sep=""))
        ORFHop <- read.table(paste(pathHop,fileHop[indexHop],sep=""))
        bDiff <- append(bDiff,length(intersect(ORFHil[,1],ORFHop[,1]))/length(union(ORFHil[,1],ORFHop[,1])))
        bootNumber <- bootNumber - 1
        print(bootNumber)
    }   
}

mean(bSame)
mean(bDiff)
wilcox.test(bSame,bDiff,alternative="greater")


#---------------
#read genes from Daniel and Dglucose
path <- "~/projects/yeastOverdominance/collection/Hoepfner/odORFs/"
g1o <- read.table(paste(path,"D-Glucose (starvation).0.25eo.dat",sep=""))
g1b <- read.table(paste(path,"D-Glucose (starvation).0.25hetBen.dat",sep=""))
g2o <- read.table(paste(path,"D-Glucose (starvation).0.5eo.dat",sep=""))
g2b <- read.table(paste(path,"D-Glucose (starvation).0.5hetBen.dat",sep=""))
g3o <- read.table(paste(path,"D-Glucose (starvation).0.75eo.dat",sep=""))
g3b <- read.table(paste(path,"D-Glucose (starvation).0.75hetBen.dat",sep=""))
ks <- read.table("~/projects/yeastOverdominance/chemostat/sequencing/colony/mutations/hapdip/systName.csv",sep="\t")
#list of mutated genes in experimental evolution
systN <- c("YDL194W","YDL190C","YDR211W","YDR443C","YEL035C","YER157W","YFR020W","YGL023","YHR120W","YJL190C","YMR240C","HXT6","HXT7")

shortN <- c("SNF3","UFD2","GCD6","SSN2","UTR5","YER157W","YFR020W","YGL023","YHR120W","YJL190C","YMR240C","HXT6","HXT7")


intersect(g1b$V1,ks$V1)
[1] "YDL194W" "YER133W" "YJR115W" "YPR049C"
intersect(g2b$V1,ks$V1)
1] "YAL056W" "YAR019C" "YGL250W" "YIR023W" "YNL098C" "YNR031C" "YOR307C"
[8] "YOR360C" "YPL016W"
intersect(g3b$V1,ks$V1)
[1] "YAR035W" "YNL098C"
also MIG2/YGL209W present

> intersect(g1o$V1,g2o$V1)
[1] "YER017C" "YIL036W" "YMR063W"
> intersect(g1o$V1,g3o$V1)
[1] "YER017C"
> intersect(g2o$V1,g3o$V1)
[1] "YER017C"
> 
intersect(systN,g1b$V1)
YDL194W #SNF3
intersect(systN,g2b$V1)
intersect(systN,g3b$V1)

#           AFG3 common in all
#  YER017C

#make plots for D-Glucose 
which(condNames=="D-Glucose (starvation)")
956 1073 2099

counter <- 2099
i<-unq[counter]
STD1<-5800
SNF3<- 1018
i

pdf("temp.pdf")
par(mfcol=c(2,2))
                                        #plot
plot(hipColumnZ,hipColumn,cex=0.5,col="#00000070",main=paste("HIP",titleString),pch=19)
points(hipColumnZ[HetBeneficial],hipColumn[HetBeneficial],col="#ff000060",pch=19)
points(hipColumnZ[HomDeleterious],hipColumn[HomDeleterious],col="#0000ff60",pch=19)
points(hipColumnZ[SNF3],hipColumn[SNF3],col="#00ff0099",pch=19)

abline(h=nss[counter,c(1,2)],col="red")
abline(v=nss[counter,c(3,4)],col="red")
abline(a=0,b=1)       

plot(hopColumnZ,hopColumn,cex=0.5,col="#00000070",main=paste("HOP",titleString),pch=19)
points(hopColumnZ[HetBeneficial],hopColumn[HetBeneficial],col="#ff000060",pch=19)
points(hopColumnZ[HomDeleterious],hopColumn[HomDeleterious],col="#0000ff60",pch=19)

points(hopColumnZ[SNF3],hopColumn[SNF3],col="#00ff0099",pch=19)

abline(h=nss[counter,c(5,6)],col="red")
abline(v=nss[counter,c(7,8)],col="red")
abline(a=0,b=1)       
dev.off()



#make  plots for CuSO4
#CuSO4 in Hoepfner
which(grepl("Copper",condNames))
709  840 1819 2154

counter <- 709
i<-unq[counter]


source("~/projects/fgmo/colors.R")
pdf("HopCuSO4.pdf",width=10,height=6)
par(mfrow=c(1,2),bty="l",lwd=2,las=1)
alpha <- "ff"
                                        #plot
plot(hipColumnZ,hipColumn,cex=1,col="#00000040",main="Heterozygous",pch=19,xlim=c(-4,4),ylim=c(-10,10),xlab="z-score", ylab="sensitivity")
points(hipColumnZ[HetBeneficial],hipColumn[HetBeneficial],col=paste(cgreen4,alpha,sep=""),pch=19)
abline(h=nss[counter,c(1,2)],col="brown",lty=2)
abline(v=nss[counter,c(3,4)],col="brown",lty=2)
abline(a=0,b=1)       
points(hipColumnZ[eo],hipColumn[eo],col=paste(dblue,"aa",sep=""),pch=19)

plot(hopColumnZ,hopColumn,cex=1,col="#00000040",main="Homozygous",pch=19,xlab="z-score", ylab="sensitivity",xlim=c(-4,4),ylim=c(-10,10))
points(hopColumnZ[HetBeneficial],hopColumn[HetBeneficial],col=paste(cgreen4,alpha,sep=""),pch=19)
abline(h=nss[counter,c(5,6)],col="brown",lty=2)
abline(v=nss[counter,c(7,8)],col="brown",lty=2)
abline(a=0,b=1)
points(hopColumnZ[eo],hopColumn[eo],col=paste(dblue,"aa",sep=""),pch=19)
dev.off()


# print summary of common conditions across datasets
#common conditions
pathHop <- "~/projects/yeastOverdominance/collection/Hoepfner/odORFs/"
pathHil <- "~/projects/yeastOverdominance/collection/analysis/odORFs/"
fileHil <- numeric()
fileHop <- numeric()

fileHil[1] <- "ACLACINOMYCINA20newscanner"
fileHil[2] <- "COCL220newscanner"
fileHil[3] <- "CUSO420newscanner"
fileHil[4] <- "CYCLOHEXIMIDE20newscanner"
fileHil[5] <- "LATRUNCULIN20newscanner"
fileHil[6] <- "MECHLORETHAMINE20oldscanner"
fileHil[7] <- "HGCL220newscanner"
fileHil[8] <- "METHOTREXATE20oldscanner"
fileHil[9] <- "MYCOPHENOLICACID20oldscanner"
fileHil[10] <- "RAPAMYCIN20newscanner"
fileHil[11] <- "ZNCL220newscanner"
compoundHil <- c("AclacinomycinA","CoCl2","CuSO4","Cycloheximide","Latrunculin","Mechlorethamine","HgCl2","Methotrexate","MycophenolicAcid","Rapamycin","ZnCl2")

fileHop[1] <- "Aclacinomycin A, Aclarubicin.4.857"
fileHop[2] <- "Cobalt(Il)-Chlorid.200"
fileHop[3] <- "Cobalt(Il)-Chlorid.92.63"
fileHop[4] <- "Copper(II)Sulfate.120"
fileHop[5] <- "Copper(II)Sulfate.200"
fileHop[6] <- "Copper(II)Sulfate.50"
fileHop[7] <- "Copper(II)Sulfate.75"
fileHop[8] <- "Cycloheximide.0.03"
fileHop[9] <- "Cycloheximide.0.05"
fileHop[10] <- "Latrunculin A.0.7"
fileHop[11] <- "Latrunculin A.0.9"
fileHop[12] <- "Mechlorethamine.95"
fileHop[13] <- "Mercury(II) chloride.35"
fileHop[14] <- "Mercury(II) chloride.50"
fileHop[15] <- "Methotrexate.200"
fileHop[16] <- "Mycophenolic Acid (Myfortic).100"
fileHop[17] <- "Mycophenolic Acid (Myfortic).120"
fileHop[18] <- "Mycophenolic Acid (Myfortic).70"
fileHop[19] <- "Rapamycin.0.001"
fileHop[20] <- "Rapamycin.2e-04"
fileHop[21] <- "Rapamycin.5e-04"
fileHop[22] <- "Zinc Chloride.20000"
fileHop[23] <- "Zinc Chloride.200"
fileHop[24] <- "Zinc Chloride.40000"
compoundHop <- c("AclacinomycinA","CoCl2","CoCl2","CuSO4","CuSO4","CuSO4","CuSO4","Cycloheximide","Cycloheximide","Latrunculin","Latrunculin","Mechlorethamine","HgCl2","HgCl2","Methotrexate","MycophenolicAcid","MycophenolicAcid","MycophenolicAcid","Rapamycin","Rapamycin","Rapamycin","ZnCl2","ZnCl2","ZnCl2")

library(binom)
for (i in c(1:length(compoundHil))){
    arpa <- try(ORFHilhb <- read.table(paste(pathHil,fileHil[i],"hetBen.dat",sep="")))
    if(inherits(arpa,"try-error")){
        HilHb <- 0
    }else{
        HilHb <- length(ORFHilhb[,1])
    }
    arpa <- try(ORFHileo <- read.table(paste(pathHil,fileHil[i],"eo.dat",sep="")))
    if(inherits(arpa,"try-error")){
        Hileo <- 0
    }else{
        Hileo <- length(ORFHileo[,1])
    }
    HilBinom <- binom.confint(Hileo,HilHb,method="exact")
    for(j in which(compoundHop==compoundHil[i])){
        arpa <- try(ORFHophb <- read.table(paste(pathHop,fileHop[j],"hetBen.dat",sep="")))
        if(inherits(arpa,"try-error")){
            HopHb <- 0
        }else{
            HopHb <- length(ORFHophb[,1])
        }
        arpa <- try(ORFHopeo <- read.table(paste(pathHop,fileHop[j],"eo.dat",sep="")))
        if(inherits(arpa,"try-error")){
            Hopeo <- 0
        }else{
            Hopeo <- length(ORFHileo[,1])
        }
        HopBinom <- binom.confint(Hopeo,HopHb,method="exact")
#        print(paste(fileHil[i],fileHop[j],Hileo/HilHb*100,HilBinom$lower*100,HilBinom$upper*100,Hopeo/HopHb*100,HopBinom$lower*100,HopBinom$upper*100),sep="\t")
        print(paste(fileHil[i],fileHop[j],HilHb,HopHb))
    }
}

#subsample ORFs
pick random 60 from neutral ORFs set
use them to calculate
1. extreme overdominance in the rest of the 60 
2. and in 60 other random ORFs
we expect that 1. will be much lower than 2.
in a fourth pass
