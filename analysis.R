## Deletion collection analysis
## -----------------------------
##     analysis of deletion collection data. Separately for Hillenmeyer etal 2008 and Hoepfner et al., 2014 but with similar methodology

## - Match homozygote and heterozygote experiments and find multidrug resistance distributions
## - Exclude deletions that seem to be beneficial across multiple conditions
## - Calculate confidence intervals for dubious ORFs and random subset of ORFs
## - Find extreme overdominant


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
                                        #each row is a bootstrap. for each integer from min to max get mean and quartiles
png("~/projects/yeastOverdominance/manuscript/figures/multidrugHillenmeyer.png")
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
plot(h$mids+0.5,log10(h$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,60)) #set axes etc
rect(h$mids,0,h$mids+1,logHeights,col="grey")
arrows(h$mids+.5,log10(ll),h$mids+.5,log10(ul),code=0,angle=90)
points(h$mids+0.5,log10(h$counts),pch=19,type="o")
dev.off()
    
notTrusted <- names(mdr)[mdr>10]

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
sum(dsdfhom[,1]-dsdfhom[,2]<0)/111*100)/4

(sum(uadfhom[,1]-uadfhom[,3]<0)/111*100 +
sum(dadfhom[,1]-dadfhom[,3]<0)/111*100 +
sum(usdfhom[,1]-usdfhom[,3]<0)/111*100 +
sum(dsdfhom[,1]-dsdfhom[,3]<0)/111*100)/4

(sum(uadfhom[,2]-uadfhom[,3]<0)/111*100 +
sum(dadfhom[,2]-dadfhom[,3]<0)/111*100 +
sum(usdfhom[,2]-usdfhom[,3]<0)/111*100 +
sum(dsdfhom[,2]-dsdfhom[,3]<0)/111*100)/4


(sum(uadfhet[,1]-uadfhet[,2]<0)/111*100 +
sum(dadfhet[,1]-dadfhet[,2]<0)/111*100 +
sum(usdfhet[,1]-usdfhet[,2]<0)/111*100 +
sum(dsdfhet[,1]-dsdfhet[,2]<0)/111*100)/4

(sum(uadfhet[,1]-uadfhet[,3]<0)/111*100 +
sum(dadfhet[,1]-dadfhet[,3]<0)/111*100 +
sum(usdfhet[,1]-usdfhet[,3]<0)/111*100 +
sum(dsdfhet[,1]-dsdfhet[,3]<0)/111*100)/4

(sum(uadfhet[,2]-uadfhet[,3]<0)/111*100 +
sum(dadfhet[,2]-dadfhet[,3]<0)/111*100 +
sum(usdfhet[,2]-usdfhet[,3]<0)/111*100 +
sum(dsdfhet[,2]-dsdfhet[,3]<0)/111*100)/4

# third pass 
#-------------
minBarNo <- 1 # >minBarNo should be available for each gene 
eoL <- numeric(0) # number of extreme overdominant in each condition
hetBenL <- numeric(0)
eoD<- numeric(0) # number of extreme overdominant in each condition
hetBenD<- numeric(0)
eoR <- numeric(0) # number of extreme overdominant in each condition
hetBenR <- numeric(0)
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

    write.table(hetBen,file=paste("analysis/odORFs/",i,"hetBen.dat",sep=""),row.names=F,col.names=F)
    write.table(eo,file=paste("analysis/odORFs/",i,"eo.dat",sep=""),row.names=F,col.names=F)

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


png("eoRatioHillenmeyer.png",height=800,width=640)
par(las=1,bty="l",mai=c(.5,3,0,0.2),mgp=c(1.5,0.5,0),cex=0.8)
hb <- data.frame(conditions$V1,eoL/hetBenL,eoL,hetBenL)
orderedhB <- hb[order(hb[,2]),]
orderedhB <- orderedhB[which(orderedhB[,2]>0),]
ci <- binom.confint(orderedhB[,3],orderedhB[,4],method="exact")
bar <- barplot(orderedhB[,2],names.arg=orderedhB[,1],horiz=TRUE,xlim=c(0,0.55),xlab="extreme overdominance")
arrows(ci$lower,bar,ci$upper,bar,code=0)
dev.off()

########################
# Hoepfner et al.,2014 #
########################

#./matchHoepfner.pl > ../Hoepfner/groups.dat

a<-read.table("../Hoepfner/groups.dat")
#a$V1 compound
#a$V2 concentration
#a$V3 z-score?
#a$V4 HIP/HOP
#a$V5 column index
unq <- unique(paste(a$V1,a$V2)) #unique combinations of compound and concentration
for(i in unq){ #for each combination.
    row.index <- which(paste(a$V1,a$V2)==i)
    hipi  <- a[row.index,4] == "HIP" & a[row.index,3] == 0
    hipZi <- a[row.index,4] == "HIP" & a[row.index,3] == 1
    hopi  <- a[row.index,4] == "HOP" & a[row.index,3] == 0
    hopZi <- a[row.index,4] == "HOP" & a[row.index,3] == 1

    hipIndex <- a[row.index,][hipi,5]
    hopIndex <- a[row.index,][hopi,5]
    hipIndexZ <- a[row.index,][hipZi,5]
    hopIndexZ <- a[row.index,][hopZi,5]
}
library(binom)
library(stringr)

#read data files
hip<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HIP_scores.txt",sep="\t")
hop<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HOP_scores.txt",sep="\t")
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs
m<-read.table("~/projects/yeastOverdominance/collection/Hoepfner/match.dat")
d<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
dI<-which(hop[,1]%in%d[,1]) # 120 neutrals
n<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/compoundNames.csv") 

#first pass find MDR
#--------------------
filteredORFs<-rep(0,length(hip[,1])) # list of counters, one for each gene, to count how many times is heterozygously beneficail
for(i in c(1:(length(m$V1)))){ #for each condition
  hipIndex <- m[i,2]+1  # HIP sensitivity
  hopIndex <- m[i,3]+1  # HOP sensitivity
  hipIndexZ <- m[i,4]+1 # HIP z-score
  hopIndexZ <- m[i,5]+1 # HOP z-score

  #calculate percentiles
  hipD<-quantile(hip[dI,hipIndex],probs=c(.05,.95),na.rm=T,names=F)
  hopD<-quantile(hop[dI,hipIndex],probs=c(.05,.95),na.rm=T,names=F)
  hipDZ<-quantile(hip[dI,hipIndexZ],probs=c(.05,.95),na.rm=T,names=F)
  hopDZ<-quantile(hop[dI,hipIndexZ],probs=c(.05,.95),na.rm=T,names=F)
  
  #find outliers for HIP experiments
  HetBen<-which(hip[,hipIndex]>hipD[2]) # Row Index of strains with HIP larger than 95% 
  HetBenZ<-which(hip[,hipIndexZ]>hipDZ[2]) # Row Index of strains with HIP larger than 95% for z-score
  #intersection
  HetResistant <- intersect(HetBen,HetBenZ)
  #find multidrug resistant ORFs
  filteredORFs[HetResistant] <- filteredORFs[HetResistant] + 1
}

#make plot of multidrug resistance
h <- hist(filteredORFs,breaks=c(0:60),plot=F)
#estimate Poisson distribution with the same mean
bootPoisson <- data.frame()
for (r in c(1:1000)){
  hb <- hist(rpois(length(filteredORFs[filteredORFs!=0]),mean(filteredORFs[filteredORFs!=0])),breaks=c(0:60))
  bootPoisson <- rbind(bootPoisson,hb$counts)
}

#each row is a bootstrap. for each integer from min to max get mean and quartiles
png("~/projects/yeastOverdominance/manuscript/figures/multidrugHoepfner.png")
par(las=1,bty="l")
##hist(filteredORFs)
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
plot(h$mids+0.5,log10(h$counts),pch=19,type="n",ylab="log(# of strains)",xlab="# of experiments resistant",xlim=c(0,23)) #set axes etc
rect(h$mids,0,h$mids+1,logHeights,col="grey")
arrows(h$mids+.5,log10(ll),h$mids+.5,log10(ul),code=0,angle=90)
points(h$mids+0.5,log10(h$counts),pch=19,type="o")
dev.off()

  #remove ORFs with multidrug resistance as possibly carrying secondary mutations
notTrusted<-which(filteredORFs>10)

#second pass: dubious and random sets
#-----------------------------------
rssE <- numeric(0)
rssO <- numeric(0)
dssE <- numeric(0)
dssO <- numeric(0)
nssE <- numeric(0)
nssO <- numeric(0)
rssEZ <- numeric(0)
rssOZ <- numeric(0)
dssEZ <- numeric(0)
dssOZ <- numeric(0)
nssEZ <- numeric(0)
nssOZ <- numeric(0)

for(i in c(1:(length(m$V1)))){ #for each condition
    cmbid <- as.numeric(strsplit(as.character(m[i,1]),"_")[[1]][2]) #parse out compound name
  name<-n[which(cmbid==n[,1]),2]
    if(length(name)==0){
        titleString<-paste("CMB",cmbid,sep="")
    }else{
        titleString<-name
    }
 
    print(i)
    hipIndex <- m[i,2]+1  # HIP sensitivity
    hopIndex <- m[i,3]+1  # HOP sensitivity
    hipIndexZ <- m[i,4]+1 # HIP z-score
    hopIndexZ <- m[i,5]+1 # HOP z-score
  #calculate percentiles
    hipD<-quantile(hip[dI,hipIndex],probs=c(.05,.95),na.rm=T,names=F)
    hopD<-quantile(hop[dI,hopIndex],probs=c(.05,.95),na.rm=T,names=F)

    hipDZ<-quantile(hip[dI,hipIndexZ],probs=c(.05,.95),na.rm=T,names=F)
    hopDZ<-quantile(hop[dI,hopIndexZ],probs=c(.05,.95),na.rm=T,names=F)

    ds.homD  <- numeric(0)
    ds.homDZ <- numeric(0)
    ds.hetD  <- numeric(0)
    ds.hetDZ <- numeric(0)
    rs.homD  <- numeric(0)
    rs.homDZ <- numeric(0)
    rs.hetD  <- numeric(0)
    rs.hetDZ <- numeric(0)

    for(boot in c(1:1000)){
        dds <- sample(allDubious$V2,120) #dubious dubious ORFs
        dboot <- which(hop[,1] %in% dds)
        d.homD<-quantile(hop[dboot,hopIndex],probs=c(.05,.95),na.rm=T,names=F)
        d.homDZ<-quantile(hop[dboot,hopIndexZ],probs=c(.05,.95),na.rm=T,names=F)
        d.hetD<-quantile(hip[dboot,hipIndex],probs=c(.05,.95),na.rm=T,names=F)
        d.hetDZ<-quantile(hip[dboot,hipIndexZ],probs=c(.05,.95),na.rm=T,names=F)

        rds <- sample(length(hip[,1]),120) #random ORFs
        r.homD  <- quantile(hop[rds,hopIndex],probs=c(.05,.95),na.rm=T,names=F)
        r.homDZ <- quantile(hop[rds,hopIndexZ],probs=c(.05,.95),na.rm=T,names=F)
        r.hetD  <- quantile(hip[rds,hipIndex],probs=c(.05,.95),na.rm=T,names=F)
        r.hetDZ <- quantile(hip[rds,hipIndexZ],probs=c(.05,.95),na.rm=T,names=F)
        
        #save
        ds.homD  <- rbind(ds.homD,d.homD)
        ds.homDZ <- rbind(ds.homDZ,d.homDZ)
        ds.hetD  <- rbind(ds.hetD,d.hetD)
        ds.hetDZ <- rbind(ds.hetDZ,d.hetDZ)

        rs.homD  <- rbind(rs.homD, r.homD)
        rs.homDZ <- rbind(rs.homDZ,r.homDZ)
        rs.hetD  <- rbind(rs.hetD, r.hetD)
        rs.hetDZ <- rbind(rs.hetDZ,r.hetDZ)
    }

     #save all
     rssE <- rbind(rssE, colMeans(rs.hetD))
     rssO <- rbind(rssO, colMeans(rs.homD))
     dssE <- rbind(dssE, colMeans(ds.hetD))
     dssO <- rbind(dssO, colMeans(ds.homD))
     nssE <- rbind(nssE, hipD)
     nssO <- rbind(nssO, hopD)
     
     rssEZ <- rbind(rssEZ, colMeans(rs.hetDZ))
     rssOZ <- rbind(rssOZ, colMeans(rs.homDZ))
     dssEZ <- rbind(dssEZ, colMeans(ds.hetDZ))
     dssOZ <- rbind(dssOZ, colMeans(ds.homDZ))
     nssEZ <- rbind(nssEZ, hipDZ)
     nssOZ <- rbind(nssOZ, hopDZ)
 }


     #compare ranges
rrangeE <- rssE[,2] - rssE[,1]
drangeE <- dssE[,2] - dssE[,1]
nrangeE <- nssE[,2] - nssE[,1]
rrangeEZ <- rssEZ[,2] - rssEZ[,1]
drangeEZ <- dssEZ[,2] - dssEZ[,1]
nrangeEZ <- nssEZ[,2] - nssEZ[,1]

rrangeO <- rssO[,2] - rssO[,1]
drangeO <- dssO[,2] - dssO[,1]
nrangeO <- nssO[,2] - nssO[,1]
rrangeOZ <- rssOZ[,2] - rssOZ[,1]
drangeOZ <- dssOZ[,2] - dssOZ[,1]
nrangeOZ <- nssOZ[,2] - nssOZ[,1]

sum(nrangeE<drangeE)/277*100
sum(nrangeEZ<drangeEZ)/277*100
sum(nrangeO<drangeO)/277*100
sum(nrangeOZ<drangeOZ)/277*100

sum(nrangeE<rrangeE)/277*100
sum(nrangeEZ<rrangeEZ)/277*100
sum(nrangeO<rrangeO)/277*100
sum(nrangeOZ<rrangeOZ)/277*100

#compare upper CI
sum(nssE[,2]>rssE[,2])/277*100
sum(nssEZ[,2]>rssEZ[,2])/277*100
sum(nssE[,2]>dssE[,2])/277*100
sum(nssEZ[,2]>dssEZ[,2])/277*100

#third pass calculate proportion of extreme overdominance
#--------------------------------------------------------
eoL <- numeric(0) # number of extreme overdominant in each condition
hetBenL <- numeric(0)
eoD<- numeric(0) # number of extreme overdominant in each condition
hetBenD<- numeric(0)
eoR <- numeric(0) # number of extreme overdominant in each condition
hetBenR <- numeric(0)
condNames <- character(0)
for(counter in c(1:(length(m$V1)))){ #for each condition
    cmbid <- as.numeric(strsplit(as.character(m[counter,1]),"_")[[1]][2]) #parse out compound name
    name<-n[which(cmbid==n[,1]),2]
    if(length(name)==0){
        titleString<-paste("CMB",cmbid,sep="")
    }else{
        titleString<-name
    }
    condNames <- append(condNames,as.character(titleString))

    hipIndex <- m[counter,2]+1  # HIP sensitivity
    hopIndex <- m[counter,3]+1  # HOP sensitivity
    hipIndexZ <- m[counter,4]+1 # HIP z-score
    hopIndexZ <- m[counter,5]+1 # HOP z-score

                                        #find outliers for HIP experiments
    HetBen<-which(hip[,hipIndex]>nssE[counter,2]) # Row Index of strains with HIP larger than 95% 
    HetBenZ<-which(hip[,hipIndexZ]>nssEZ[counter,2]) # Row Index of strains with HIP larger than 95% for z-score
                                        #intersection
    HetBeneficial <- intersect(HetBen,HetBenZ)
                                        #remove ORFs with multidrug resistance as possibly carrying secondary mutations
    HetBeneficial <- setdiff(HetBeneficial,notTrusted)

                                        #find homozygously deleterious (outliers from HOP experiment)
    HomDel<-which(hop[,hopIndex]<nssO[counter,1])
    HomDelZ<-which(hop[,hopIndexZ]<nssOZ[counter,1])
    HomDeleterious <- intersect(HomDel,HomDelZ)

                                        # extreme overdominant
    eo <- intersect(HetBeneficial,HomDeleterious)
    write.table(hip[HetBeneficial,1],file=paste("Hoepfner/odORFs/",titleString,"hetBen.dat",sep=""),row.names=F,col.names=F)
    write.table(hip[eo,1],file=paste("Hoepfner/odORFs/",titleString,"eo.dat",sep=""),row.names=F,col.names=F)

    eoL <- append(eoL,length(eo))
    hetBenL <- append(hetBenL,length(HetBeneficial))

    pdf("temp.pdf")
    par(mfcol=c(2,2))
    plot(hip[,hipIndexZ],hip[,hipIndex],cex=0.5,col="#00000070",main=paste("HIP",titleString),pch=19)
    points(hip[HetBeneficial,hipIndexZ],hip[HetBeneficial,hipIndex],col="red",main=paste("HIP",titleString),pch=19)
    abline(h=nssE[counter,])
    abline(v=nssEZ[counter,])

    plot(hop[,hopIndexZ],hop[,hopIndex],cex=0.5,col="#00000070",pch=19)
    points(hop[HetBeneficial,hopIndexZ],hop[HetBeneficial,hopIndex],col="red",pch=19)
    abline(h=nssO[counter,])
    abline(v=nssOZ[counter,])
dev.off()
stop()
    #na do an brisko ta outliers pou entopizoun sto paper
    #na allaxo to match.pl oste na kanei match pio broadly kai na do an brisko ta outliers pou entopizoun sto paper
    
  #repeat for dubious and random ORF sets
  HetBen<-which(hip[,hipIndex]>dssE[counter,2])
  HetBenZ<-which(hip[,hipIndexZ]>dssEZ[counter,2])
  HetBeneficial <- intersect(HetBen,HetBenZ)
  HetBeneficial <- setdiff(HetBeneficial,notTrusted)
  HomDel<-which(hop[,hopIndex]<dssO[counter,1])
  HomDelZ<-which(hop[,hopIndexZ]<dssOZ[counter,1])
  HomDeleterious <- intersect(HomDel,HomDelZ)
  eo <- intersect(HetBeneficial,HomDeleterious)
  eoD <- append(eoD,length(eo))
  hetBenD <- append(hetBenD,length(HetBeneficial))

  HetBen<-which(hip[,hipIndex]>rssE[counter,2])
  HetBenZ<-which(hip[,hipIndexZ]>rssEZ[counter,2])
  HetBeneficial <- intersect(HetBen,HetBenZ)
  HetBeneficial <- setdiff(HetBeneficial,notTrusted)
  HomDel<-which(hop[,hopIndex]<rssO[counter,1])
  HomDelZ<-which(hop[,hopIndexZ]<rssOZ[counter,1])
  HomDeleterious <- intersect(HomDel,HomDelZ)
  eo <- intersect(HetBeneficial,HomDeleterious)
  eoR <- append(eoR,length(eo))
  hetBenR <- append(hetBenR,length(HetBeneficial))
  print(counter)
    
}

#compare across conditions
ho1 <- read.table("./Hoepfner/odORFs/ClotrimazolehetBen.dat")
hi1 <- read.table("./analysis/odORFs/CLOTRIMAZOLE20oldscannerhetBen.dat")
ho2 <- read.table("./Hoepfner/odORFs/MethotrexatehetBen.dat")
hi2 <- read.table("./analysis/odORFs/METHOTREXATE20oldscannerhetBen.dat")
ho3 <- read.table("./Hoepfner/odORFs/StaurosporinhetBen.dat")
hi3 <- read.table("./analysis/odORFs/STAUROSPORINE20newscannerhetBen.dat")

intersect(ho1$V1,hi1$V1)
intersect(ho2$V1,hi2$V1)
intersect(ho3$V1,hi3$V1)
ho1
hi1
counter<-173
hipIndex <- m[counter,2]+1  # HIP sensitivity
hopIndex <- m[counter,3]+1  # HOP sensitivity
hipIndexZ <- m[counter,4]+1 # HIP z-score
hopIndexZ <- m[counter,5]+1 # HOP z-score

hi1Index<- which(hip[,1]%in%hi1$V1)
ho1Index<- which(hip[,1]%in%ho1$V1)

png("temp.png")
plot(hip[,hipIndexZ],hip[,hipIndex],cex=0.5,col="#00000070")
points(hip[hi1Index,hipIndexZ],hip[hi1Index,hipIndex],col="#ff000099")
points(hip[ho1Index,hipIndexZ],hip[ho1Index,hipIndex],col="#0000ff99")

dev.off()

eoo1 <- read.table("./Hoepfner/odORFs/Clotrimazoleeo.dat")
eoi1 <- read.table("./analysis/odORFs/CLOTRIMAZOLE20oldscannereo.dat")
eoo2 <- read.table("./Hoepfner/odORFs/Methotrexateeo.dat")
eoi2 <- read.table("./analysis/odORFs/METHOTREXATE20oldscannereo.dat")
eoo3 <- read.table("./Hoepfner/odORFs/Staurosporineo.dat")
eoi3 <- read.table("./analysis/odORFs/STAUROSPORINE20newscannereo.dat")
eoo2
eoi2

png("temp.png",height=800,width=640)
par(las=1,bty="l",mai=c(.5,3,0,0.2),mgp=c(1.5,0.5,0),cex=0.8)
hb <- data.frame(condNames,eoL/hetBenL,eoL,hetBenL)
orderedhB <- hb[order(hb[,2]),]
orderedhB <- orderedhB[which(orderedhB[,2]>0),]
ci <- binom.confint(orderedhB[,3],orderedhB[,4],method="exact")
bar <- barplot(orderedhB[,2],names.arg=orderedhB[,1],horiz=TRUE,xlim=c(0,0.55),xlab="extreme overdominance")
arrows(ci$lower,bar,ci$upper,bar,code=0)
dev.off()

mean(eoL/hetBenL)
mean(eoD/hetBenD)
mean(eoR/hetBenR)

