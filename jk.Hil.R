# jackknife resampling of Dataset I

# read required libraries
library(stringr)
# read required files
conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")
dubious<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
YPDdel <- read.csv("~/projects/yeastOverdominance/collection/Deutschbauer_etal/TableS2.csv")
YPDben <- read.csv("~/projects/yeastOverdominance/collection/Sliwa_Korona/adaptive",header=FALSE)

minBarNo <- 1 # >minBarNo should be available for each gene 
maxJ <- 1000 # number of jackknife resampling

notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)

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
        
                                        #filter out not Trusted genes
        hetBen <- setdiff(hetBen,therio)
        eo <- setdiff(eo,therio)
        
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

# save
write.table(data.frame(conditions=conditions$V1,Weo,Whb),file="~/projects/yeastOverdominance/collection/data/jk.Hil.dat",quote=FALSE,row.names=FALSE, col.names=FALSE)
