#read required libraries
library(stringr)

#read required files
conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")
dubious<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs

# remove from dubious not trusted

notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)
neutral <- dubious[which(!dubious[,1] %in% therio),1] # neutral set without not trusted

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
    ua.ind<-which(ua$V1%in%neutral)
    da.ind<-which(da$V1%in%neutral)
    us.ind<-which(us$V1%in%neutral)
    ds.ind<-which(ds$V1%in%neutral)

                                        # confidence intervals from the manually curated set of 120 dubious ORFs
    ua.ind <- which(ua$V1%in%neutral)
    da.ind <- which(da$V1%in%neutral)
    us.ind <- which(us$V1%in%neutral)
    ds.ind <- which(ds$V1%in%neutral)

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

                                        # save data in files
save(uahet,file="~/projects/yeastOverdominance/collection/data/uahet.dat")
save(dahet,file="~/projects/yeastOverdominance/collection/data/dahet.dat")
save(ushet,file="~/projects/yeastOverdominance/collection/data/ushet.dat")
save(dshet,file="~/projects/yeastOverdominance/collection/data/dshet.dat")
                 
save(duahet,file="~/projects/yeastOverdominance/collection/data/duahet.dat")
save(ddahet,file="~/projects/yeastOverdominance/collection/data/ddahet.dat")
save(dushet,file="~/projects/yeastOverdominance/collection/data/dushet.dat")
save(ddshet,file="~/projects/yeastOverdominance/collection/data/ddshet.dat")
                 
save(ruahet,file="~/projects/yeastOverdominance/collection/data/ruahet.dat")
save(rdahet,file="~/projects/yeastOverdominance/collection/data/rdahet.dat")
save(rushet,file="~/projects/yeastOverdominance/collection/data/rushet.dat")
save(rdshet,file="~/projects/yeastOverdominance/collection/data/rdshet.dat")

save(uahom,file="~/projects/yeastOverdominance/collection/data/uahom.dat")
save(dahom,file="~/projects/yeastOverdominance/collection/data/dahom.dat")
save(ushom,file="~/projects/yeastOverdominance/collection/data/ushom.dat")
save(dshom,file="~/projects/yeastOverdominance/collection/data/dshom.dat")
                 
save(duahom,file="~/projects/yeastOverdominance/collection/data/duahom.dat")
save(ddahom,file="~/projects/yeastOverdominance/collection/data/ddahom.dat")
save(dushom,file="~/projects/yeastOverdominance/collection/data/dushom.dat")
save(ddshom,file="~/projects/yeastOverdominance/collection/data/ddshom.dat")
                 
save(ruahom,file="~/projects/yeastOverdominance/collection/data/ruahom.dat")
save(rdahom,file="~/projects/yeastOverdominance/collection/data/rdahom.dat")
save(rushom,file="~/projects/yeastOverdominance/collection/data/rushom.dat")
save(rdshom,file="~/projects/yeastOverdominance/collection/data/rdshom.dat")

