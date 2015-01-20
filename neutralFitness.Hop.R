#calculate and save neutral fitness bootstrap values for Hop data

#read data files
hip<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HIP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5912)))
hop<-read.csv("~/projects/yeastOverdominance/collection/Hoepfner/HOP_scores.txt",sep="\t",nrows=6681,comment.char="",colClasses=c("character",rep("numeric",5846)))
allDubious <-read.table("~/projects/yeastOverdominance/collection/sgd/dubious.tsv",sep="\t") # all dubious ORFs
d<-read.table("~/projects/yeastOverdominance/collection/Agrawal.Whitlock/dubious120.csv",header=T)

notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)
neutral <- d[which(!d[,1]%in%therio),1] # neutral set without not trusted
dI<-which(hop[,1]%in%neutral) # neutrals
#neutral set without not trusted

m<-read.table("~/projects/yeastOverdominance/collection/Hoepfner/groups.dat")
unq <- unique(paste(m$V1,m$V2)) #unique combinations of compound and concentration

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

save(nss,file="~/projects/yeastOverdominance/collection/data/nss.dat")
save(dss,file="~/projects/yeastOverdominance/collection/data/dss.dat")
save(rss,file="~/projects/yeastOverdominance/collection/data/rss.dat")


