#calculate correlations across barcodes
#read all barcodes for each condition
#find common genes
#calculate pairwise correlation
library(stringr)
notZero <- function(x){
    if(length(x) == 0){
        return("NA")
    }else{
        return(x)
    }
}

conditions <- read.table("~/projects/yeastOverdominance/collection/analysis/conditions.dat")

allHomCors <- matrix(nrow=length(conditions$V1),ncol=6)
allHetCors <- matrix(nrow=length(conditions$V1),ncol=6)
allHomPs <- matrix(nrow=length(conditions$V1),ncol=6)
allHetPs <- matrix(nrow=length(conditions$V1),ncol=6)
conditionCounter <- 1
for(i in conditions$V1){
    print(i)
    ua <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ua.dat",sep=""),header=F)
    da <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".da.dat",sep=""),header=F)
    us <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".us.dat",sep=""),header=F)
    ds <- read.table(paste("~/projects/yeastOverdominance/collection/analysis/fitComp/",i,".ds.dat",sep=""),header=F)
#str_extract(ds$V1,"") 
    ua$V1 <- str_extract(ua$V1,"[^:]*::[^:]*")
    us$V1 <- str_extract(us$V1,"[^:]*::[^:]*")
    da$V1 <- str_extract(da$V1,"[^:]*::[^:]*")
    ds$V1 <- str_extract(ds$V1,"[^:]*::[^:]*") 

    ## ua$V1 <- str_extract(ua$V1,"[^:]*")
    ## us$V1 <- str_extract(us$V1,"[^:]*")
    ## da$V1 <- str_extract(da$V1,"[^:]*")
    ## ds$V1 <- str_extract(ds$V1,"[^:]*") 
    allORFs <- unique(c(ua$V1,da$V1,us$V1,ds$V1))
    hets <- matrix(nrow=length(allORFs),ncol=4)
    homs <- matrix(nrow=length(allORFs),ncol=4)
    counter <- 1
    for(j in allORFs){
        uaSlice <- ua[which(ua$V1==j),c(2:3)] #hom, het
        daSlice <- da[which(da$V1==j),c(2:3)]
        usSlice <- us[which(us$V1==j),c(2:3)]
        dsSlice <- ds[which(ds$V1==j),c(2:3)]
                                        #ua, da, us, ds
        homs[counter,1] <- notZero(uaSlice$V2)
        hets[counter,1] <- notZero(uaSlice$V3)
        homs[counter,2] <- notZero(daSlice$V2)
        hets[counter,2] <- notZero(daSlice$V3)
        homs[counter,3] <- notZero(usSlice$V2)
        hets[counter,3] <- notZero(usSlice$V3)
        homs[counter,4] <- notZero(dsSlice$V2)
        hets[counter,4] <- notZero(dsSlice$V3)
        counter <- counter + 1
    }
    mode(homs)<-"numeric" #make numeric
    mode(hets)<-"numeric"
    hom1 <- cor.test(homs[,1],homs[,2])
    hom2 <- cor.test(homs[,1],homs[,3])
    hom3 <- cor.test(homs[,1],homs[,4])
    hom4 <- cor.test(homs[,2],homs[,3])
    hom5 <- cor.test(homs[,1],homs[,4])
    hom6 <- cor.test(homs[,3],homs[,4])
    het1 <- cor.test(hets[,1],hets[,2])
    het2 <- cor.test(hets[,1],hets[,3])
    het3 <- cor.test(hets[,1],hets[,4])
    het4 <- cor.test(hets[,2],hets[,3])
    het5 <- cor.test(hets[,1],hets[,4])
    het6 <- cor.test(hets[,3],hets[,4])
    allHomCors[conditionCounter,] <- c(hom1$estimate,hom2$estimate,hom3$estimate,hom4$estimate,hom5$estimate,hom6$estimate)
    allHetCors[conditionCounter,] <- c(het1$estimate,het2$estimate,het3$estimate,het4$estimate,het5$estimate,het6$estimate)
    allHomPs[conditionCounter,] <- c(hom1$p.value,hom2$p.value,hom3$p.value,hom4$p.value,hom5$p.value,hom6$p.value)
    allHetPs[conditionCounter,] <- c(het1$p.value,het2$p.value,het3$p.value,het4$p.value,het5$p.value,het6$p.value)
    ## hom.cors <- cor(homs,use="complete.obs")              #save all pairwise correlations
    ## het.cors <- cor(hets,use="complete.obs")
    ## mask <- upper.tri(hom.cors,diag=F)
    ## allHomCors[conditionCounter,] <- hom.cors[mask]
    ## allHetCors[conditionCounter,] <- het.cors[mask]
    conditionCounter <- conditionCounter + 1
}
homMean <- mean(allHomCors)
#[1] 0.6375129
homSD <- sd(allHomCors)
#[1] 0.1744695
hetMean <- mean(allHetCors)
#[1] 0.4141954
hetSD <- sd(allHetCors)
#[1] 0.2007997

#find how many are significantly positive
homSign <- sum(p.adjust(allHomPs)<0.05)*100/length(allHomPs)
hetSign <- sum(p.adjust(allHetPs)<0.05)*100/length(allHetPs)

#save
write.table(list(homMean = homMean, homSD = homSD, hetMean = hetMean, hetSD = hetSD, homSign = homSign, hetSign = hetSign),file="~/projects/yeastOverdominance/collection/data/Hil.barCor.dat",quote=FALSE,row.names=FALSE, col.names=TRUE)

