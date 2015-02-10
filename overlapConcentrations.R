# find overlap in het ben in conditions with different concentrations in Hoepfner et al., 2014
load("~/projects/yeastOverdominance/collection/data/condID.dat")
load("~/projects/yeastOverdominance/collection/data/condNames.dat")
load("~/projects/yeastOverdominance/collection/data/concs.dat")
load("~/projects/yeastOverdominance/collection/data/nsum.dat")

df <- data.frame(condID,condNames,as.numeric(concs),nsum[,1:5]) #hetben, hetdel, homben, homdel

#compare overlap of heterozygous beneficial across compounds and across conditions
path <- "~/projects/yeastOverdominance/collection/Hoepfner/odORFs/"
files <- dir(path=path, pattern="*hetBen.dat")
bootstraps <- 10000
#find conditions for which we have more than one concentration
conditions <- character(0)
for(f in files){
  conditions <- append(conditions,strsplit(f,".",fixed=T)[[1]][1])
}

filesD <- which(duplicated(conditions))

# calculate overlap in all possible combinations of files with shared compound
#for each duplicate conditions
bcomSame <- numeric(0)
for(duplCond in unique(conditions[filesD])){
    pairs <- combn(which(conditions==duplCond),2) # find all combinations
                                        #for each combination
    for(combination in c(1:ncol(pairs))){
        pair <- pairs[,combination]
        fileA <- read.table(paste(path,files[pair[1]],sep=""))
        fileB <- read.table(paste(path,files[pair[2]],sep=""))
        bcomSame <- append(bcomSame,length(intersect(fileA[,1],fileB[,1]))/length(union(fileA[,1],fileB[,1])))
   }
}

#for all combinations combn(x,2)
# calculate overlap

                                        #perform 10000 bootstrap measurements, to measure overlap across measurements with the same compound
bcomDiff <- numeric(0)
bootNumber <- 10000
while(bootNumber > 0){
  pair <- sample(filesD,2) #pick a pair of files
  compoundA <- strsplit(files[pair[1]],".",fixed=T)[[1]][1]
  compoundB <- strsplit(files[pair[2]],".",fixed=T)[[1]][1]
  if(compoundA != compoundB){
      fileA <- read.table(paste(path,files[pair[1]],sep=""))
      fileB <- read.table(paste(path,files[pair[2]],sep=""))
      bcomDiff <- append(bcomDiff,length(intersect(fileA[,1],fileB[,1]))/length(union(fileA[,1],fileB[,1])))
      bootNumber <- bootNumber - 1
      print(bootNumber)
  }
}

overlapConcentrations <- list("bootstraps" = bootstraps,"bcomSame" = bcomSame,"bcomDiff" = bcomDiff)
save(overlapConcentrations,file="~/projects/yeastOverdominance/collection/data/overlapConcentrations.dat")
