#find overlap in hetBen in conditions with multiple timepoints from Hilenmeyer et al., 2008
#calculate overlap
path <- "~/projects/yeastOverdominance/collection/analysis/odORFs/"
bootstraps <- 10000
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
bootNumber <- bootstraps
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
bootNumber <- bootstraps
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
overlapGenerations <- list("bootstraps" = bootstraps,"bcomSame" = bcomSame,"bcomDiff" = bcomDiff) #result
save(overlapGenerations,file="~/projects/yeastOverdominance/collection/data/overlapGenerations.dat")
