# jackknife resampling of Dataset II

maxJ <- 1000 # number of jackknife resampling

# read required libraries

# read required data files

# load not trusted genes
notTrusted.Hil <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hil.dat")
notTrusted.Hop <- read.table("~/projects/yeastOverdominance/collection/data/notTrusted.Hop.dat")
therio <- union(notTrusted.Hil$V1,notTrusted.Hop$V1)

# for each condition
#   for maxJ
#     half <- sample 60 from the neutral set
#     rest 60
#     random 60 not in half
#     calculate hetBen, hetBenZ, homDel, homDelZ
#     exclude therio
#     calculate eo with rest (N) and rrest (R)
#     append in data structure
#   compare N and R
#adjust for multiple testing and report p.values
