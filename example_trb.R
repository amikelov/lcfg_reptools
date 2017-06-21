library(tcR)
library(reshape2)
source("repGrep.R")

hc.54456<- parse.folder("data/E-GEOD-54456/TCRBs/hc/","mixcr")
ls.54456<- parse.folder("data/E-GEOD-54456/TCRBs/ps/","mixcr")

##basic statistics - number of unique clones derived from RNA-seq
mean(sapply(hc.54456,function(x){length(x[,1])}))
mean(sapply(ls.54456,function(x){length(x[,1])}))

#### intersections:

## shared inside subsets:
shared.hc.54456<-shared.counts(hc.54456)
shared.ls.54456<-shared.counts(ls.54456)
