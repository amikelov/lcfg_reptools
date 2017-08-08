library(tcR)
library(reshape2)
source("repGrep.R")

# rn1 = 54456 - большой RNA-seq
hc.rn1<- parse.folder("data/E-GEOD-54456/TCRBs/hc/","mixcr")
ls.rn1<- parse.folder("data/E-GEOD-54456/TCRBs/ps/","mixcr")

# rn2 = 47944 - маленкьий RNA-seq
hc.rn2<- parse.folder("data/E-GEOD-47944/TCRBs/hc/","mixcr")
ls.rn2<- parse.folder("data/E-GEOD-47944/TCRBs/ls/","mixcr")
nls.rn2<- parse.folder("data/E-GEOD-47944/TCRBs/nls/","mixcr")

# rep = Harden 2015
hc.rep<- parse.folder("data/Harden/beta/hc/","immunoseq3")
ls.rep<- parse.folder("data/Harden/beta/ls/","immunoseq3")
nls.rep<- parse.folder("data/Harden/beta/nonls/","immunoseq3")

##basic statistics - number of unique clones derived from RNA-seq
mean(sapply(hc.rn1,function(x){length(x[,1])}))
mean(sapply(ls.rn1,function(x){length(x[,1])}))

mean(sapply(hc.rn2,function(x){length(x[,1])}))
mean(sapply(ls.rn2,function(x){length(x[,1])}))
mean(sapply(nls.rn2,function(x){length(x[,1])}))

#### intersections:

## shared inside subsets:
shared.hc.rn<-shared.counts(c(hc.rn1,hc.rn2), V="family")
shared.ls.rn<-shared.counts(c(ls.rn1,ls.rn2), V="family")

shared.hc.rn2<-shared.counts(hc.rn2)
shared.ls.rn2<-shared.counts(ls.rn2)
shared.nls.rn2<-shared.counts(nls.rn2)





## get fuzzy matches from harden
fromls<-repGrep(c(ls.rn1,ls.rn2), V="family")[,c(3,6)]


lsinharden <- repGrep(ls.rep,tgs = fromls,V= "family")
nlsinharden <- repGrep(nls.rep,tgs = fromls,V="family")


# get counts for fuzzy matches samples
shared.counts(nls.rep,as.character(shared.ls.rn1[,1]))
shared.counts(ls.rep,as.character(shared.ls.rn1[,1]))

# get counts for fuzzy matches clones
shared.counts(ls.rep,as.character(shared.ls.rn1[,1]),mod="clones")
shared.counts(nls.rep,as.character(shared.ls.rn1[,1]),mod="clones")


## calculate distance matrices
d<-lapply(unique(lsinharden$target_cdr), FUN = function(x) {stringdistmatrix(lsinharden[lsinharden$target_cdr==x,]$match_ncdr,method = "hamming", nthread =20 )})
names(d)<-unique(lsinharden$target_cdr)




ageing <- parse.folder("~/datasets/aging","mitcr")
public_match_cdr <- lapply(ageing,FUN=function(x) merge(x= data.frame(match_cdr=unique(lsinharden$match_cdr)),y= x, by.x = "match_cdr",by.y="CDR3.amino.acid.sequence", all.x =F, all.y =F))
shared.counts(ageing, as.character(shared.hc.rn1[,1]))

