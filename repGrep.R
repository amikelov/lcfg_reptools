getCdrs<-function(cls_list)melt(lapply(cls_list,"[[",6))

repGrep <- function(cls_list,tgs, nsubs) {
  
  melted<- getCdrs(cls_list)
  if (missing(tgs)) tgs <- as.character(melted[,1])
  if (missing(nsubs)) nsubs <-1
  match.list <- sapply(tgs,
                       FUN = function(x){
                         agrep(pattern =x,
                               x= melted[,1],
                               value = F,
                               max = list(cost=1, all=1, sub=nsubs, del =0,ins=0)
                               )
                       }
  )
  
  match.list<-match.list[sapply(match.list,length)>1]
  melted.match.list<-melt(match.list)
  colnames(melted.match.list)[1]<-"match_cdr_num"
  colnames(melted.match.list)[2]<-"target_cdr"
  melted.match.list$match_cdr<- melted[melted.match.list[,1],1]
  melted.match.list$match_sample<- melted[melted.match.list[,1],2]
  
  melted.match.list$exact<- melted[melted.match.list[,1],2]
  
  melted.match.list<-melted.match.list[,c(2,1,3,4)]
  
  return(melted.match.list)  ##[melted.match.list$ind>0,])
}


shared.counts<-function(cls_list,tgs,nsubs,mod){
  if(missing(mod)) mod <-"samples"
  if(missing(nsubs)) nsubs<-1
  if(missing(tgs)){
    melted.match.list<-repGrep(cls_list,nsubs=nsubs)
  } else {
    melted.match.list<-repGrep(cls_list,tgs,nsubs=nsubs)
  }
  
  if (mod=="samples"){
    cdr.counts<-aggregate(melted.match.list, by = list(melted.match.list$match_sample, melted.match.list$target_cdr),FUN = length)
    cdr.counts<-aggregate(cdr.counts,by = list(cdr.counts$Group.2), FUN = length)
    colnames(cdr.counts)<-c("cdr","number_of_samples")
  } else{
    cdr.counts<-aggregate(melted.match.list, by = list(melted.match.list$target_cdr),FUN = length)
    colnames(cdr.counts)<-c("cdr","number_of_clones")
  }
  cdr.counts<-cdr.counts[,1:2]
  cdr.counts<-cdr.counts[cdr.counts[,2]>1,]
  cdr.counts<-cdr.counts[order(-cdr.counts[,2]),]
  return(cdr.counts)
}