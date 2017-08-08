library(tcR)
library(reshape2)
library(parallel)

getCdrs<-function(cls_list) melt(lapply(cls_list,"[",c(4,5,6,7,8)))

repGrep <- function(cls_list,tgs, nsubs = 1, minlen = 6,V = "family") {

  melted<- getCdrs(cls_list)
  melted<- melted[,c(2,7,1,3,4,6)]
  
  self<- F
  
  if (missing(tgs)) {
    tgs <- melted[,c(grep("CDR3.amino.acid.sequence", colnames(melted)),grep("V.gene", colnames(melted))) ]
    self<-T
  } else colnames(tgs)<-c ("CDR3.amino.acid.sequence","V.gene")
  
  tgs<-tgs[ nchar(tgs[,1]) > minlen,]
      
  i<- 0
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
 
  print(dim(tgs))
  
  match.list <- parSapply(cl = cl,
                          tgs$CDR3.amino.acid.sequence,
                          FUN = function(x){
                          i<<-i+1
                          
                          ind <-  agrep(pattern =x,                      #ищем индексы строк, совпадающих по аа сиквенсу с nsubs замен
                               x= as.character(melted[,1]),
                               value = F,
                               max = list(cost=1, all=1, sub=nsubs, del =0,ins=0) )
                          
                          #if (self) ind<-ind[ind!=i]    
                          
                          
                          
                            if  ( V == "family"| V == "exact" ) {         #сравниваем по V гену если указаны соотв. аргументы
                              V.genes <- unlist(strsplit(tgs$V.gene[i],", "))     
                                
                                if (V == "family")  V.genes <- sapply(V.genes,               #если нужно, обрезаем V ген до семейства
                                                                       function(x){substr(x, 
                                                                                          1, 
                                                                                          regexpr(pattern = "-",text = x)[1])
                                                                         }
                                                                      )
                                
                              m <- sapply(V.genes,                                    # ищем каждый из V.genes хотя бы одно пересечение среди V генов мэтча 
                                          function(x) grepl(pattern = x,                   
                                                            x = melted$V.gene[ind]))
                              ind<-ind[apply(m,1,FUN = any)]
                            } else return(ind)
                        
                          }
  )
  stopCluster(cl)
  match.list<-match.list[sapply(match.list,length)>1]
  melted.match.list<-melt(match.list)
  if (nrow(melted.match.list)==0) melted.match.list$target_cdr <- character(0)
  colnames(melted.match.list)[1]<-"match_cdr_num"
  colnames(melted.match.list)[2]<-"target_cdr"
  melted.match.list$match_cdr<- melted[melted.match.list[,1],1]
  melted.match.list$match_sample<- melted[melted.match.list[,1],2]
  melted.match.list$match_ncdr<- melted[melted.match.list[,1],3]
  melted.match.list$match_v<- melted[melted.match.list[,1],4]
  melted.match.list$match_j<- melted[melted.match.list[,1],5]
  melted.match.list$match_p<- melted[melted.match.list[,1],6]
  #melted.match.list$exact<- melted[melted.match.list[,1],2]

  melted.match.list<-melted.match.list[,c(2,1,3:ncol(melted.match.list))]
  
  return(melted.match.list)  ##[melted.match.list$ind>0,])
}


shared.counts<-function(cls_list,tgs,nsubs = 1,mod = "samples", minlen = 6,V = "family"){
 
  if(missing(tgs)){
    melted.match.list<-repGrep(cls_list,nsubs=nsubs,minlen=minlen, V=V)
  } else {
    melted.match.list<-repGrep(cls_list,tgs,nsubs=nsubs,minlen=minlen,V=V)
  }
  if (nrow(melted.match.list) == 0) {
    print(paste("Detected no clones with shared CDR3 amino-acid sequences with", nsubs, "substitution(s)"))
    return(0)
  }
  
  if (mod=="samples"){
    cdr.counts<-aggregate(melted.match.list, by = list(melted.match.list$match_sample, melted.match.list$target_cdr,melted.match.list$match_v),FUN = length)
    cdr.counts<-aggregate(cdr.counts,by = list(cdr.counts$Group.2), FUN = length)
    colnames(cdr.counts)<-c("cdr","number_of_samples")
  } else if (mod=="clones") {
    cdr.counts<-aggregate(melted.match.list, by = list(melted.match.list$target_cdr,melted.match.list$match_v),FUN = length)
    colnames(cdr.counts)<-c("cdr","number_of_clones")
  } 
  cdr.counts<-cdr.counts[,1:3]
  #cdr.counts<-cdr.counts[cdr.counts[,2]>1,]
  #cdr.counts<-cdr.counts[order(-cdr.counts[,2]),]
  return(cdr.counts)
}