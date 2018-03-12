library(dplyr)
library(data.table)
library(stringr)
library(magrittr)

#################### MIXCR TRB ###############

mixcr_TRB<- function (samples,outputFolder, new_meta_features,return_cloneset_list=T, threads=4){
  mixcr_file<-"/raid/users/Artem/.linuxbrew/bin/mixcr"
  if(str_sub(outputFolder,nchar(outputFolder),nchar(outputFolder)) == "/"){
    outputFolder<-str_sub(outputFolder,1,nchar(outputFolder)-1)
  }
  
  if (class(samples)=="character") {
    # if (str_sub(samples,start = str_length(samples)-3,end=str_length(samples)) ==".txt"){
    meta<-fread(samples,sep="\t",header = T)
    # } else {
    #   write("provided string is not a name of a .txt file, trying to read from folder", stderr())
    #   R1<-system(paste0("ls ", samples, "*R1*"),intern = T)
    #   meta<-data.table(R1=R1, R2=str_replace(R1,"R1","R2")) %>%
    #     filter(R2 %in% system(paste0("ls ", samples, "*R2*"),intern = T))
    # }
  } else if ("data.frame" %in% class(samples)) {
    meta<-samples
  }
  
  if(missing(new_meta_features)){
    new_meta_features<-colnames(meta)
  }
  
  system(paste0('mkdir -p ',outputFolder,'/alignments'))
  system(paste0('mkdir -p ',outputFolder,'/clonesets'))
  system(paste0('mkdir -p ',outputFolder,'/export'))
  
  mapply(FUN = function(id,r1,r2){
    system(paste0(mixcr_file,
                  ' align -f -OvParameters.geneFeatureToAlign=VTranscript --report ',
                  outputFolder,"/alignments/", id,'_report.log ',
                  "-t ", threads, " ",
                  r1," ",r2," ",
                  outputFolder,"/alignments/", id,".vdjca" ))
  },
  id=meta$id,r1=meta$R1,r2=meta$R2
  )
  
  sapply(meta$id, function(id){
    system(paste0(mixcr_file,
                  ' assemble -f -OseparateByV=true -OaddReadsCountOnClustering=true --report ',
                  outputFolder,"/clonesets/", id,'_report.log ',
                  "-t ", threads, " ",
                  outputFolder,"/alignments/", id,".vdjca ",
                  outputFolder,"/clonesets/", id,".clns"))
  })
  
  sapply(meta$id, function(id){
    system(paste0(mixcr_file,
                  ' exportClones -f -s -c TRB -o -t -count -fraction -nFeature CDR3 -aaFeature CDR3 -vGene -dGene -jGene -positionOf VEndTrimmed -positionOf DBeginTrimmed -positionOf DEndTrimmed -positionOf JBeginTrimmed ',
                  outputFolder,"/clonesets/", id,".clns ",
                  outputFolder,"/export/", id,".txt"))
  })
  
  clonesets_list<- lapply(meta$id, function(x) {
    cloneset<-fread (paste0(outputFolder,"/export/", x,".txt"),
                     sep = "\t",
                     col.names=c("count",
                                 "frequency", 
                                 "CDR3nt",
                                 "CDR3aa",
                                 "V","D","J",
                                 "Vend",
                                 "Dstart",
                                 "Dend",
                                 "Jstart") )
    fwrite(cloneset,paste0(outputFolder,"/export/", x,".txt"),quote = F,sep="\t",col.names = T)
    return(cloneset)
  })
  meta$file<-paste0(outputFolder,"/export/", meta$id,".txt")
  names(clonesets_list)<-meta$id
  
  new_meta_features<-new_meta_features[new_meta_features %in% colnames(meta)]
  
  meta%<>% select(c("id","file",new_meta_features),-R1,-R2)
  fwrite(meta,paste0(outputFolder,"/meta_clonesets.txt"),sep="\t")
  if(return_cloneset_list){
    return(clonesets_list)
  }
}




mixcr_RNAseq<- function (samples,outputFolder, new_meta_features,return_cloneset_list=T){
  mixcr_file<-"/raid/users/Artem/.linuxbrew/bin/mixcr"
  if(str_sub(outputFolder,nchar(outputFolder),nchar(outputFolder)) == "/"){
    outputFolder<-str_sub(outputFolder,1,nchar(outputFolder)-1)
  }
  
  if (class(samples)=="character") {
    # if (str_sub(samples,start = str_length(samples)-3,end=str_length(samples)) ==".txt"){
    meta<-fread(samples,sep="\t",header = T)
    # } else {
    #   write("provided string is not a name of a .txt file, trying to read from folder", stderr())
    #   R1<-system(paste0("ls ", samples, "*R1*"),intern = T)
    #   meta<-data.table(R1=R1, R2=str_replace(R1,"R1","R2")) %>%
    #     filter(R2 %in% system(paste0("ls ", samples, "*R2*"),intern = T))
    # }
  } else if ("data.frame" %in% class(samples)) {
    meta<-samples
  }
  
  if(missing(new_meta_features)){
    new_meta_features<-colnames(meta)
  }
  
  system(paste0('mkdir -p ',outputFolder,'alignments'))
  system(paste0('mkdir -p ',outputFolder,'clonesets'))
  system(paste0('mkdir -p ',outputFolder,'export'))
  
  mapply(FUN = function(id,r1,r2){
    system(paste0(mixcr_file,
                  ' -p rna-seq -OallowPartialAlignments=true -f --report ',
                  outputFolder,"/alignments/", id,'_report.log ',
                  r1," ",r2," ",
                  outputFolder,"/alignments/", id,".vdjca" ))
  },
  id=meta$id,r1=meta$R1,r2=meta$R2
  )
  
  sapply(meta$id, function(id){
    system(paste0(mixcr_file,
                  ' assemblePartial --report ',
                  outputFolder,"/alignments/", id,'_report.log ',
                  outputFolder,"/alignments/", id,".vdjca ",
                  outputFolder,"/alignments/", id,"Rescued_1.vdjca"))
  })
  
  sapply(meta$id, function(id){
    system(paste0(mixcr_file,
                  ' assemblePartial --report ',
                  outputFolder,"/alignments/", id,'_report.log ',
                  outputFolder,"/alignments/", id,"Rescued_1.vdjca ",
                  outputFolder,"/alignments/", id,"Rescued_2.vdjca"))
  })
  
  sapply(meta$id, function(id){
    system(paste0(mixcr_file,
                  ' assemble -f -OseparateByV=true -OaddReadsCountOnClustering=true --report ',
                  outputFolder,"/clonesets/", id,'_report.log ',
                  outputFolder,"/alignments/", id,".vdjca ",
                  outputFolder,"/clonesets/", id,".clns"))
  })
  
  sapply(meta$id, function(id){
    system(paste0(mixcr_file,
                  ' exportClones -f -s -c TRB -o -t -count -fraction -nFeature CDR3 -aaFeature CDR3 -vGene -dGene -jGene -positionOf VEndTrimmed -positionOf DBeginTrimmed -positionOf DEndTrimmed -positionOf JBeginTrimmed ',
                  outputFolder,"/clonesets/", id,".clns ",
                  outputFolder,"/export/", id,".txt"))
  })
  
  clonesets_list<- lapply(meta$id, function(x) {
    cloneset<-fread (paste0(outputFolder,"/export/", x,".txt"),
                     sep = "\t",
                     col.names=c("count",
                                 "frequency", 
                                 "CDR3nt",
                                 "CDR3aa",
                                 "V","D","J",
                                 "Vend",
                                 "Dstart",
                                 "Dend",
                                 "Jstart") )
    fwrite(cloneset,paste0(outputFolder,"/export/", x,".txt"),quote = F,sep="\t",col.names = T)
    return(cloneset)
  })
  meta$file<-paste0(outputFolder,"/export/", meta$id,".txt")
  names(clonesets_list)<-meta$id
  
  new_meta_features<-new_meta_features[new_meta_features %in% colnames(meta)]
  
  meta%<>% select(c("id","file",new_meta_features),-R1,-R2)
  fwrite(meta,paste0(outputFolder,"/meta_clonesets.txt"),sep="\t")
  if(return_cloneset_list){
    return(clonesets_list)
  }
}
