library('R.matlab')

#Inserir pasta
setwd('')

folder<-'DATA_sep'
files<-list.files(folder)
data<-list()

for (i in files){
  parts<-unlist(strsplit(unlist(strsplit(i,'.mat')),'__'))
  if(length(parts)==2){
    data[[parts[1]]][[parts[2]]]<-readMat(paste(folder,i,sep='/'))$fi
  }
  else if(length(parts)==3){
    data[[parts[1]]][[parts[2]]][[parts[3]]]<-readMat(paste(folder,i,sep='/'))$fi
  }
  else{
    data[[parts[1]]][[parts[2]]][[parts[3]]][[parts[4]]]<-readMat(paste(folder,i,sep='/'))$fi
  }
}

saveRDS(data,file='MetroAreasData.rds')
asa<-readRDS('MetroAreasData.R')
