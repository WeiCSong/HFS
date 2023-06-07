

library(R.utils)
library(data.table)
library(dplyr)
library(tidyr)

loci=fread("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/data/loci.bed",data.table=F)
base=fread("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/UKB/Table S1.csv",data.table=F)
base=base[which(is.na(base$nhap)),4]
loci=loci[which(loci[,4] %in% base),1:5]
loci=loci[order(loci[,2]),]

stat=data.frame()
for (i in 1:(nrow(loci))){  
  if(loci[i,1]=="chr22"){next()}
  file=paste(c(loci[i,1],loci[i,4],"sei"),collapse="/")
  if(!file.exists(file)){next()}
  sei=fread(file,data.table=F,sep=" ",fill=T)
  sei=separate(sei,V1,c("x","y"),sep="\t")
  colnames(sei)=paste("V",1:ncol(sei),sep="")

  sei[which(sei[,1]=="int"),1]=1 
  sei[]=apply(sei,2,as.numeric)
  sei=sei[order(sei[,1]),]
  sei=sei[match(1:sei[nrow(sei),1],sei[,1]),-1,drop=F] #just to avoid missing haplotype ids.
  rownames(sei)=1:nrow(sei)
  if(is.na(sei[1,1])){sei[1,]=sei[2,]}
  nhap=nrow(sei)
  class=colnames(sei)[which.max(sei[1,])]
  refvec=sei[1,]
  convert=length(which(apply(sei,1,which.max)!=which.max(sei[1,])))
  sei=sei[,class,drop=F]
  lower=min(sei)
  upper=max(sei)
  ref=sei[1,1]
  
  int=data.frame(locus=loci[i,4],class=class,nhap=nhap,convert=convert,ref=ref,lower=lower,upper=upper,refvec)
  stat=rbind(stat,int)
}

fwrite(stat,"HFSstat_miss.txt",sep="\t")



