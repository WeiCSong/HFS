library(R.utils)
library(data.table)
library(susieR)
library(doParallel)
library(arrow)
library(xgboost)
library(tidyr)
library(dplyr)
library(Rfast)
#ind=fread(ind,data.table=F)
#ind=ind[,1]
#fwrite(as.list(c("iid",ind)),"efile",sep=" ",col.names=F)
#d=d[which(d[,6]<0.01),]

#pheno=fread("~/FLAT/data/pheno.txt",data.table=F)
#ind=match(ind,pheno[,1])
#locilist=list.dirs(recursive=F)
#locilist=gsub("./","",locilist)

#for (loci in locilist){
#  setwd(loci)
#  if(!file.exists("sei")){
#    setwd("..")
#    next()
#  }
#  system("sed 's/\t/ /g' sei | awk 'length($0)<300'   > sei1")
#  sei=fread("sei1",data.table=F,fill=TRUE,header=F)
#  unlink("sei1")
#  hapsgeno=fread("hapsgeno.gz",data.table=F,fill=TRUE)
#  miss=502410-nrow(hapsgeno)
#  if(miss>0){hapsgeno=rbind(data.frame(V1=rep(NA,miss),V2=rep(NA,miss)),hapsgeno)} 
#  sei[which(sei[,1]=="int"),1]=1
#  sei[,1]=as.numeric(sei[,1])
#  sei=sei[order(sei[,1]),]
#  if(sei[1,1]!=1){
#    x=sei[1,]
#    x[1,1]=1
#    sei=rbind(x,sei)
#  }
#  sei=sei[match(1:sei[nrow(sei),1],sei[,1]),-1] #just to avoid missing haplotype ids.
#  rownames(sei)=1:nrow(sei)
#  sei=sei[,which.max(sei[1,])]
#  m1=sei[hapsgeno[,1]]
#  m2=sei[hapsgeno[,2]]
#  m=(m1+m2)/2
#  m=c(loci,m[ind])
#  setwd("..")
#  fwrite(as.list(m),"efile",sep=" ",append=T,col.names=F,na="NA")
#}

#system("~/gcta/osca --tefile efile --gene-expression --make-bod --thread-num 4 --no-fid --out hfs")

loci=fread("~/FLAT/data/loci.bed",data.table=F)
lociinfo=loci[,1:5]
pheno=fread("~/FLAT/data/pheno.txt",data.table=F)



get_locus_value=function(x){
  loci=gsub(" ","",x[1])
  var=gsub(" ","",x[2])
  chr=lociinfo[which(lociinfo[,4]==loci),1]
  print(loci)
  setwd(paste(chr,loci,sep="/"))
  
  system("sed 's/\t/ /g' sei | awk 'length($0)<300'   > sei1")
  sei=fread("sei1",data.table=F,fill=TRUE,header=F)
  unlink("sei1")
  colnames(sei)=paste("V",1:ncol(sei),sep="")
  hapsgeno=fread("hapsgeno.gz",data.table=F,fill=TRUE)
  miss=502410-nrow(hapsgeno)
  if(miss>0){hapsgeno=rbind(data.frame(V1=rep(NA,miss),V2=rep(NA,miss)),hapsgeno)}  #fread ignore blank lines at the beginning, 

  VAR=var
  if (grepl("max", var, fixed = TRUE)){
    type="max"
    VAR=gsub("max","",var)
  } else if (grepl("min", var, fixed = TRUE)) {
    type="min"
    VAR=gsub("min","",var)
  } else {type="add"}   
  
  sei=sei[,c("V1",VAR)]
  ref=sei[which(sei[,1]=="int")[1],2]
  if(is.na(ref)){m=c()}else{
    sei[which(sei[,1]=="int"),1]=1 
    sei[]=apply(sei,2,as.numeric)
    sei=sei[order(sei[,1]),]
    sei=sei[match(1:sei[nrow(sei),1],sei[,1]),-1,drop=F] #just to avoid missing haplotype ids.
    rownames(sei)=1:nrow(sei)
    m1=sei[hapsgeno[,1],,drop=F]
    m2=sei[hapsgeno[,2],,drop=F]
    m1[which(is.na(m1[,1])),1]=ref
    m2[which(is.na(m2[,1])),1]=ref
    m1=m1-ref
    m2=m2-ref
    if(type=="max"){
      m=Reduce(pmax, list(m1,m2))
    } else if (type=="min") {
      m=Reduce(pmin, list(m1,m2)) 
    } else {m=(m1+m2)/2}   
    colnames(m)=loci
    rownames(m)=1:nrow(m)
  }
  setwd("../..")
  return(m)
  gc()
}

  d=fread("BMD.summary",data.table=F)
  d=d[order(d$PIPfull,decreasing=T),]
  intd=d[1:2500,1:2]
  intd=as.list(as.data.frame(t(intd)))
  dat=mclapply(intd,get_locus_value,mc.cores=64)
  dat=dplyr::bind_cols(dat)
  rownames(dat)=pheno[,1]
  write_feather(data.frame(dat),"BMD.0.4.feather",compression="zstd")
  
  d=fread("BMI.summary",data.table=F)
  d=d[order(d$PIPfull,decreasing=T),]
  intd=d[1:2500,1:2]
  intd=as.list(as.data.frame(t(intd)))
  dat=mclapply(intd,get_locus_value,mc.cores=64)
  dat=dplyr::bind_cols(dat)
  rownames(dat)=pheno[,1]
  write_feather(data.frame(dat),"BMI.0.4.feather",compression="zstd") 
  
  d=fread("FVC.summary",data.table=F)
  d=d[order(d$PIPfull,decreasing=T),]
  intd=d[1:2500,1:2]
  intd=as.list(as.data.frame(t(intd)))
  dat=mclapply(intd,get_locus_value,mc.cores=64)
  dat=dplyr::bind_cols(dat)
  rownames(dat)=pheno[,1]
  write_feather(data.frame(dat),"FVC.0.4.feather",compression="zstd")  
   