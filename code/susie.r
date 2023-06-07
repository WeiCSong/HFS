param <- commandArgs(trailingOnly=T)
trait= eval(paste(text=param[1]))
chr=eval(paste(text=param[2]))
dir.create(trait)

library(R.utils)
library(data.table)
library(susieR)
library(doParallel)
library(arrow)
library(xgboost)
library(tidyr)
library(dplyr)
library(Rfast)

#cl <- makeCluster(4)
#registerDoParallel(cl)



d=fread("d",data.table=F)
d=d[grep("V",d[,3]),]
#d=d[which(d[,6]<0.01),]
loci=fread(paste("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/UKB/int/",chr,sep=""),data.table=F)
lociinfo=loci[,1:5]
pheno=fread("~/FLAT/data/pheno.txt",data.table=F)

covar=colnames(pheno)[2:16]
cdat=pheno[,c(trait,covar)]
rownames(cdat)=pheno[,1]

covarmod=lm(paste(trait,"~.",sep=""),data=cdat)
resid=resid(covarmod)
resid=na.omit(resid)

get_locus_value=function(x){
  loci=gsub(" ","",x[1])
  var=gsub(" ","",x[2])
  print(loci)
  setwd(loci)
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
  setwd("..")
  return(m)
  gc()
}
#imputecol=function(x){
#  x[which(is.na(x))]=median(x,na.rm=T)
#  return(x)
#}


susie_block=function(trait,block){
  locilist=lociinfo[which(lociinfo[,5]==block),4]
  intd=d[which(d[,1]==trait & d[,2] %in% locilist),2:3]
  intd=as.list(as.data.frame(t(intd)))
  dat=mclapply(intd,get_locus_value,mc.cores=20)
  dat=dplyr::bind_cols(dat)
  rownames(dat)=pheno[,1]
  l=which(pheno$group=="Caucasian" & pheno[,1] %in% names(resid))
  mat=dat[l,]
  #mat=as(apply(mat,2,imputecol),"matrix")
	fitted <- susie(data.matrix(mat),resid[rownames(mat)],L = 10,verbose = FALSE,coverage=0.3)  #has Rfast really used multithreading? Not sure how to check.
  coef=susie_get_posterior_mean(fitted)
  pip=fitted$pip
  pip=data.frame(pip=pip,beta=coef)
  dat=data.matrix(dat) %*% coef
  write_feather(data.frame(dat),paste(c(trait,"/",block,".","feather"),collapse=""),compression="zstd")
  return(pip)
  gc()
}

susieres=c()
for (block in unique(lociinfo[,5])){
  int=susie_block(trait=trait,block=block)
  susieres=rbind(susieres,int)
}
save(susieres,file=paste(trait,"susieres.RData",sep="/"))



