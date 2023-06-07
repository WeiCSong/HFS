param <- commandArgs(trailingOnly=T)
path= eval(paste(text=param[1]))
chr=eval(paste(text=param[2]))
trait=eval(paste(text=param[3]))
h2=eval(paste(text=param[4]))
n=eval(paste(text=param[5]))
ncore=eval(paste(text=param[6]))
h2=as.numeric(h2)
n=as.numeric(n)
ncore=as.numeric(ncore)

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

setwd(paste("~/FLAT/UKB/",chr,sep=""))
d=fread(paste(c("~/FLAT/UKB/",trait,".summary"),collapse=""),data.table=F)
d[,2]=gsub("max","",d[,2])
d[,2]=gsub("min","",d[,2])

loci=fread(paste("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/UKB/int/",chr,sep=""),data.table=F)
lociinfo=loci[which(loci[,1]==chr),1:5]
fid=fread("~/FLAT/data/pheno.txt",data.table=F)
fid=fid[,1,drop=F]
d=d[which(d[,1] %in% lociinfo[,4]),]

get_locus_value=function(x){
  loci=gsub(" ","",x)
  
  print(loci)
  setwd(loci)
  system("sed 's/\t/ /g' sei | awk 'length($0)<300'   > sei1")
  sei=fread("sei1",data.table=F,fill=TRUE,header=F)
  unlink("sei1")
  if(nrow(sei)==0){VAR=NA}else{
    colnames(sei)=paste("V",1:ncol(sei),sep="")
    hapsgeno=fread("hapsgeno.gz",data.table=F,fill=TRUE)
    miss=502410-nrow(hapsgeno)
    if(miss>0){hapsgeno=rbind(data.frame(V1=rep(NA,miss),V2=rep(NA,miss)),hapsgeno)}  #fread ignore blank lines at the beginning, 
    VAR=colnames(sei)[which.max(as.numeric(sei[which(sei[,1]=="int")[1],2:40]))+1]
  }
  if(is.na(VAR)){ref=NA}else{
    sei=sei[,c("V1",VAR)]
    ref=sei[which(sei[,1]=="int")[1],2]
  }
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
    m=(m1+m2)/2   
    colnames(m)=loci
    rownames(m)=1:nrow(m)
  }
  setwd("..")
  return(m)
  gc()
}

pheno=fread(path,data.table=F)
if(ncol(pheno)<3){
  pheno=data.frame(FID=pheno[,1],IID=pheno[,1])
  for (i in 1:n){
    useid=sample(d[,1],343,replace=F,prob=d[,6])
    dat=mclapply(useid,get_locus_value,mc.cores=ncore)
    dat=dplyr::bind_cols(dat)    
    dat[]=apply(dat,2,scale)
    beta=rnorm(ncol(dat),mean=0,sd=1)
    dat=dat[match(pheno[,1],pheno[,1]),]
    dat[is.na(dat)]=0
    x=as.matrix(dat) %*% beta
    ssr <- sum((x-mean(x))^2) # sum of squared residuals
    e <- rnorm(length(x))
    e <- resid(lm(e ~ x))
    e <- e*sqrt((1-h2)/h2*ssr/(sum(e^2)))
    y <- x + e
    simupheno=paste("HFS",i,sep="")
    pheno[,simupheno]=y
  }
  fwrite(pheno,"~/gcta/simusnp/hfssimu",sep="\t")
} #simulate phenotype
pheno=pheno[,1:(n+2)]

susie_block=function(simutrait,block){
  locilist=lociinfo[which(lociinfo[,5]==block),4]
  dat=mclapply(locilist,get_locus_value,mc.cores=ncore)
  dat=dplyr::bind_cols(dat)
  rownames(dat)=fid[,1]
  l=pheno[,1]
  mat=dat[as.character(l),]
  #mat=as(apply(mat,2,imputecol),"matrix")
	fitted <- susie(data.matrix(mat),pheno[,simutrait],L = 10,verbose = FALSE,coverage=0.3)  #has Rfast really used multithreading? Not sure how to check.
  coef=susie_get_posterior_mean(fitted)
  pip=fitted$pip
  pip=data.frame(pip=pip,beta=coef,locus=names(pip),trait=simutrait)
  return(pip)
  gc()
}

for(simutrait in colnames(pheno)[3:(n+2)]){
  for (block in unique(lociinfo[,5])){
    if(file.exists(paste(c("snpsimu/",simutrait,"_",block),collapse=""))){next()}
    int=susie_block(simutrait=simutrait,block=block)
    fwrite(int,paste(c("snpsimu/",simutrait,"_",block),collapse=""),sep="\t")
  }
}






