

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

l=list.files(pattern="*summary")
l=gsub(".summary","",l)


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
#imputecol=function(x){
#  x[which(is.na(x))]=median(x,na.rm=T)
#  return(x)
#}


for(trait in l){
  if(file.exists(paste(trait,"tagging",sep="."))){next()}
  d=fread(paste(trait,"summary",sep="."),data.table=F)
  intd=d[which(d$PIPfull>0.95),1:2]
  intd=as.list(as.data.frame(t(intd)))
  dat=mclapply(intd,get_locus_value,mc.cores=20)
  dat=dplyr::bind_cols(dat)
  rownames(dat)=pheno[,1]
  
  covar=colnames(pheno)[2:16]
  cdat=pheno[,c(trait,covar)]
  rownames(cdat)=pheno[,1]
  covarmod=lm(paste(trait,"~.",sep=""),data=cdat)
  resid=resid(covarmod)
  resid=na.omit(resid)
  
  res=data.frame()
  for (ethic in unique(pheno$group)){
    ind=pheno[which(pheno$group==ethic),1]
    ind=intersect(ind,names(resid))
    data=data.frame(y=resid[ind],dat[ind,])
    if(length(which(!is.na(data$y)))==0){next()}
    model=lm(y~.,data=data)
    R2=summary(model)$adj.r.squared
    aic=AIC(model)
    n=ncol(dat)
    N=length(ind)
    res=rbind(res,c(ethic=ethic,R2=R2,AIC=aic,n=n,N=N))
  }
  res=data.frame(res)
  colnames(res)=c("ethic","R2","AIC","n","N")
  fwrite(res,paste(trait,"tagging",sep="."),sep="\t")
  gc()
}





