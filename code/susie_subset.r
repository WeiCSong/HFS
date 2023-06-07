
###usage: Rscript susie_subset.r $trait $chr 0.1 or Rscript susie_subset.r $trait $chr nc
###0.1: extract phs for loci pip>0.1
###nc: extract phs for non-coding region


param <- commandArgs(trailingOnly=T)
trait= eval(paste(text=param[1]))
chr=eval(paste(text=param[2]))
listname=eval(paste(text=param[3]))


library(R.utils)
library(data.table)
library(doParallel)
d=fread("d",data.table=F)
d=d[grep("V",d[,3]),]
#cl <- makeCluster(4)
#registerDoParallel(cl)
load(paste(trait,"susieres.RData",sep="/"))
pheno=fread("~/FLAT/data/pheno.txt",data.table=F)
if(listname=="nc"){
  LIST=fread("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/data/gene/nc.bed",data.table=F)
  lociinfo=LIST[which(LIST[,1]==chr),]}else{
  listname=as.numeric(listname)
  loci=fread(paste("/dssg/home/acct-bmelgn/bmelgn-3/FLAT/UKB/int/",chr,sep=""),data.table=F)
  lociinfo=loci[,1:5]
  lociname=rownames(susieres)[which(susieres[,1]>listname)]
  lociinfo=lociinfo[which(lociinfo[,4] %in% lociname),]
}




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
    m[which(is.na(m[,1])),1]=ref
    rownames(m)=1:nrow(m)
  }
  setwd("..")
  coef=susieres[which(rownames(susieres)==loci),2]
  m=m*coef
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
  dat=Reduce("+", dat) 
  return(dat)
}

phs=rep(0,nrow(pheno))
for (block in unique(lociinfo[,5])){    #parallel computation implemented in susie_block function, so here use a for loop instead of apply function.
  int=susie_block(trait=trait,block=block)
  phs=phs+int
}
fwrite(data.frame(phs=phs),paste(c(trait,"/",listname,".txt.gz"),collapse=""),col.names=F)



