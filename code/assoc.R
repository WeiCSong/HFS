param <- commandArgs(trailingOnly=T)
chunk= eval(paste(text=param[1]))

phenofile="~/FLAT/data/pheno.txt"
is.binary=function(x){length(unique(na.omit(x))) ==2}
library(R.utils)
library(data.table)
#library(xgboost)
#library(rhdf5)
#library(arrow)
library(susieR)

chunk=fread(chunk,data.table=F)
chunk=chunk[,4]

phenodata=fread(phenofile,data.table=F)
#phenodata$sex=ifelse(phenodata$sex=="Male",1,0)
covar=colnames(phenodata)[2:16]
PHENO=colnames(phenodata)[17:30]
train=which(phenodata$group=="Caucasian" & phenodata$relate==0) #all individual files will have the same order as phenofile, so numeric index is used, not character individual ids.


for(locus in chunk){
  print(locus)
  setwd(as.character(locus))
  if(!file.exists("sei")){
    setwd("..")
    next
  }
 
  system("sed 's/\t/ /g' sei > sei1")
  sei=fread("sei1",data.table=F,header=F)
  unlink("sei1")
  hapsgeno=fread("hapsgeno.gz",data.table=F,fill=TRUE)
  miss=502410-nrow(hapsgeno)
  if(miss>0){hapsgeno=rbind(data.frame(V1=rep(NA,miss),V2=rep(NA,miss)),hapsgeno)}  #fread ignore blank lines at the beginning, 
  sei[which(sei[,1]=="int"),1]=1
  sei[,1]=as.numeric(sei[,1])
  sei=sei[order(sei[,1]),]
  sei=sei[match(1:sei[nrow(sei),1],sei[,1]),-1] #just to avoid missing haplotype ids.
  rownames(sei)=1:nrow(sei)
  m1=sei[hapsgeno[,1],]
  m2=sei[hapsgeno[,2],]
  m=(m1+m2)/2
#  SD=apply(m[train,],2,function(x){sd(na.omit(x))})
#  SD=SD[order(SD,decreasing=T)]
#  SD=names(SD)[1:11]
  ref=colnames(sei)[which.max(sei[1,])]
  m=m[,ref,drop=F]
#  SD=unique(c(SD,ref))
#  m=m[,SD]
#  m1=m1[,SD]
#  m2=m2[,SD]
#  m=cbind(m,Reduce(pmax, list(m1,m2)))
#  m=cbind(m,Reduce(pmin, list(m1,m2))) #m contains all sample
#  rownames(m)=rownames(phenodata)
#  colnames(m)=c(SD,paste("max",SD,sep=""),paste("min",SD,sep=""))
  rm(m1)
  rm(m2)
  rownames(m)=1:nrow(m)
  for(pheno in PHENO){
  	pd=phenodata[train,c("relate",pheno,covar)]
  	pd=na.omit(pd)
    pd=pd[which(pd$relate==0),-1]
    colnames(pd)[1]="trait"
#    cmod=lm(trait~.,data=pd)
#    resid=resid(cmod)
#    
#  	#SUSIE on m[rownames(pd),]
#    mat=as(na.omit(m[names(resid),]),"matrix")
#  	fitted <- susie(mat,resid[rownames(mat)],L = 10,verbose = FALSE,coverage=0.3)
#    #if(is.null(maxprob)){next}
#  	#sets <- susie_get_cs(fitted,X = m[rownames(pd),],coverage = 0.5,min_abs_corr = 0.1)
#  	index=names(fitted$pip)[which.max(fitted$pip)]
#    if(max(fitted$pip)<0.5){index=ref}
#    pip=fitted$pip[index]
#    names(pip)=c()
    val=m[rownames(pd),1]
    int=data.frame(val=val,pd)
    form=paste(covar,collapse="+")
    form=paste(c("trait~val+",form),collapse="")
    bin=is.binary(phenodata[,pheno])
    if(bin){model=glm(form,data=int,family = binomial)}else{model=lm(form,data=int)}
    d=summary(model)$coefficients
    coef=d[match("val",rownames(d)),1]
    pval=d[match("val",rownames(d)),4]
    res=data.frame(trait=pheno,id=locus,var=ref,PIPin=NA,h=coef,p=pval)
    fwrite(res,"~/FLAT/UKB/newsummary",sep="\t",col.names=F,append=T)
  }
  setwd("..")

}















