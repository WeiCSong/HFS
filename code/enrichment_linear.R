library(R.utils)
library(data.table)
library(Rfast)
library(speedglm)
library(mgcv)

baseline=fread("baseline.txt.gz",data.table=F)
#load("allpath.RData")
baseline[which(is.na(baseline$range)),"range"]=0
loci=fread("../loci.bed",data.table=F)
#term=colnames(allpath)[1:3874]
#baseline$all=allpath$all
#baseline$block=loci[match(baseline[,1],loci[,4]),5]


#
#for(trait in c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis","alcohol","smoke")){
#  res=data.frame()
#  summary=fread(paste(c("~/FLAT/UKB/",trait,".summary"),collapse=""),data.table=F)
#  y=summary[match(baseline[,1],summary[,1]),"PIPfull"]
#  y[is.na(y)]=0
#  for (T in term){
#    dat=cbind(baseline[,-1],allpath[,T])
#    colnames(dat)[ncol(dat)]=T
#	  colnames(dat)=gsub("_/_","_",colnames(dat))
#	  colnames(dat)=gsub("-","_",colnames(dat))
#	  colnames(dat)=gsub("+","p",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub(",","p",fixed=TRUE,colnames(dat))
#	  l=setdiff(colnames(dat),c("block","pip"))
#	  l=paste(l,collapse=" + ")
#	  form=paste(c("pip~",l),collapse="")
#    dat$pip=y
#    model=speedglm(as.formula(form),data=dat)
#		sum=summary(model)$coefficients[T,]
#    res=rbind(res,sum)
#  }
#  fwrite(res,paste(trait,"separate.pathway",sep="_"),row.names=TRUE,sep="\t")
#}


#load("atac.RData")
#atac=atac[match(baseline[,1],atac[,1]),-1]
#cell=colnames(atac)
#
#for(trait in c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis","alcohol","smoke")){
#  res=data.frame()
#  summary=fread(paste(c("~/FLAT/UKB/",trait,".summary"),collapse=""),data.table=F)
#  y=summary[match(baseline[,1],summary[,1]),"PIPfull"]
#  y[is.na(y)]=0
#  for (C in cell){
#    dat=cbind(baseline[,-1],atac[,C])
#    colnames(dat)[ncol(dat)]=C    
#	  colnames(dat)=gsub("(","",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub("_/_","_",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub(")","",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub("-","_",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub("+","p",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub(",","p",fixed=TRUE,colnames(dat))
#	  l=setdiff(colnames(dat),c("block","pip"))
#	  l=paste(l,collapse=" + ")
#	  form=paste(c("pip~",l),collapse="")
#    dat$pip=y
#    model=speedglm(as.formula(form),data=dat)
#		sum=summary(model)$coefficients[C,]
#    res=rbind(res,sum)
#  }
#  fwrite(res,paste(trait,"separate.atac",sep="_"),row.names=TRUE,sep="\t")
#}
#
#load("chromatin.RData")
#chromatin=chromatin[match(baseline[,1],chromatin[,1]),-1]
#cell=colnames(chromatin)
#
#for(trait in c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis","alcohol","smoke")){
#  res=data.frame()
#  summary=fread(paste(c("~/FLAT/UKB/",trait,".summary"),collapse=""),data.table=F)
#  y=summary[match(baseline[,1],summary[,1]),"PIPfull"]
#  y[is.na(y)]=0
#  for (C in cell){
#    dat=cbind(baseline[,-1],chromatin[,C])
#    colnames(dat)[ncol(dat)]=C    
#	  colnames(dat)=gsub("(","",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub("_/_","_",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub(")","",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub("-","_",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub("+","p",fixed=TRUE,colnames(dat))
#	  colnames(dat)=gsub(",","p",fixed=TRUE,colnames(dat))
#	  l=setdiff(colnames(dat),c("block","pip"))
#	  l=paste(l,collapse=" + ")
#	  form=paste(c("pip~",l),collapse="")
#    dat$pip=y
#    model=speedglm(as.formula(form),data=dat)
#		sum=summary(model)$coefficients[C,]
#    res=rbind(res,sum)
#  }
#  fwrite(res,paste(trait,"separate.chromatin",sep="_"),row.names=TRUE,sep="\t")
#}

load(".RData")
cell=colnames(allcell)
l=rowSums(allcell)
baseline$allcell=ifelse(l>0,1,0)

for(trait in c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis","alcohol","smoke")){
  res=data.frame()
  summary=fread(paste(c("~/FLAT/UKB/",trait,".summary"),collapse=""),data.table=F)
  y=summary[match(baseline[,1],summary[,1]),"PIPfull"]
  y[is.na(y)]=0
  for (C in cell){
    dat=cbind(baseline[,-1],allcell[,C])
    colnames(dat)[ncol(dat)]=C    
	  colnames(dat)=gsub("(","",fixed=TRUE,colnames(dat))
	  colnames(dat)=gsub("_/_","_",fixed=TRUE,colnames(dat))
	  colnames(dat)=gsub(")","",fixed=TRUE,colnames(dat))
	  colnames(dat)=gsub("-","_",fixed=TRUE,colnames(dat))
	  colnames(dat)=gsub("+","p",fixed=TRUE,colnames(dat))
	  colnames(dat)=gsub(",","p",fixed=TRUE,colnames(dat))
	  l=setdiff(colnames(dat),c("block","pip"))
	  l=paste(l,collapse=" + ")
	  form=paste(c("pip~",l),collapse="")
    dat$pip=y
    model=speedglm(as.formula(form),data=dat)
		sum=summary(model)$coefficients[C,]
    res=rbind(res,sum)
  }
  fwrite(res,paste(trait,"separate.cell",sep="_"),row.names=TRUE,sep="\t")
}







