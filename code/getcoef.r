library(data.table)
library(R.utils)

trait=c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis","alcohol","smoke")
for(chr in 1:22){
  coef=c()
  for (t in trait){
    polyfun=fread(paste(c(t,"/polyfun_",t,".",chr,".agg.txt.gz"),collapse=""),data.table=F)
    polyfun=polyfun[,c(3,4,12)]
    colnames(polyfun)=c("SNP","A1","coef")
    prscs=fread(paste(c("~/PRScs-master/",t,"/",t,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"),collapse=""),data.table=F)
    prscs=prscs[,c(2,4,6)]
    colnames(prscs)=c("SNP","A1","coef")
    sbrc=fread(paste(c("~/sbayesrc/",t,"/chr",chr,".sbrc"),collapse=""),data.table=F)
    sbrc=sbrc[,1:3]
    colnames(sbrc)=c("SNP","A1","coef")
    d=rbind(polyfun,prscs)
    d=rbind(d,sbrc)
    d=d[!duplicated(d[,1]),1:2]
    polyfun=polyfun[match(d[,1],polyfun[,1]),]
    polyfun[is.na(polyfun)]=0
    polyfun[,3]=ifelse(polyfun[,2]==d[,2],polyfun[,3],-polyfun[,3])
    prscs=prscs[match(d[,1],prscs[,1]),]
    prscs[is.na(prscs)]=0
    prscs[,3]=ifelse(prscs[,2]==d[,2],prscs[,3],-prscs[,3])    
    sbrc=sbrc[match(d[,1],sbrc[,1]),]
    sbrc[is.na(sbrc)]=0
    sbrc[,3]=ifelse(sbrc[,2]==d[,2],sbrc[,3],-sbrc[,3])      
    d=data.frame(d,polyfun[,3],prscs[,3],sbrc[,3])
    colnames(d)[3:5]=paste(t,c("polyfun","prscs","sbrc"),sep="_")
    if(length(coef)==0){coef=d}else{coef=merge(coef,d,by=c("SNP","A1"))}
  }
  fwrite(coef,paste(chr,".cpef",sep=""),sep="\t")
}    
    