library(data.table)
l=c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis")
pheno=fread("~/FLAT/data/pheno.txt",data.table=F)
loci=fread("~/FLAT/data/loci.bed",data.table=F)
classvec=paste("V",2:40,sep="")
hm3=fread("/dssg/home/acct-bmelgn/bmelgn-3/ldsc/w_hm3.snplist",data.table=F)
hm3=hm3[,1]
d=fread("newsummary",data.table=F)
d=split(d,f=d[,1])


summ=c()
for(trait in l){
	susie=d[[trait]]
	l=match(susie[,2],loci[,4])
	susie$chr=loci[l,1]
	susie$mid=(loci[l,2]+loci[l,3])/2
	type=ifelse(length(unique(na.omit(pheno[,trait])))==2,"binary","continuous")
	int=pheno[which(pheno$group=="Caucasian" & !is.na(pheno[,trait])),]
	N=nrow(int)
	prop=NA
	if(type=="binary"){prop=length(which(int[,trait]==1))/N}
	int=susie[which(susie[,6]<5e-8),]
	Nsig=nrow(int)
	allclass=int[,2]
#	dom=length(grep("m",allclass))
#	dom=dom/Nsig
#	allclass=gsub("min","",allclass)
#	allclass=gsub("max","",allclass)
#	allclass=table(allclass)
#	allclass=allclass[classvec]
	int=susie[which(susie[,6]<5e-8),]
	int=int[order(int[,6]),]
	int$mid=as.numeric(int$mid)
	k=nrow(int)
	for(i in 1:k){
		if(nrow(int)==i-1){break}
		chr=int[i,"chr"]
		mid=int[i,"mid"]
		l=which(int$chr==chr & abs(int$mid-mid)<100000)
		if(length(l)==1){next()}
		l=setdiff(l,1:i)
		int=int[-l,]
		
	}
	Nind=nrow(int)
	GWAS=fread(paste(c("~/regenie/",trait,".regenie"),collapse=""),data.table=F)
	GWAS=GWAS[which(GWAS[,3] %in% hm3),]
	t_statistic <- qt(1 - GWAS$LOG10P, N, lower.tail = F)
	lambdaSNP=median(na.omit(t_statistic)^2)
	t_statistic <- qt(1 - susie[,6], N, lower.tail = F)
	lambdaHFS=median(na.omit(t_statistic)^2)
	row=c(trait=trait,type=type,N=N,prop=prop,Nsig=Nsig,Nind=Nind,lambdaSNP=lambdaSNP,lambdaHFS=lambdaHFS)
	summ=rbind(summ,row)
}
summ=data.frame(summ)
fwrite(summ,"newtraitstat.csv")