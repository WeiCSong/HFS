library(arrow)
library(data.table)
library(r2redux)
l=c("height","BMI","BMD","platelet","intelligence","SBP","chronotype","insomnia","menarche","FVC","asthma","rhinitis")
pheno=fread("~/FLAT/data/pheno.txt",data.table=F)

res=data.frame()
for(trait in l ){
pattern=paste(c("chr*/",trait,"/*",".feather"),collapse="")
file=system(paste(c("ls ",pattern),collapse=""),intern=T)
d=sapply(file,read_feather)
d=as.data.frame(d)
D=apply(d,2,function(x){x-mean(x)})
D=rowSums(D)
ldak=fread(paste(c("~/ldak/",trait,".sscore"),collapse=""),data.table=F)

d=data.frame(hfs=D,ldak=ldak[match(pheno[,1],ldak[,1]),5],group=pheno$group)

covar=colnames(pheno)[2:16]
cdat=pheno[,c(trait,covar)]
rownames(cdat)=pheno[,1]
covarmod=lm(paste(trait,"~.",sep=""),data=cdat)
if(length(unique(na.omit(cdat$trait)))==2){covarmod=glm(paste(trait,"~.",sep=""),data=cdat,family = "binomial")}
resid=resid(covarmod)
d$trait=resid[match(as.character(pheno[,1]),as.character(names(resid)))]	
val=c()
for(group in c("Ctest","African","Asian","South Asian")){
int=d[which(d$group==group),]
int=na.omit(int)
if(nrow(int)<5){next()}
model=lm(trait~ldak+hfs,data=int)
int$pred=predict(model,int)
int=int[,c(4,5,1,2)]
nv=nrow(int)
v1=c(1)
v2=c(2)
output=r2_diff(na.omit(int),v1,v2,nv)
output=unlist(output)
row=c(tool="hfs",group=group,output)
val=rbind(val,row)
v1=c(1)
v2=c(3)
output=r2_diff(na.omit(int),v1,v2,nv)
output=unlist(output)
row=c(tool="ldak",group=group,output)
val=rbind(val,row)
}
val=data.frame(val,trait=trait)
res=rbind(res,val)
}

fwrite(res,"hfs_ldak.r2",sep="\t")