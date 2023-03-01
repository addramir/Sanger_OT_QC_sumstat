setwd("~/Projects/Sanger_OT_QC_sumstat/02_qc_draft/")
library("data.table")
x=fread("QC_results.csv",data.table=F)
x=x[,-7]

for (i in (2:ncol(x))){
  hist(x[,i],n=30,main=colnames(x)[i])
}

library(corrplot)
corrplot(cor(x[,-1]),method = "square",hclust.method = "ward.D2")


for (i in (2:ncol(x))){
  l=x[,i]
  
  iqr=IQR(l)
  q25=quantile(l,probs = 0.25)
  q75=quantile(l,probs = 0.75)
  qmin=q25-3*iqr
  qmax=q75+3*iqr
  ind=which(l<qmin | l>qmax)
  print(paste(colnames(x)[i],length(ind)))
}

L=NULL
#####
i=3
print(colnames(x)[i])
ind=x[,i]>=5
table(ind)
L=c(L,x[ind,1])
#####
i=4
print(colnames(x)[i])
ind=x[,i]<100
table(ind)
L=c(L,x[ind,1])
#####
i=6
print(colnames(x)[i])
ind=x[,i]<(-0.05) | x[,i]>(0.05)
table(ind)
L=c(L,x[ind,1])
#####
i=7
print(colnames(x)[i])
ind=x[,i]>(0.3)
table(ind)
L=c(L,x[ind,1])
#####
i=9
print(colnames(x)[i])
ind=x[,i]>(1.05) | x[,i]<(0.95)
table(ind)
L=c(L,x[ind,1])
#####
i=11
print(colnames(x)[i])
ind=x[,i]>(0.05) | x[,i]<(-0.05)
table(ind)
L=c(L,x[ind,1])
######
length(unique(L))

write.table(x=unique(L),file="filtered_study_ids.txt",col.names=F,row.names=F,quote=F,sep="\t")

table(grepl("GCST",x[,1]))
table(grepl("FINNGEN",x[,1]))
table(grepl("NEALE",x[,1]))
table(grepl("SAIGE",x[,1]))


i=9
print(colnames(x)[i])
ind=(abs(x[,i]-1)<0.05 & abs(x[,i]-1)>0.01)
table(ind)
l=x[ind,1]
L=unique(L)
l=l[!l%in%L]
write.table(x=l,file="pvalues_from_b_se_study_ids.txt",col.names=F,row.names=F,quote=F,sep="\t")
