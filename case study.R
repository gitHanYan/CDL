library(ggplot2)
library(ggpubr)
source("function.r")
####import datasets
patient=read.csv("LUAD_mRNA_patient.csv",row.names=1)
normal=read.csv("LUAD_mRNA_normal.csv",row.names=1)
cov_P=cov(patient)
cov_N=cov(normal)
ob=nrow(patient)+nrow(normal)
dim=72
#tuning parameters by BIC
LUAD=list(mRNA,miRNA,methy)
AIC_list=list()
aic_f_norm=vector()
aic_i_norm=vector()
####tuning parameter selection
for (i in 1:length(lambda)) {
  delta=CDL(patient,normal,K,lambda[i],rho)
  aic_f[i]=ob*norm((0.5*(cov_P%*%delta%*%cov_N+cov_N%*%delta%*%cov_P)-cov_P+cov_N),"F")+
    2*sum(sign(abs(delta)))
  aic_i[i]=ob*norm((0.5*(cov_P%*%delta%*%cov_N+cov_N%*%delta%*%cov_P)-cov_P+cov_N),"I")+
    2*sum(sign(abs(delta)))
}

after_tuning=CDL(patient,normal,K,lambda[i],rho)

pic=sign(abs(after_tuning))

###picture drawing
l1=which(cut_AB==1)
l2=which(cut_AB==2)
l3=which(cut_AB==3)
vcolor=c("#7bbfea","#ed1941","#ffd400")

jpeg(file="c1.jpeg",width=3900,height=3900,res=650)
data_stru_1=graph_from_adjacency_matrix(pic,mode = "undirected")
V(data_stru_1)[l1]$color=vcolor[1]
V(data_stru_1)[l2]$color=vcolor[2]
V(data_stru_1)[l3]$color=vcolor[3]
#set.seed(1)
plot(data_stru_1,vertex.size=5,vertex.label.cex=0.7,
     vertex.label.dist=1.2,edge.width=1.5)
dev.off()
