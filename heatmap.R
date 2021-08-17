#heatmap
library("ggplot2")
library("ggpubr")
install.packages("pheatmap")
library("pheatmap")
#######################
ttt<-cbind(c(1,2,3),c(5,6,7),c(1,5,9))
pheatmap(ttt)
apm_cancertype_tmb_heat<-matrix(nrow=32,ncol=18)
for(i in 31:33)
{
  for(j in 1:18)
  {
    a<-as.numeric(OS_TMB_apmgene_all[(which(OS_TMB_apmgene_all[,1]==cancer_type[i]&OS_TMB_apmgene_all[,6]=="1")),j+6])
    b<-as.numeric(OS_TMB_apmgene_all[(which(OS_TMB_apmgene_all[,1]==cancer_type[i]&OS_TMB_apmgene_all[,6]=="0")),j+6])
    if(mean(a)>mean(b))
    {
      apm_cancertype_tmb_heat[i-1,j]<--log(t.test(a,b)$p.value,2)
    }else
    {
      apm_cancertype_tmb_heat[i-1,j]<-log(t.test(a,b)$p.value,2)
    }
  }
}
bk <- c(seq(-20,-0.1,by=0.01),seq(0,20,by=0.01))
pheatmap(t(apm_cancertype_tmb_heat),breaks=seq(-5,5,1),color = colorRampPalette(c("blue", "white", "red"))(10),cluster_rows=FALSE,fontsize_col=12)
?pheatmap
colnames(apm_cancertype_tmb_heat)<-apmgene_esemble[,1]
rownames(apm_cancertype_tmb_heat)<-cancer_type[c(1:29,31,32,33)]
#################################
cal_deathandalive<-function(ma,day,tmb_high)                                          #计算os，如计算5年生存率，day=365x5，用(followup>500+death>500)的除以(death>500+followup>500+death<500)
{                                                                               #因为有很多na老报错，所以我这里写的特别麻烦，你有时间可以改一下                                                                                                                                                                 
  a=0                                                                           #参数tmb_high，如果是1，就只计算tmb为1的，如果是0就只计算tmb是0的
  b=0 
  if(tmb_high==2)
  {
    for(i in 1:nrow(ma))
    {
      if(!is.na(ma[i,5]))
      {
        if(ma[i,2]=='')
        {
          if(ma[i,3]!='')
          {
            if(as.integer(ma[i,3])>day)
              a<-a+1
          }
        } else
        {
          if(as.integer(ma[i,2])>day)
            a<-a+1 else
              b<-b+1
        }
      }
    }
  }else{                                                                                          #这里ma是一个5列的数据，第一列编号，第二列death，第三列followup，第四列tmb，第五列是tmb高或低（0or1）
    for(i in 1:nrow(ma))
    {
      if(!is.na(ma[i,5]))
      {
        if(ma[i,5]==tmb_high)
        {
          if(ma[i,2]=='')
          {
            if(ma[i,3]!='')
            {
              if(as.integer(ma[i,3])>day)
                a<-a+1
            }
          } else
          {
            if(as.integer(ma[i,2])>day)
              a<-a+1 else
                b<-b+1
          }
        }
      }
    }
  }
  ret<-c(a,b)
  return(ret)
}
tigs_hr<-function(os_tmb,day)
{
  sign<-"-"
  daa0<-cal_deathandalive(os_tmb,day,0)
  daa1<-cal_deathandalive(os_tmb,day,1)
  hr<-(daa1[2]/daa1[1])/(daa0[2]/daa0[1])
  var<-sqrt((1/daa1[1])+(1/daa1[2])+(1/daa0[1])+(1/daa0[2]))
  ci_95<-paste("(",as.character(round(exp(log(hr)-1.96*var),digits=3)),",",as.character(round(exp(log(hr)+1.96*var),digits=3)),")")
  p<-pnorm(0,log(hr),var)
  if(p>0.5)
  {
    p<-1-p
  }
  p<-format(signif(p,3),scientific = TRUE)
  if(as.numeric(p)<0.05)
  {sign<-"*"}
  if(as.numeric(p)<0.01)
  {sign<-"**"}
  hr<-round(hr,3)
  ret<-c(hr,daa0[1]+daa0[2]+daa1[1]+daa1[2],ci_95,p,sign)
  return(ret)
}
for(i in 31:33)
{
  tigs<-OS_TMB_apmgene_all[which(OS_TMB_apmgene_all[,1]==cancer_type[i]),c(1:5,27)]
  tigs_med<-median(as.numeric(tigs[,6]))
  for(i in 1:nrow(tigs))
  {
    if(as.numeric(tigs[i,6])>tigs_med)
    {
      tigs[i,6]="1"
    }else
    {
      tigs[i,6]="0"
    }
  }
  print(tigs[1,1])
  print(tigs_hr(tigs[,c(2:6)],3*365))
}
##########################
pass<-0
for(i in 1:2578)
{
  if(shapiro.test(as.numeric(svm_all_c[c(1:5000),100]))$p.value>0.05)
  {
    pass<-pass+1
  }
  
  
}

cor_TMB_APM<-na.omit(OS_TMB_apmgene_all[,c(4)])
for(i in 1:18)
{
  print(cor(as.numeric(OS_TMB_apmgene_all[,5]),as.numeric(OS_TMB_apmgene_all[,6+i])))
}
print(apmgene_esemble[,1])

#y=apm,x=type,box
median(OS_TMB_geneall_all[,])
apm_type_box<-data.frame(OS_TMB_apmgene_all[,c(1,6,25)],stringsAsFactors = FALSE)
apm_type_box[,3]<-as.numeric(apm_type_box[,3])
colnames(apm_type_box)<-c("Type","TMB","APMscore")
apm_type_box[,1]<-toupper(apm_type_box[,1])
apm_type_box[which(apm_type_box[,2]==0),2]<-"low"
apm_type_box[which(apm_type_box[,2]==1),2]<-"high"
ggboxplot(apm_type_box, x = "Type", y = "APMscore",color="TMB",palette=c("skyblue","#ff0000"))
?ggboxplot

gene_tigs<-OS_TMB_geneall_all[,c(6:2584)]
save(gene_tigs,file="C:/Users/11569/Desktop/毕设/gene-tigs.RData")

