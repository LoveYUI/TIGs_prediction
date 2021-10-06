library("e1071")

#svm 5-cross-validation 
load("./ml_train_c.Rdata")
load("./ml_test_ope.Rdata")
conclusion<-function(model,test)
{
  pre_svm<- predict(model,newdata =  test,probability = TRUE)
  tn<-0
  tp<-0
  fn<-0
  fp<-0
  for(i in 1:nrow(test))
  {
    if(test$target[i]=="1")
    {
      if(pre_svm[i]=="1"){tp<-tp+1}else{fp<-fp+1}
    }else{
      if(pre_svm[i]=="1"){fn<-fn+1}else{tn<-tn+1}
    }
  }
  conc_svm<-cbind(c("accuracy","recall","precision","F_measure"),c((tp+tn)/(tp+tn+fp+fn),tp/(tp+fn),tp/(tp+fp),2*(tp/(tp+fp))*(tp/(tp+fn))/(tp/(tp+fp)+tp/(tp+fn))))
  print(conc_svm)
}
gridsearch<-function(gm,c,cr,df)
{ 
  best_acc<-0
  for(i in gm)
  {
    for(j in c)
    {
      print(paste("start"," ",as.character(i)," ",as.character(j)))
      model<-svm(target~.,data=df,type="C-classification",probability=TRUE,cross=cr,gamma=i,cost=j)
      print(model$accuracies)
      if(mean(model$accuracies)>best_acc)
      {
        best_model<-model
        best_acc<-mean(model$accuracies)
      }
    }
  }
  return(best_model)
}
model_svm<-gridsearch(c(2^-9,2^-10,2^-10.5,2^-11,2^-12),c(1,1.25,1.5,1.75,2,2.25),5,ml_train_c)
conclusion(model_svm,ml_test_ope)
summary(model_svm)
save(model_svm,file="./model_svm.Rdata",version = 2)
