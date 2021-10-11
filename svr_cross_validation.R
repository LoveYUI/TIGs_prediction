#regression in ftp
library("e1071")
load("./ml_train_r.Rdata")
load("./ml_test_ope_r.Rdata")
conclusion<-function(model,test)
{
  pre_svm<- predict(model,newdata =  test)
  rss<-(pre_svm-test[,ncol(test)])%*%(pre_svm-test[,ncol(test)])
  mean_test<-mean(test[,ncol(test)])
  tss<-(mean_test-test[,ncol(test)])%*%(mean_test-test[,ncol(test)])
  print(1-rss/tss)
}
gridsearch<-function(gm,c,cr,df)
{ 
  best_mse<-100
  for(i in gm)
  {
    for(j in c)
    {
      print(paste("start"," ",as.character(i)," ",as.character(j)))
      model<-svm(target~.,data=df,type="nu-regression",probability=TRUE,cross=cr,gamma=i,cost=j)
      print(mean(as.numeric(model$MSE)))
      
      conclusion(model,ml_test_ope_r)
      if(mean(model$MSE)<best_mse)
      {
        best_model<-model
        best_mse<-mean(model$MSE)
      }
    }
  }
  return(best_model)
}
model_svm_r<-gridsearch(2^-10.5,c(3,3.5,4,4.5,5),5,ml_train_r)
save(model_svm_r,file="./model_svm_r.Rdata",version = 2)
