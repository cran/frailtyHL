lognorm1.BIC <- function(X,B,Z,v,penalty,alpha,tun1,tun2,di,Mi,idx2,del,varfixed){

if((penalty=="lasso")|(penalty=="scad")){
bicrep<-length(tun1)
bic<-c()
for (g in 1:bicrep){
result_bic<-lognorm1.iter(X,B,Z,v,penalty,alpha,tun1=tun1[g],tun2,di,Mi,idx2,del,varfixed)
bic[g]<-result_bic$BIC
}
sel_tun1<-tun1[which.min(bic)]
result<-lognorm1.iter(X,B,Z,v,penalty,alpha,tun1=sel_tun1,tun2,di,Mi,idx2,del,varfixed)
}

if(penalty=="hl"){
bicrep<-length(tun1)*length(tun2)
tunset1<-rep(tun1,each=length(tun2))
tunset2<-rep(tun2,length(tun1))

bic<-c()
for (g in 1:bicrep){
result_bic<-lognorm1.iter(X,B,Z,v,penalty,alpha,tun1=tunset1[g],tunset2[g],di,Mi,idx2,del,varfixed)
bic[g]<-result_bic$BIC
}
sel_tun1<-tunset1[which.min(bic)]
sel_tun2<-tunset2[which.min(bic)]
result<-lognorm1.iter(X,B,Z,v,penalty,alpha,tun1=sel_tun1,tun2=sel_tun2,di,Mi,idx2,del,varfixed)
}

return(result)
}
