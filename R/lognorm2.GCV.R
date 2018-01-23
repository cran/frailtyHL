lognorm2.GCV <-
function(X,B,Z,v,penalty,alpha,tun1,tun2,di,Mi,idx2,del){

if((penalty=="lasso")|(penalty=="scad")){
gcvrep<-length(tun1)
gcv<-c()
for (g in 1:gcvrep){
result_gcv<-lognorm2.iter(X,B,Z,v,penalty,alpha,tun1=tun1[g],tun2,di,Mi,idx2,del)
gcv[g]<-result_gcv$BIC
}
sel_tun1<-tun1[which.min(gcv)]
result<-lognorm2.iter(X,B,Z,v,penalty,alpha,tun1=sel_tun1,tun2,di,Mi,idx2,del)
}

if(penalty=="hl"){
gcvrep<-length(tun1)*length(tun2)
tunset1<-rep(tun1,each=length(tun2))
tunset2<-rep(tun2,length(tun1))

gcv<-c()
for (g in 1:gcvrep){
result_gcv<-lognorm2.iter(X,B,Z,v,penalty,alpha,tun1=tunset1[g],tunset2[g],di,Mi,idx2,del)
gcv[g]<-result_gcv$BIC
}
sel_tun1<-tunset1[which.min(gcv)]
sel_tun2<-tunset2[which.min(gcv)]
result<-lognorm2.iter(X,B,Z,v,penalty,alpha,tun1=sel_tun1,tun2=sel_tun2,di,Mi,idx2,del)
}

return(result)
}
