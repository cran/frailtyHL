lognorm2.iter <-
function(X,B,Z,v,penalty,alpha,tun1=NULL,tun2=NULL,di,Mi,idx2,del){

n<-nrow(X)
p<-ncol(X)


iB<-B
iv<-v
ialpha<-alpha

rep<-1
err<-3
convergeind<-c(0)

while((rep<300)&(err>0.0000001)&(err<5)){

result<-lognorm2(X=X,B=iB,Z=Z,v=iv,penalty=penalty,alpha=ialpha,tun1=tun1,tun2=tun2,di,Mi,idx2,del)

err<-sum(abs(ialpha-result[[3]]))

#print("==== rep & err ===")
#print(c(rep,err))
#print(ialpha)

if (err<=0.0000001) {convergeind<-c(1)} else {convergeind<-c(-1)}

iB<-result[[1]]
iv<-result[[2]]
ialpha<-result[[3]]
mat<-result[[4]]
WL<-result[[5]]
Bi<-result[[6]]
hft<-result[[7]]
rep<-rep+1

#print(ialpha)
#print(convergeind)
#print(iB)
#print(iv)



}

#print(iB)

Zs<-Z
nrand<-length(Zs)
qs<-rep(0,nrand)

z<-NULL;v_h1<-NULL; alpha_h1<-NULL; rand.idx<-NULL

for (i in 1:nrand){
	z<-cbind(z,Zs[[i]])
	qs[i]<-dim(Zs[[i]])[2]
	v_h1<-c(v_h1,iv[[i]])
	alpha_h1<-c(alpha_h1,rep(ialpha[i],qs[i]))
	rand.idx<-c(rand.idx, rep(i, qs[i]))	
}

q<-sum(qs)
n.alpha<-length(ialpha)
## correlated model (rho)
## qs should be identical
if (n.alpha>nrand) {
  alpha.new<-ialpha[n.alpha]
  J<-matrix(1,nrand,nrand)
  Sig<-kronecker(J,diag(rep(alpha.new,qs[1])))
}

diag(Sig)<-alpha_h1
iSig<-solve(Sig)


beta_h1<-iB

x<-X

mat<-mat

u_h1=exp(v_h1)

mat11=t(x)%*%mat%*%x
mat12=t(x)%*%mat%*%z
mat13=t(z)%*%mat%*%z

#U=(1/alpha_h1)*diag(rep(1,q))
U=iSig
#print(det(mat13+U))

if ((err>5)|(det(mat13+U)<0.0000001)){GCV<-100000; se_beSp<-rep(0,len=p)}else{mmat=mat11-mat12%*%solve(mat13+U)%*%t(mat12)
hminv=solve(mmat+n*WL)
se_bem=sqrt(diag(hminv))
t=beta_h1/se_bem

#################################################

eta=as.numeric(x%*%beta_h1 + z%*%v_h1)
expeta=exp(eta)
Wi=diag(expeta)

minv=solve(-mmat-n*WL)
cov=mmat
iHs=(-minv)%*%cov%*%(-minv)
se_beS=sqrt(diag(iHs))
se_beSp=se_beS[1:p]
t1=beta_h1/se_beSp

eta=as.numeric(x%*%beta_h1 + z%*%v_h1)
expeta=exp(eta)
one=rep(1,len=n)
done=rep(1,len=idx2)
oq=rep(1,q)

pi=3.14159265359
term0=as.numeric(t(Mi)%*%expeta)

hlike1= (t(one)%*%(del*eta))-(t(done)%*%(di*log(term0)))
#hlike2=-(q/2)*log(2*pi)-(1/2)*sum(log(alpha_h1))
hlike2=-(q/2)*log(2*pi)-(1/2)*log(det(Sig))
#hlike3=-0.5*(t(v_h1)%*% diag(1/alpha_h1) %*%v_h1)
hlike3=-0.5*(t(v_h1)%*% iSig %*%v_h1)
hlike=hlike1+hlike2+hlike3

mat0=(Wi%*%Bi)
H22=t(z)%*%mat%*%z+U

oq=rep(1,len=q)
eiv=eigen(H22)$value
adj2=as.numeric(-0.5*(t(oq)%*%log(eiv/(2*pi))))
pvh=hlike+adj2

if(length(beta_h1)!=1){
WL=diag(hft/(abs(beta_h1)))
}else{
WL=hft/(abs(beta_h1))
}
iH=solve(mmat+n*WL)%*%(mmat)
elam=sum(diag(iH))
 
GCV=-2*pvh+log(n)*elam
}

if((penalty=="lasso")|(penalty=="scad")){
result<-list(round(iB,digits=5),round(se_beSp,digits=5),round(v_h1,digits=5),round(ialpha,digits=5),tun1,round(GCV,digits=5),rep,convergeind)
names(result)<-c("beta","se","v","alpha","tuning_parameter","BIC","repetition","convergence")
return(result)
}
if(penalty=="hl"){
result<-list(round(iB,digits=5),round(se_beSp,digits=5),round(v_h1,digits=5),round(ialpha,digits=5),c(tun1,tun2),round(GCV,digits=5),rep,convergeind)
names(result)<-c("beta","se","v","alpha","tuning_parameter","BIC","repetition","convergence")
return(result)
}

}
