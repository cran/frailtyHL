lognorm<-function(X,B,Z,v,penalty,alpha,tun1=NULL,tun2=NULL,di,Mi,idx2,del, varfixed){

n<-nrow(X)
p<-ncol(X)

Zs<-Z; vs<-v
nrand<-length(Zs)

qs<-rep(0,nrand)
alphas<-alpha
n.alpha<-length(alphas)

Z<-NULL;v<-NULL; alpha<-NULL; rand.idx<-NULL

for (i in 1:nrand){
	qs[i]<-dim(Zs[[i]])[2]
	Z<-cbind(Z,Zs[[i]])
	v<-c(v,vs[[i]])
	alpha<-c(alpha,rep(alphas[i],qs[i]))
	rand.idx<-c(rand.idx, rep(i, qs[i]))	
}

q<-sum(qs)

## covariance matrix of v
Sig<-Sig.idx<-matrix(0,q,q)

## correlated model (rho)
## qs should be identical
if (n.alpha>nrand) {
  alpha.new<-alphas[n.alpha]
  J<-matrix(1,nrand,nrand)
  Sig<-kronecker(J,diag(rep(alpha.new,qs[1])))
  Sig.idx<-kronecker(J,diag(rep(n.alpha,qs[1])))
}

diag(Sig)<-alpha
iSig<-solve(Sig)
diag(Sig.idx)<-rand.idx




eta<-X%*%B+Z%*%v
expeta<-exp(eta)

cla0<-di/(t(Mi)%*%expeta)

Wi<-diag(as.numeric(expeta))
Ai<-diag(as.numeric(cla0))
done<-rep(1,len=idx2)

clam0<-as.numeric(Mi%*%Ai%*%done)

expclm<-expeta*clam0
wei<-diag(as.numeric(expclm))

one<-rep(1,len=n)
oq<-rep(1,len=q)
u_h0<-exp(v)

if (penalty=="hl"){
tau<-tun1
sig<-tun2
hft<-rep(0,len=p)
u<-rep(0,len=p)
delu<-rep(0,p)
for (i in 1:p){
u[i]<-(sqrt(8*tau*((B[i])^2)/sig+((2-tau)^2))+2-tau)/4
delu[i]<-(2*tau*abs(B[i])/sig)/(sqrt(8*tau*((abs(B[i]))^2)/sig+((2-tau)^2)))
hft[i]<-(((tau-2)/(2*tau)*delu[i]/(u[i]+0.00000001))+2*delu[i]/tau)/n
}
}

if (penalty=="lasso"){
lam<-tun1
hft=rep(0,len=p)
for (i in 1:p){
hft[i]<-lam
}
}

if (penalty=="scad"){
lam<-tun1
a=3.7
 hft=rep(0,len=p)
for (i in 1:p){
if ((0<=abs(B[i]))&(lam>abs(B[i]))){hft[i]=lam}
if ((lam<=abs(B[i]))&(a*lam>abs(B[i]))){hft[i]=(a*lam -abs(B[i]))/(a-1)}
}
}

if(length(as.numeric(hft/(abs(B)+0.00000001)))!=1){
WL<-diag(as.numeric(hft/(abs(B)+0.00000001)))
}else{
WL<-as.numeric(hft/(abs(B)+0.00000001))
}

dft1<-t(X)%*%(del-Wi%*%clam0)
dft1<-dft1-n*WL%*%B
#dft2<-t(Z)%*%(del-Wi%*%clam0)-(v/alpha)
dft2<-t(Z)%*%(del-Wi%*%clam0)-iSig%*%v

dft<-c(as.numeric(dft1),as.numeric(dft2))

Bi<-diag(as.numeric(clam0))
As<-diag(as.numeric(cla0^2/di))
mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%As%*%(t(Mi)%*%Wi)

#U<-diag(1/alpha)
U<-iSig

H<-rbind(cbind(t(X)%*%mat%*%X+n*WL,t(X)%*%mat%*%Z),cbind(t(Z)%*%mat%*%X,t(Z)%*%mat%*%Z+U))
Hinv<-solve(H)
H222<-Hinv[(p+1):(p+q),(p+1):(p+q)]

be_h0<-c(B,v)
be_h<-as.numeric(be_h0+(Hinv%*%dft))

B_h<-be_h[1:p]
v_h<-be_h[(p+1):(p+q)]

#expeta=as.numeric(exp(X%*%B_h+Z%*%v_h))



alpha_h<-alphas

if (varfixed==FALSE){
H22<-solve(t(Z)%*%mat%*%Z+U)

dk<-rep(0,n.alpha); 
ddk<-matrix(0,n.alpha,n.alpha)

for (i in 1:n.alpha){

dSig<-Sig*0
dSig[Sig.idx==i]<-1
  
#iB<-diag(1/alpha^2)
#diag(iB)[rand.idx!=i]<-0

iB<- iSig%*%dSig%*%iSig

c_vh<-iB%*%v_h
dv<-H22%*%c_vh

dexpeta<-as.numeric(expeta*as.numeric(Z%*%dv))
dcla0<--(as.numeric((di/((t(Mi)%*%expeta)^2)))*as.numeric(t(Mi)%*%dexpeta))

dWi<-diag(as.numeric(dexpeta))
dAi<-diag(as.numeric(dcla0))
dBi<-diag(as.numeric(Mi%*%dAi%*%done))
dvec<-2*(cla0*dcla0)
dAs<-diag(as.numeric(dvec/di))

dmat<-(dWi%*%Bi)+(Wi%*%dBi)-(dWi%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi)

dia1<- -iB
	
Hd = rbind(cbind(t(X)%*%dmat%*%X,t(X)%*%dmat%*%Z),cbind(t(Z)%*%dmat%*%X,t(Z)%*%dmat%*%Z+dia1))

#sgamma<--alphas[i]*sum(diag(Hinv%*%Hd))
#alpha_h[i]<-(as.numeric((t(v_h[rand.idx==i])%*%v_h[rand.idx==i])/(qs[i]-sgamma)))

tr1<-sum(diag(iSig%*%dSig))
#tr1<-tr(iSig, dSig)
tr2<-t(v_h)%*%iB%*%v_h
tr3<-sum(diag(Hinv%*%Hd))
#tr3<-tr(Hinv, Hd)
dk[i]<- -0.5*tr1+0.5*tr2-0.5*tr3
iSig.dSig<-iSig%*%dSig
H222.iB<-H222%*%iB

ddk[i,i]<- -0.5*sum(diag(iSig%*%dSig%*%iSig%*%dSig)) + t(v_h)%*%(iSig%*%dSig%*%iB)%*%v_h - 0.5*sum(diag(H222%*%iB%*%H222%*%iB)) +
      sum(diag(H222%*%iSig%*%dSig%*%iSig%*%dSig%*%iSig))
#ddk[i,i]<- -0.5*tr(iSig.dSig, iSig.dSig) + t(v_h)%*%(iSig.dSig%*%iB)%*%v_h - 0.5*tr(H222.iB, H222.iB) +tr(H222%*%iSig.dSig, iSig.dSig%*%iSig)

if (i>1){
  for (j in 1:(i-1)){
    dSigj<-Sig*0
    dSigj[Sig.idx==j]<-1
    
    iBj<- iSig%*%dSigj%*%iSig
    
    ddk[i,j]<-ddk[j,i]<- -0.5*sum(diag(iSig%*%dSig%*%iSig%*%dSigj)) + 0.5*t(v_h)%*%(iSig%*%dSigj%*%iB + iSig%*%dSig%*%iBj)%*%v_h +
                         -0.5*sum(diag(H222%*%iBj%*%H222%*%iB)) + sum(diag(H222%*%iSig%*%dSig%*%iSig%*%dSigj%*%iSig))
#    ddk[i,j]<-ddk[j,i]<- -0.5*tr(iB,dSigj) + 0.5*t(v_h)%*%(iSig%*%dSigj%*%iB + iSig%*%dSig%*%iBj)%*%v_h -0.5*tr(H222%*%iBj, H222%*%iB) + tr(H222%*%iB, dSigj%*%iSig)    
    
  }
}

}
## update of alpha 

alpha_h = alphas + solve(ddk, dk)

}

v_hs<-v_h
v_h<-NULL
for (i in 1:nrand) {v_h[[i]]<-v_hs[rand.idx==i]}




#alpha_h<-c(0.001,0.937)
return(list(B_h,v_h,alpha_h,mat,WL,Bi,hft))
}
