mlmfit1 <-
function(jm1,data_surv,Maxiter=200) {
 mc <- match.call()


 Vfunc <- function(x) {
   return(dnorm(x)/(1-pnorm(x)))
}

xifunc <- function(x) {
   return(Vfunc(x)*(Vfunc(x)-x))
}


 res1<-design_frailtyHL(jm1[[3]],data=data_surv)


yy<-log(res1[[1]])

### dispersion parameters #
phi<-1
alpha<-0.05

## x matrix
xx<-res1[[2]]

## z matrix
zz<-res1[[4]][[1]]

nn<-nrow(yy)
pp<-ncol(xx)+1
qq<-ncol(zz)

## delta

delta<-res1[[3]]

## 
beta_h0<-rbind(0,res1[[9]])
#beta_h0<-rbind(4.9,-0.2)

##### initial values by using linear mixed model
#res_lm<-lmer(yy~gp+(1|patient),skin_pair1)
#res_lm2<-summary(res_lm)
#beta_h0<-matrix(res_lm2[[10]][,1],pp,1)
#phi<-res_lm2[[11]]^2

x<-cbind(rep(1,nn),xx)
xx<-x
z<-zz
v_h0<-res1[[11]]
delta<-res1[[3]]
y<-yy

convergence<-1
i<-0

iteration<-1
while (convergence>0.000001 && iteration<=Maxiter) {
iteration<-iteration+1

## weight matrix
muij<-x %*% beta_h0 + z %*% v_h0
mij<-1/sqrt(phi)*(y-muij)
weight1<-(delta+(1-delta)*xifunc(mij))
ystar<-y*delta+(muij+sqrt(phi)*Vfunc(mij))*(1-delta)
adj_yy<-muij+1/weight1*(ystar-muij)
adj_zz2<-matrix(rep(0,qq),qq,1)
adj_zz<-rbind(adj_yy,adj_zz2)

oo<-matrix(0,qq,pp)
ii<-diag(rep(1,qq))
TT<-rbind(cbind(xx,zz),cbind(oo,ii))

weight11<-diag(as.vector(weight1/phi))
weight2<-diag(rep(1/alpha,qq))
weight<-dbind(weight11,weight2)
## if (i==0) weight<-diag(rep(1,nrow(TT)))

hinv<-solve(t(TT)%*%weight%*%TT)
old_beta<-beta_h0
beta_vv<-hinv%*%t(TT)%*%weight%*%adj_zz
temp1<-1
temp2<-pp
beta_h0<-matrix(beta_vv[temp1:temp2,1],pp,1)
temp1<-1+pp
temp2<-pp+qq
v_h0<-matrix(beta_vv[temp1:temp2,1],qq,1)

se_beta_vv<-sqrt(diag(hinv))

#################### alpha_h #############
old_alpha<-alpha

We<-diag(as.vector(weight1))
lam<-as.vector(phi/alpha)
#lam=1
i_q<-diag(rep(1,qq))
U<-lam*i_q



iW<- solve(t(z)%*%We%*%z +U)  

iD<- -diag(rep(1/alpha^2,qq))
#iD<- -(1/alpha^2)*diag(rep(1,qq))
dva<- -phi*iW%*%iD%*%v_h0

da_weight1<- -1/sqrt(phi)*z %*% dva
da_weight2 <- (2*Vfunc(mij)*xifunc(mij)-Vfunc(mij)-mij*xifunc(mij))
da_weight<- da_weight1*da_weight2 
d_weight11<- diag(as.vector(1/phi*(1-delta)*da_weight))

weight_star2<- dbind(d_weight11,-1*diag(rep(1/alpha^2,qq)))


#weight_star2<- dbind(0*weight11,-1*diag(rep(1/alpha^2,qq)))

gamma1<- -alpha*sum(hinv*(t(TT)%*%weight_star2%*%TT))
alpha<-t(v_h0)%*%(v_h0)/(qq-gamma1)

#################### phi_h #############
old_phi<-phi

a=(delta*mij)+(1-delta)*Vfunc(mij);
term=t(z)%*%(a+We%*%mij);

dvp<- -0.5*sqrt(1/phi)*iW%*%term
dp_weight1<- -0.5*mij/phi - 1/sqrt(phi)*z %*% dvp

dp_weight<- dp_weight1*da_weight2 
#d_weight<- -0.5*mij/phi*(2*Vfunc(mij)*xifunc(mij)-Vfunc(mij)-mij*xifunc(mij))


d_weight_phi<-1/phi*(1-delta)*dp_weight
weight111<-diag(-rep(1,nn)/phi^2)*diag(as.vector(weight1))+diag(as.vector(d_weight_phi))
weight_star<- dbind(weight111,0*weight2)

gamma0<- phi*sum(hinv*(t(TT)%*%weight_star%*%TT))
#phi<-t(ystar-muij)%*%(ystar-muij)/(sum(weight1)-pp-qq+gamma0)
phi<-t(ystar-muij)%*%(ystar-muij)/(sum(weight1)+gamma0)

phi<-as.vector(phi)
alpha<-as.vector(alpha)

beta_h<-beta_vv[1:pp]
v_h<-beta_vv[pp+1:qq]

convergence<-sum(abs(old_beta-beta_h0))+sum(abs(old_alpha-alpha))+sum(abs(old_phi-phi))

i<-i+1
##############################################################

se_beh<-se_beta_vv[1:pp]
t_value=beta_h/se_beh
p_value=2*(1-pnorm(abs(t_value)))

phi_h<-phi
alpha_h<-alpha

#print(alpha_h)


#print("beta : ")
#print(beta_vv[1:pp])
#print("SE of beta : ")
#print(se_beta_vv[1:pp])
#print("alpha : ")
#print(alpha)
#print("phi : ")
#print(phi)
}

print("iterations : ")
print(i)
print("convergence : ")
print(convergence)

F.Est<-cbind(beta_h,se_beh, t_value, p_value)  # Est. of fixed effects
D.Est<-cbind(alpha_h,phi_h)     # Est. of dispersion paras
V.Est<-as.vector(v_h)

#print(F.Est)
#print(D.Est)
#print("V.Est")
 #print(V.Est)


ress1<-list(F.Est=F.Est,D.Est=D.Est,V.Est=V.Est)

print("Estimates for fixed effects")
print_mat1<-ress1$F.Est
colnames(print_mat1)<-c("Estimate","Std. Error", "t_value", "p_value")
#rownames(print_mat1) <- c(namesX1,namesX2)
print(round(print_mat1,digits=5))
print("Estimates for dispersion parameters")

print_mat2<-ress1$D.Est
colnames(print_mat2)<-c("alpha_h", "phi_h" )
print(round(print_mat2,digits=5))
 
# print("Estimates for random effects")
# print(V.Est)


return(ress1)
}
