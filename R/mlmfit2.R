mlmfit2 <-
function(jm1, data, Maxiter=300) {
 mc <- match.call()
 data_surv=data


 Vfunc <- function(x) {
   return(dnorm(x)/(1-pnorm(x)))
}

xifunc <- function(x) {
   return(Vfunc(x)*(Vfunc(x)-x))
}

 res1<-design_frailtyHL(jm1[[3]],data=data_surv)

yy<-log(res1[[1]])

### dispersion parameters #
phi<- 1
alpha1<-0.05
alpha2<-0.1

## x matrix
xx<-res1[[2]]

## z matrix
zz1<-res1[[4]][[1]]
zz2<-res1[[4]][[2]]
zz<-cbind(zz1,zz2)

#zz<-cbind(res1[[4]][[1]],res1[[4]][[2]])
#zz<-res1[[4]][[1]]

nn<-nrow(yy)
pp<-ncol(xx)+1
qq<-ncol(zz)
qq1<-ncol(zz1)
qq2<-ncol(zz2)

#qq1=13
#qq2=128

## delta

delta<-res1[[3]]

## 
 beta_h0<-rbind(0,res1[[9]])
 #beta_h0<-rbind(0,0)

##### initial values by using linear mixed model
#res_lm<-lmer(yy~treat+(1|hospital)+(1|patient),cgd1)
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
#while (convergence>0.000001) {
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

weight21<-diag(rep(1/alpha1,qq1))
weight22<-diag(rep(1/alpha2,qq2))
weight<-dbind(dbind(weight11,weight21),weight22)

#weight2<-diag(rep(1/alpha1,qq))
#weight<-dbind(weight11,weight2)

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

v_h01<- v_h0[1:qq1]
v_h02<- v_h0[qq1+1:qq2]


#################### alpha1_h #############
old_alpha1<-alpha1

We<-diag(as.vector(weight1))

lam1<-as.vector(phi/alpha1)
i_q1<-diag(rep(1,qq1))
U1<-lam1*i_q1

lam2<-as.vector(phi/alpha2)
i_q2<-diag(rep(1,qq2))
U2<-lam2*i_q2

U=dbind(U1,U2)

iW<- solve(t(z)%*%We%*%z +U)  
iD10<- -diag(rep(1/alpha1^2,qq1))
iD11<- -diag(rep(0,qq2))
iD1<-dbind(iD10,iD11)
dva1<- -phi*iW%*%iD1%*%v_h0

da_weight1<- -1/sqrt(phi)*z %*% dva1
da_weight2 <- (2*Vfunc(mij)*xifunc(mij)-Vfunc(mij)-mij*xifunc(mij))
da_weight<- da_weight1*da_weight2 
d_weight11<- diag(as.vector(1/phi*(1-delta)*da_weight))

weight_star2<- dbind(d_weight11,iD1)

#weight_star2<- dbind(dbind(d_weight11, iD10), iD11)
#weight_star2<- dbind(d_weight11,-1*diag(rep(1/alpha^2,qq)))

gamma1<- -alpha1*sum(hinv*(t(TT)%*%weight_star2%*%TT))
alpha1<-t(v_h01)%*%(v_h01)/(qq1-gamma1)

#alpha<-t(v_h0)%*%(v_h0)/(qq-gamma1)

#################### alpha2_h #############
old_alpha2<-alpha2

#U=dbind(U1,U2)
#iW<- solve(t(z)%*%We%*%z +U) 

iD20<- -diag(rep(0,qq1))
iD21<- -diag(rep(1/alpha2^2,qq2))
iD2<-dbind(iD20,iD21)
dva2<- -phi*iW%*%iD2%*%v_h0

da_weight1_2<- -1/sqrt(phi)*z %*% dva2
#da_weight2 <- (2*Vfunc(mij)*xifunc(mij)-Vfunc(mij)-mij*xifunc(mij))
da_weight_2<- da_weight1_2*da_weight2 
d_weight11_2<- diag(as.vector(1/phi*(1-delta)*da_weight_2))

weight_star2_2<- dbind(d_weight11_2, iD2)
#weight_star2_2<- dbind(dbind(d_weight11_2, iD20), iD21)

gamma2<- -alpha2*sum(hinv*(t(TT)%*%weight_star2_2%*%TT))
alpha2<-t(v_h02)%*%(v_h02)/(qq2-gamma2)

#################### phi_h #############
old_phi<-phi

a=(delta*mij)+(1-delta)*Vfunc(mij);
term=t(z)%*%(a+We%*%mij);

dvp<- -0.5*sqrt(1/phi)*iW%*%term
dp_weight1<- -0.5*mij/phi - 1/sqrt(phi)*z %*% dvp

dp_weight<- dp_weight1*da_weight2 

d_weight_phi<-1/phi*(1-delta)*dp_weight
weight111<-diag(-rep(1,nn)/phi^2)*diag(as.vector(weight1))+diag(as.vector(d_weight_phi))
weight_star<- dbind(dbind(weight111,0*weight21),0*weight22)

#weight_star<- dbind(weight111,0*weight2)
#weight<-dbind(dbind(weight11,weight21),weight22)

gamma0<- phi*sum(hinv*(t(TT)%*%weight_star%*%TT))
phi<-t(ystar-muij)%*%(ystar-muij)/(sum(weight1)+gamma0)

phi<-as.vector(phi)
alpha1<-as.vector(alpha1)
alpha2<-as.vector(alpha2)

beta_h<-beta_vv[1:pp]
v_h<-beta_vv[pp+1:qq]

convergence<-sum(abs(old_beta-beta_h0))+sum(abs(old_alpha1-alpha1))
       +sum(abs(old_alpha2-alpha2))+sum(abs(old_phi-phi))

i<-i+1

se_beh<-se_beta_vv[1:pp]
t_value=beta_h/se_beh
p_value=2*(1-pnorm(abs(t_value)))

phi_h<-phi
alpha1_h<-alpha1
alpha2_h<-alpha2

#print(alpha_h)

#print("beta : ")
#print(beta_vv[1:pp])
#print("SE of beta : ")
#print(se_beta_vv[1:pp])
#print("alpha1 : ")
#print(alpha1)
#print("alpha2 : ")
#print(alpha2)
#print("phi : ")
#print(phi)
}
print("iterations : ")
print(i)
print("convergence : ")
print(convergence)

F.Est<-cbind(beta_h,se_beh, t_value, p_value)  # Est. of fixed effects
D.Est<-cbind(alpha1_h,alpha2_h, phi_h)     # Est. of dispersion paras
V.Est<-as.vector(v_h)

#print(F.Est)
#print(D.Est)
#print(V.Est)


res_u<-list(F.Est=F.Est,D.Est=D.Est,V.Est=V.Est)

print("Estimates for fixed effects")
print_mat1<-res_u$F.Est
colnames(print_mat1)<-c("Estimate","Std. Error", "t_value", "p_value")
print(round(print_mat1,digits=5))

print("Estimates for dispersion parameters")
print_mat2<-res_u$D.Est
colnames(print_mat2)<-c("alpha1_h", "alpha2_h", "phi_h" )
print(round(print_mat2,digits=5))

# print("Estimates for random effects")
# print(V.Est)

return(res_u)

}
