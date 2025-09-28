library(survival)

Sepable_cox=function(p.time,z.level,data,setting="median",level.set=NULL){
#############################################
z1=z.level[1]
z2=z.level[2]

pLA=function(yy,tt,LL){
pla=function(yi){
loc=sum(yi>=tt)
if(loc==0)ans=0
if(loc>0)ans=LL[loc]
ans
}
apply(matrix(yy),1,pla)
}
vLA=function(yy,tt,LL){
pla=function(yi){
loc=sum(yi>=tt)-1
if(loc<=0)ans=0
if(loc>0)ans=LL[loc]
ans
}
apply(matrix(yy),1,pla)
}
############################################
MSCOX=function(data){

T1=data$X1
T2=data$X2
d2=data$D
d1=(T1<T2)*1
Z=data$S
X=NULL
PP=0
if(length(data$covariate)>0 )PP=NCOL(data$covariate)
if(PP>0) X=as.matrix(data$covariate)
m=length(T1)

########################################
MM=as.matrix(cbind(Z,X))
Xo2=list()
for(i in 1:m){Xo2[[i]]=MM[i,]%*%t(MM[i,])}

BB1=coxph(Surv(T1,d1)~MM)
st1=sort(T1,index.return=TRUE)
sd1=d1[st1$ix]
MMst1=MM[st1$ix,]
st1=st1$x

exB1=exp(apply(t(MMst1)*BB1$coef,2,sum))
dLA1=sd1/(cumsum(exB1[m:1])[m:1])
LA1=cumsum(dLA1)

LA1_i=pLA(T1,st1,LA1)
exB1_i=exp(apply(t(MM)*BB1$coef,2,sum))


A1i=LA1_i*exB1_i
UB11=0

for(i in 1:m){
UB11=UB11-as.matrix(Xo2[[i]])*A1i[i]
}

UB1L1=NULL
et1=st1[sd1>0]
met1=length(et1)
for(i in 1:met1){
UB1L1_i= -apply(MM*exB1_i*(T1>=et1[i]),2,sum)
UB1L1=cbind(UB1L1,UB1L1_i)
}

########################################
T02=ifelse(d1==1,T1,T2)
D02=ifelse(d1==1,0,d2)
BB2=coxph(Surv(T02,D02)~MM)

st2=sort(T02,index.return=TRUE)
sd2=D02[st2$ix]
MMst2=MM[st2$ix,]
st2=st2$x

exB2=exp(apply(t(MMst2)*BB2$coef,2,sum))
dLA2=sd2/(cumsum(exB2[m:1])[m:1])
LA2=cumsum(dLA2)

LA2_i=pLA(T02,st2,LA2)
exB2_i=exp(apply(t(MM)*BB2$coef,2,sum))

A2i=LA2_i*exB2_i
UB22=0

for(i in 1:m){
UB22=UB22-as.matrix(Xo2[[i]])*A2i[i]
}

UB2L2=NULL
et2=st2[sd2>0]
met2=length(et2)
for(i in 1:met2){
UB2L2_i= -apply(MM*exB2_i*(T02>=et2[i]),2,sum)
UB2L2=cbind(UB2L2,UB2L2_i)
}

########################################
LN1=T2>T1
M3=cbind(Z,X)[LN1,]
D03=(d1*d2)[LN1]
T3=T2[LN1]
R3=T1[LN1]
BB3=coxph(Surv(R3,T3,D03)~M3)
LA3=basehaz(BB3,centered=FALSE)

st3=sort(T2,index.return=TRUE)
sd3=(d1*d2)[st3$ix]
sr3=T1[st3$ix]
MMst3=MM[st3$ix,]
st3=st3$x
m3=length(st3)
exB3=exp(apply(t(MMst3)*BB3$coef,2,sum))
R3m=matrix(sr3,m3,m3)
T3m=matrix(st3,m3,m3,byrow=TRUE)
Risk3=ifelse((R3m<=T3m)&(T3m<=t(T3m)),1,0)

dLA3=sd3/apply( Risk3*exB3,2,sum)
LA3=cumsum(dLA3)

LA3_2i=pLA(T2,st3,LA3)
LA3_1i=pLA(T1,st3,LA3)
LA3_i=LA3_2i-LA3_1i
exB3_i=exp(apply(t(MM)*BB3$coef,2,sum))

A3i=LA3_i*exB3_i
UB33=0

for(i in 1:m){
UB33=UB33-as.matrix(Xo2[[i]])*A3i[i]
}

UB3L3=NULL
et3=st3[sd3>0]
met3=length(et3)
for(i in 1:met3){
UB3L3_i= -apply(MM*exB3_i*(T1<=et3[i]&et3[i]<=T2),2,sum)
UB3L3=cbind(UB3L3,UB3L3_i)
}


########################################

FIM=matrix(0,(PP+1)*3+met1+met2+met3,(PP+1)*3+met1+met2+met3)
FIM[1:(PP+1),((PP+1)*3+1):((PP+1)*3+met1)]=UB1L1
FIM[(PP+1+1):(2*(PP+1)),((PP+1)*3+met1+1):((PP+1)*3+met1+met2)]=UB2L2
FIM[(2*(PP+1)+1):(3*(PP+1)),(3*(PP+1)+1+met1+met2):(3*(PP+1)+met1+met2+met3)]=UB3L3

FIM=FIM+t(FIM)

FIM[1:(PP+1),1:(PP+1)]=UB11
FIM[(PP+1+1):(2*(PP+1)),(PP+1+1):(2*(PP+1))]=UB22
FIM[(2*(PP+1)+1):(3*(PP+1)),(2*(PP+1)+1):(3*(PP+1))]=UB33

LL=-c(rep(0,(PP+1)*3),1/dLA1[sd1>0],1/dLA2[sd2>0],1/dLA3[sd3>0])^2
FIM=FIM+diag(LL)

FIM=-(FIM/m)
JM=solve(FIM)

#diag(JM)[1:12]^0.5/sqrt(m)
EX=MX=NULL
if(PP>0){
EX=apply(MM,2,mean)[-1]
MX=apply(MM,2,median)[-1]
}
}
report=list(
dLA1=dLA1,
dLA2=dLA2,
dLA3=dLA3,
st1=st1,st2=st2,st3=st3,
sd1=sd1,sd2=sd2,sd3=sd3,JM=JM,
###################
B1=BB1$coef,
B2=BB2$coef,
B3=BB3$coef,
EX=EX,
MX=MX
##############
)
return(report) 
}

#######################################################

ans=MSCOX(data)
PX=NULL #setting="median"
if(setting=="center")PX=ans$EX
if(setting=="median")PX=ans$MX
if(setting=="zero")PX=ans$MX*0
if(setting=="level")PX=level.set

eb1=exp( sum(PX*ans$B1[-1]))
eb2=exp( sum(PX*ans$B2[-1]))
eb3=exp( sum(PX*ans$B3[-1]))

dLA1=ans$dLA1*eb1
LA1=cumsum(dLA1)
dLA2=ans$dLA2*eb2
LA2=cumsum(dLA2)
dLA3=ans$dLA3*eb3

LA1_st2=pLA(ans$st2,ans$st1 ,LA1)


Comp=function(z1,z2){
eb1_z1=exp(z1*ans$B1[1])
eb2_z2=exp(z2*ans$B2[1])
cumsum(dLA2*eb2_z2*exp(-LA2*eb2_z2-LA1_st2*eb1_z1))
}

LA2_st1=pLA(ans$st1,ans$st2 ,LA2)
LA3=cumsum(dLA3)
LA3_st1=pLA(ans$st1,ans$st3 ,LA3)


In_semi=function(z1,z2){
eb1_z1=exp(z1*ans$B1[1])
eb2_z2=exp(z2*ans$B2[1])
eb3_z2=exp(z2*ans$B3[1])

cumsum(
dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2)
)
}

Semi=function(z1,z2){
Semi_inn=In_semi(z1,z2)
eb3_z2=exp(z2*ans$B3[1])
cumsum(
dLA3*eb3_z2* vLA(ans$st3,ans$st1,Semi_inn)* exp(-LA3*eb3_z2)
)
}

CountF=function(z1,z2){
C02=Comp(z1,z2)
S03=Semi(z1,z2)
pLA(p.time,ans$st2,C02)+pLA(p.time,ans$st3,S03)
}
S10=CountF(z1,z2)
DE=S10-CountF(z1,z1)
IE=CountF(z2,z2)-S10
infB=list(B1=ans$B1,B2=ans$B2,B3=ans$B3)

report=list(DE=DE,IE=IE,infB)

return(report)
}

#Sepable_cox(p.time=c(0,0.1,0.2),c(0,1),data)

