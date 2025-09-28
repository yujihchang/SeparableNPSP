
#data=meta.gen(n=500,L1=1,L2=1,L3=3,b01=0,b02=0,b03=0.5,cc1=2,cc2=2)
#z.level=c(0,1);p.time=seq(0,1.5,by=0.01);setting="median";type="CP"
#z.level=c(0,1);data=data_F;setting="user";level=lv_F;type="SP"


Sepable_cox=function(p.time,z.level,data,setting="median",level=NULL,type="SP"){
#############################################
BB3=sd3=st3=dLA3=NULL


######################
z1=z.level[1]
z2=z.level[2]
m=dim(data)[1]

type_SC=(type=="SP") *1

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

pLA_M=function(yy,tt,LL){
if(NCOL(LL)==1)LL=matrix(LL)
pla=function(yi){
loc=sum(yi>=tt)
if(loc==0)ans=0*LL[1,]
if(loc>0)ans=LL[loc,]
as.numeric(ans)
}
apply(matrix(yy),1,pla)
}

vLA_M=function(yy,tt,LL){
if(NCOL(LL)==1)LL=matrix(LL)
pla=function(yi){
loc=sum(yi>=tt)-1
if(loc<=0)ans=0*LL[1,]
if(loc>0)ans=LL[loc,]
as.numeric(ans)
}
apply(matrix(yy),1,pla)
}

#rcumsum=function(x){ x_L=length(x);cumsum(x[x_L:1])[x_L:1]  }

cut_cum=function(tt1,tt2,LL){
M1=matrix(tt1,length(tt1),length(tt2))
M2=matrix(tt2,length(tt1),length(tt2),byrow=TRUE)
cut.LL=( (M2<=M1) )*LL
apply(cut.LL,2,cumsum) }

pin_time=function(p.time,tt,LL){
vv=rep(0,length(LL))
Mp.time=matrix(p.time,length(tt),length(p.time),byrow=TRUE)
Mtt=matrix(tt, length(tt),length(p.time))
(Mp.time>=Mtt)*LL
}



############################################
MSCOX=function(data){

T1=data$X1
T2=data$X2
d2=data$D2
d1=data$D1
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
met3=0
if(sum(LN1)>0){
M3=cbind(Z,X)[LN1,]
D03=(d1*d2)[LN1]* type_SC
T3=T2[LN1]
R3=T1[LN1]
BB3=coxph(Surv(R3,T3,D03)~M3)
LA3=basehaz(BB3,centered=FALSE)

st3=sort(T2,index.return=TRUE)
sd3=(d1*d2)[st3$ix]* type_SC

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
if(met3>0){
for(i in 1:met3){
UB3L3_i= -apply(MM*exB3_i*(T1<=et3[i]&et3[i]<=T2),2,sum)
UB3L3=cbind(UB3L3,UB3L3_i)
}
}
}

########################################


FIM=matrix(0,(PP+1)*3+met1+met2+met3,(PP+1)*3+met1+met2+met3)
FIM[1:(PP+1),((PP+1)*3+1):((PP+1)*3+met1)]=UB1L1
FIM[(PP+1+1):(2*(PP+1)),((PP+1)*3+met1+1):((PP+1)*3+met1+met2)]=UB2L2

if(met3>0){
FIM[(2*(PP+1)+1):(3*(PP+1)),(3*(PP+1)+1+met1+met2):(3*(PP+1)+met1+met2+met3)]=UB3L3
}

FIM=FIM+t(FIM)

FIM[1:(PP+1),1:(PP+1)]=UB11
FIM[(PP+1+1):(2*(PP+1)),(PP+1+1):(2*(PP+1))]=UB22
if(met3>0){
FIM[(2*(PP+1)+1):(3*(PP+1)),(2*(PP+1)+1):(3*(PP+1))]=UB33
}

if(met3>0) LL=-c(rep(0,(PP+1)*3),1/dLA1[sd1>0],1/dLA2[sd2>0],1/dLA3[sd3>0])^2  else LL=-c(rep(0,(PP+1)*3),1/dLA1[sd1>0],1/dLA2[sd2>0])^2

FIM=FIM+diag(LL)
Det_FIM= which(apply(FIM,2,sum)==0)
if(length(-Det_FIM)>0){
FIM=FIM[-Det_FIM,-Det_FIM]
}

FIM=-(FIM/m)
JM=solve(FIM)

EX=MX=NULL
if(PP>0){
EX=apply(MM,2,mean)[-1]
MX=apply(MM,2,median)[-1]
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
MX=MX,pp=PP
##############
)
return(report) 
}

#######################################################

ans=MSCOX(data)
PP=ans$pp
PX=NULL #setting="median"
if(setting=="center")PX=ans$EX
if(setting=="median")PX=ans$MX
if(setting=="zero")PX=ans$MX*0
if(setting=="user")PX=level

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
if(type=="SP") S03=Semi(z1,z2) else S03=0*(ans$st3)
pLA(p.time,ans$st2,C02)+pLA(p.time,ans$st3,S03)
}
S10=CountF(z1,z2)
DE=S10-CountF(z1,z1)
IE=CountF(z2,z2)-S10
########################################################
# Delta method
########################################################

#######################################
# competing risks
#######################################

DeltaF_Comp_B=function(z1,z2){
pp=PP
eb1_z1=exp(z1*ans$B1[1])
eb2_z2=exp(z2*ans$B2[1])
#beta's
####B1####
pB1= - matrix(c(z1,PX),length(ans$st2),pp+1,byrow=TRUE)*cumsum( LA1_st2*eb1_z1*  dLA2*eb2_z2*exp(-LA2*eb2_z2-LA1_st2*eb1_z1) )

####B2####
pB2=
 matrix( c(z2,PX),length(ans$st2),pp+1,byrow=TRUE)*cumsum(                  dLA2*eb2_z2*exp(-LA2*eb2_z2-LA1_st2*eb1_z1) )-
 matrix( c(z2,PX),length(ans$st2),pp+1,byrow=TRUE)*cumsum( LA2*eb2_z2*      dLA2*eb2_z2*exp(-LA2*eb2_z2-LA1_st2*eb1_z1) )

####B3####
#pB3=0
if(type=="SP") pB3=pB2*0 else pB3=NULL
cbind(pB1,pB2,pB3)

}

#######################################
# semi-competing risks
#######################################

DeltaF_semi_B=function(z1,z2){
eb1_z1=exp(z1*ans$B1[1])
eb2_z2=exp(z2*ans$B2[1])
eb3_z2=exp(z2*ans$B3[1])
###########
pB1_Inn=
t( c(z1,PX) )%x%cumsum(             dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2))-
t( c(z1,PX) )%x%cumsum( LA1*eb1_z1* dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2))

if(NCOL(pB1_Inn)==1){
pB1=cumsum(pLA_M(ans$st3,ans$st1,pB1_Inn)* dLA3*eb3_z2 * exp(-LA3*eb3_z2))
}

if(NCOL(pB1_Inn)>1){
pB1=apply(t(pLA_M(ans$st3,ans$st1,pB1_Inn))* dLA3*eb3_z2 * exp(-LA3*eb3_z2),2,cumsum)
}

############
############

pB2_Inn=
 - t( c(z2,PX))%x%cumsum( LA2_st1*eb2_z2* dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2))

if(NCOL(pB2_Inn)==1){
pB2=cumsum(pLA_M(ans$st3,ans$st1,pB2_Inn)*dLA3*eb3_z2* exp(-LA3*eb3_z2))
}
if(NCOL(pB2_Inn)>1){
pB2=apply(t(pLA_M(ans$st3,ans$st1,pB2_Inn))*dLA3*eb3_z2* exp(-LA3*eb3_z2),2,cumsum)
}
##############
##############
##############

pB3_Inn_0=
 t( c(z2,PX))%x%cumsum(                dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2))
pB3_Inn_1=
 t( c(z2,PX))%x%cumsum(LA3_st1*eb3_z2* dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2))

if(NCOL(pB3_Inn_0)==1){
pB3=
cumsum( 
  pLA_M(ans$st3,ans$st1,pB3_Inn_0)*           dLA3*eb3_z2 * exp(-LA3*eb3_z2)+
 -pLA_M(ans$st3,ans$st1,pB3_Inn_0)*LA3*eb3_z2*dLA3*eb3_z2 * exp(-LA3*eb3_z2)+
 +pLA_M(ans$st3,ans$st1,pB3_Inn_1)*           dLA3*eb3_z2 * exp(-LA3*eb3_z2)
)}

if(NCOL(pB3_Inn_0)>1){
pB3=
apply( 
  t(pLA_M(ans$st3,ans$st1,pB3_Inn_0))*           dLA3*eb3_z2 * exp(-LA3*eb3_z2)+
 -t(pLA_M(ans$st3,ans$st1,pB3_Inn_0))*LA3*eb3_z2*dLA3*eb3_z2 * exp(-LA3*eb3_z2)+
 +t(pLA_M(ans$st3,ans$st1,pB3_Inn_1))*           dLA3*eb3_z2 * exp(-LA3*eb3_z2), 
 2,cumsum)}
###########
cbind(pB1,pB2,pB3)
}


#####################
# delta method for lambda	
#####################

DeltaF_Comp_L=function(p.time,z1,z2){
eb1_z1=exp(z1*ans$B1[1])
eb2_z2=exp(z2*ans$B2[1])
eb3_z2=exp(z2*ans$B3[1])

LL=dLA2*eb2_z2*exp(-LA2*eb2_z2-LA1_st2*eb1_z1)
LL_1=eb2_z2*exp(-LA2*eb2_z2-LA1_st2*eb1_z1)
####################

PL1_0=-cut_cum(ans$st2[ans$sd2>0],ans$st1[ans$sd1>0],(LL[ans$sd2>0] * eb1_z1*eb1) ) #dim(PL1)
PL1=pLA_M(p.time,ans$st2[ans$sd2>0],PL1_0)
####################

PL2_1=pin_time(p.time,ans$st2[ans$sd2>0],LL_1[ans$sd2>0]*eb2) 

PL2_20=-cut_cum(ans$st2[ans$sd2>0],ans$st2[ans$sd2>0],(LL[ans$sd2>0]*eb2_z2*eb2 ))  #dim(PL2_2)
PL2_2=pLA_M(p.time,ans$st2[ans$sd2>0],PL2_20)
PL2=PL2_1+PL2_2
####################

PL3=matrix(0,sum(ans$sd3),length(p.time)) #dim(PL3)
rbind(PL1,PL2,PL3)
}
    #dim(DeltaF_Comp_L(p.time,1,2))

########################################################
########################################################
########################################################

DeltaF_semi_L=function(p.time,z1,z2){

eb1_z1=exp(z1*ans$B1[1])
eb2_z2=exp(z2*ans$B2[1])
eb3_z2=exp(z2*ans$B3[1])

LL=dLA1*eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2)
LL_1=   eb1_z1*exp(-LA2_st1*eb2_z2-LA1*eb1_z1+LA3_st1*eb3_z2)

############################################

PL1_Inn_10=pin_time(ans$st3,ans$st1[ans$sd1>0],LL_1[ans$sd1>0]*eb1)
PL1_Inn_11=t(PL1_Inn_10)*dLA3*eb3_z2* exp(-LA3*eb3_z2)
PL1_Inn_12=apply(PL1_Inn_11,2,cumsum)
PL1_1=pLA_M(p.time,ans$st3,PL1_Inn_12)

PL1_Inn_20=-cut_cum(ans$st1[ans$sd1>0],ans$st1[ans$sd1>0],( LL[ans$sd1>0] * eb1_z1*eb1) )
PL1_Inn_21=t( pLA_M(ans$st3[ans$sd3>0],ans$st1[ans$sd1>0],PL1_Inn_20) )
PL1_Inn_22=PL1_Inn_21 * (dLA3*eb3_z2* exp(-LA3*eb3_z2))[ans$sd3>0]
PL1_Inn_23=apply(PL1_Inn_22,2,cumsum)
PL1_2=pLA_M(p.time,ans$st3[ans$sd3>0],PL1_Inn_23)

PL1=PL1_1+PL1_2
############################################

PL2_Inn_20=-cut_cum(ans$st1[ans$sd1>0],ans$st2[ans$sd2>0],( LL[ans$sd1>0] * eb2_z2*eb2) )
PL2_Inn_21=t(pLA_M(ans$st3[ans$sd3>0],ans$st1[ans$sd1>0],PL2_Inn_20))
PL2_Inn_22=PL2_Inn_21*( (dLA3*eb3_z2* exp(-LA3*eb3_z2) )[ans$sd3>0] )
PL2_Inn_23=apply(PL2_Inn_22,2,cumsum)
PL2=pLA_M(p.time,ans$st3[ans$sd3>0],PL2_Inn_23)

##############################################

PL3_Inn_10=cumsum(LL)
PL3_Inn_11=pLA(ans$st3,ans$st1[ans$sd1>0],PL3_Inn_10[ans$sd1>0])
PL3_Inn_12=(PL3_Inn_11*eb3_z2*eb3* exp(-LA3*eb3_z2))[ans$sd3>0]
PL3_1=pin_time(p.time,ans$st3[ans$sd3>0],PL3_Inn_12)

PL3_Inn_20=cumsum(LL)
PL3_Inn_21=pLA(ans$st3[ans$sd3>0],ans$st1[ans$sd1>0],PL3_Inn_20[ans$sd1>0])
PL3_Inn_22=PL3_Inn_21*dLA3[ans$sd3>0] * eb3_z2*eb3* exp(-LA3[ans$sd3>0]*eb3_z2)
PL3_Inn_23=-cut_cum(ans$st3[ans$sd3>0],ans$st3[ans$sd3>0],PL3_Inn_22)
PL3_2=pLA_M(p.time,ans$st3[ans$sd3>0],PL3_Inn_23)

PL3_Inn_30=cut_cum(ans$st1[ans$sd1>0],ans$st3[ans$sd3>0],(LL[ans$sd1>0] * eb3_z2*eb3))
PL3_Inn_31=t(pLA_M(ans$st3[ans$sd3>0],ans$st1[ans$sd1>0],PL3_Inn_30))
PL3_Inn_32=PL3_Inn_31*dLA3[ans$sd3>0] * eb3_z2* exp(-LA3[ans$sd3>0]*eb3_z2)
PL3_Inn_33=apply(PL3_Inn_32,2,cumsum)
PL3_3=pLA_M(p.time,ans$st3[ans$sd3>0],PL3_Inn_33)

PL3=PL3_1+PL3_2+PL3_3
rbind(PL1,PL2,PL3)
}

DeltaF_B=function(p.time,z1,z2){  
if(type=="SP") S_part=pLA_M(p.time,ans$st3,DeltaF_semi_B(z1,z2)) else S_part=0
                           pLA_M(p.time,ans$st2,DeltaF_Comp_B(z1,z2))+  S_part
                        
                           }


DeltaF_L=function(p.time,z1,z2){
if(type=="SP") S_part=DeltaF_semi_L(p.time,z1,z2) else S_part=0
S_part+DeltaF_Comp_L(p.time,z1,z2)
}


Delta_allbind=function(p.time,z1,z2){
rbind(DeltaF_B(p.time,z1,z2),DeltaF_L(p.time,z1,z2))
}
DD10=Delta_allbind(p.time,z1,z2)
DD00=Delta_allbind(p.time,z1,z1)
DD11=Delta_allbind(p.time,z2,z2)

np=length(p.time)
DE.var=IE.var=rep(NA,np)
VM=ans$JM
for(pptt in 1:np){
dDE=DD10[,pptt]-DD00[,pptt]
dIE=DD11[,pptt]-DD10[,pptt]
DE.var[pptt]=t(dDE)%*%VM%*%dDE/m
IE.var[pptt]=t(dIE)%*%VM%*%dIE/m
 }


report=list(
DE=DE,
IE=IE,
DE.sd=DE.var^0.5,
IE.sd=IE.var^0.5
)


return(report)
}

#Sepable_cox(p.time=c(0,0.1,0.2),c(0,1),data)

