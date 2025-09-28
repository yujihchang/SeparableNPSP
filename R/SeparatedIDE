## Function Name
# SeparatedIDE

## Purpose
# This function analyzes differences in the cumulative distribution function (difference in cdf) under semicompeting risks or competing risks settings.
# It utilizes the Cox proportional hazards model to calculate cumulative risks under different exposure levels.

## Arguments

### p.time
# - Description: A vector of time points at which the user is interested in calculating the difference in cdf.
# - Format: A numeric vector. For example, `c(12, 24, 36)` specifies calculations at 12, 24, and 36 months.
# - Default: None; the user must provide this input.


### data
# - Description: A matrix containing the following columns:
#   - `X1`: Observed intermediate event times.
#   - `X2`: Observed primary event times.
#   - `D1`: Status indicator for the intermediate event (1 if the event occurred, 0 otherwise).
#   - `D2`: Status indicator for the primary event (1 if the event occurred, 0 otherwise).
#   - `S`: Exposure levels; setting za=1;zb=2.


### type
# - Description: Specifies the analysis method.
#   - `"SP"`: Semicompeting risks method (default).
#   - `"CP"`: Competing risks method.
# - Default: `"SP"`.



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


############################################################################
############################################################################
############################################################################
SeparatedIDE=function(P.time,data,type){
#P.time=sort(data$X1);type="SP"
typeS=names(table(data$S))
Z=data$S
Z[data$S==typeS[1]]=0
Z[data$S==typeS[2]]=1
#####################################
#####################################
x.obs_0=data$X1[Z==0]
y.obs_0=data$X2[Z==0]
dy_0=data$D2[Z==0]
dx_0=data$D1[Z==0]
#---------------------------------------------------
x.obs_1=data$X1[Z==1]
y.obs_1=data$X2[Z==1]
dy_1=data$D2[Z==1]
dx_1=data$D1[Z==1]

z.obs=ifelse(data$X1<=data$X2,data$X1,data$X2)
dz=data$D
z.obs_0=z.obs[Z==0]
z.obs_1=z.obs[Z==1]
dz_0=dz[Z==0]
dz_1=dz[Z==1]


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


tpye_CS=(type=="SP")*1

library("expm")
Xu2010_0=function(T1,T2,d1,d2){
#T1=x.obs_1;T2=y.obs_1;d2=dy_1



m=length(T1)

T2M=matrix(T2,m,m)
#######for U2
ST1=sort(T1,index.return=TRUE)
st1=ST1$x
T1M=matrix(T1,m,m)
ST1M=t(matrix(st1,m,m))
Risk1=ifelse(T1M>=ST1M,1,0)
sd1=d1[ST1$ix]
d1i=ifelse(T1M==ST1M,1,0)*matrix(d1,m,m)
###########for U3
ST23=sort(T2,index.return=TRUE)
sd23=( d2*(1-d1) )[ST23$ix]
st23=ST23$x
ST23M=t(matrix(st23,m,m))
Risk2=ifelse(T1M>=ST23M,1,0)
d23i=ifelse(T2M==ST23M,1,0)*matrix(d2*(1-d1),m,m)
##############for U4
ST24=sort(T2,index.return=TRUE)

sd24=(d1*d2)[ST24$ix]*tpye_CS
st24=ST24$x

ST24M=t(matrix(st24,m,m))
Risk12=ifelse((T2M>=ST24M) & (ST24M>T1M),1,0)
d24i=ifelse(T2M==ST24M,1,0)*matrix(d1*d2,m,m)
pLA=function(yy,tt,LL){
pla=function(yi){
loc=sum(yi>=tt)
if(loc==0)ans=0
if(loc>0)ans=LL[loc]
ans
}
apply(matrix(yy),1,pla)
}

Nelson_Aalen=function(yy,dd,risk_na){
####################################
#LA1=st1
#LA2=st2
#dL1=diff( c(0,LA1) )*sd1
#dL2=diff( c(0,LA2) )*sd23
#dL3=diff( c(0,LA2) )*sd24
####################################
down_na=apply(risk_na,2,sum)
down_na=ifelse(down_na==0,1,down_na)
dL=dd/down_na
dL
}
dL1=Nelson_Aalen(st1,sd1,Risk1)
dL2=Nelson_Aalen(st23,sd23,Risk2)
dL3=Nelson_Aalen(st24,sd24,Risk12)*tpye_CS

############################################
#estimating 
############################################

LA1=pLA(T1,st1,cumsum(dL1))
LA2=pLA(T1,st23,cumsum(dL2))
LA32=pLA(T2,st24,cumsum(dL3))
LA31=pLA(T1,st24,cumsum(dL3))
LA3=LA32-LA31
Ai=LA1+LA2+LA3

ETA=c(dL1,dL2,dL3) 
#################################################################

down1=apply(Risk1,2,sum)
down2=apply(Risk2,2,sum)
down3=apply(Risk12,2,sum)

#################################################################
Score=function(ETA){
dL1=ETA[1               :m           ] 
dL2=ETA[(m+1)        :(m+m)]
dL3=ETA[(m+m+1):(m+m+m)]*tpye_CS
u2=sd1/dL1-down1
u3=sd23/dL2-down2
u4=sd24/dL3-down3
ans=c(u2,u3,u4*tpye_CS)
ans
}

##################################################################
SS=Score(ETA)
LPC=which(ETA!=0)
eee=(1/m)^2
ee=diag(eee,length(ETA))
UU=NULL
app_score=function(ii){( Score(ETA+ee[ii,])-SS )/eee}
uu=apply(matrix(LPC),1,app_score)
selectNA=apply(uu,1,sum)
UU=uu[selectNA !="NaN",]
UU=0.5*(UU+t(UU))
JM=solve(-UU)


########################################
report=list(
dL1=dL1,
dL2=dL2,
dL3=dL3*tpye_CS,
st1=st1,st2=st23,
d1=sd1,d23=sd23,d24=sd24,JM=JM,
###################
dL1_z=dL1,
dL2_z=dL2,
dL3_z=dL3,
t1_z=st1,
t2_z=st23,
t3_z=st24,
d1_z=sd1,
d2_z=sd23,
d3_z=sd24
##############
)
return(report)
}

#T1=x.obs_1;T2=y.obs_1;d2=dy_1
#T1=x.obs_0;T2=y.obs_0;d2=dy_0
Rcumsum=function(x){rev(cumsum(rev(x)))}

Ans_1=ans_Z2=Xu2010_0(x.obs_1,T2=y.obs_1,d1=dx_1,d2=dy_1)
Ans_0=ans_Z1=Xu2010_0(x.obs_0,T2=y.obs_0,d1=dx_0,d2=dy_0)

dLa1_z1=ans_Z1$dL1
dLa2_z1=ans_Z1$dL2
dLa3_z1=ans_Z1$dL3
dLa1_z2=ans_Z2$dL1
dLa2_z2=ans_Z2$dL2
dLa3_z2=ans_Z2$dL3


t1_z1=ans_Z1$st1
t2_z1=ans_Z1$st2

t1_z2=ans_Z2$st1
t2_z2=ans_Z2$st2

d1_z1=ans_Z1$d1
d1_z2=ans_Z2$d1
d2_z1=ans_Z1$d23
d2_z2=ans_Z2$d23
d3_z1=ans_Z1$d24
d3_z2=ans_Z2$d24
################################################################
################################################################
################################################################

# competed part_10
lambda_02_z1=Ans_1$dL2_z
LA1_j_z0=pLA(Ans_1$t2_z,Ans_0$t1_z,cumsum(Ans_0$dL1_z))
LA2_j_z1=pLA(Ans_1$t2_z,Ans_1$t2_z,cumsum(Ans_1$dL2_z))

d_competed_part_10=lambda_02_z1*exp(-(LA1_j_z0+LA2_j_z1))
S_competed_part_10=pLA(P.time,Ans_1$t2,cumsum(d_competed_part_10))

# competed part_11
#
LA1_j_z1=pLA(Ans_1$t2_z,Ans_1$t1_z,cumsum(Ans_1$dL1_z))

d_competed_part_11=lambda_02_z1*exp(-(LA1_j_z1+LA2_j_z1))
S_competed_part_11=pLA(P.time,Ans_1$t2,cumsum(d_competed_part_11))

# competed part_00
lambda_02_z0=Ans_0$dL2_z
LA1_j_z0=pLA(Ans_0$t2_z,Ans_0$t1_z,cumsum(Ans_0$dL1_z))
LA2_j_z0=pLA(Ans_0$t2_z,Ans_0$t2_z,cumsum(Ans_0$dL2_z))

d_competed_part_00=lambda_02_z0*exp(-(LA1_j_z0+LA2_j_z0))
S_competed_part_00=pLA(P.time,Ans_0$t2,cumsum(d_competed_part_00))
#plot(P.time,S_competed_part_00)
###########################################################

#semi_part_A_10_j
lambda_021_z1=Ans_1$dL3_z
LA3_j_z1=pLA(Ans_1$t3_z,Ans_1$t3_z,cumsum(Ans_1$dL3_z))
dsemi_part_A_10=lambda_021_z1*exp(-LA3_j_z1)

#semi_part_B_10_v
LA3_v0_z1=pLA(Ans_0$t1_z,Ans_1$t3_z,cumsum(Ans_1$dL3_z)) 
LA2_v0_z1=pLA(Ans_0$t1_z,Ans_1$t2_z,cumsum(Ans_1$dL2_z)) 
LA1_v0_z0=pLA(Ans_0$t1_z,Ans_0$t1_z,cumsum(Ans_0$dL1_z)) 
lambda_01_z0=Ans_0$dL1_z
S_B_10_v=exp(-(LA1_v0_z0+LA2_v0_z1-LA3_v0_z1)) * tpye_CS
semi_part_B_10_v=pLA(Ans_1$t3_z,Ans_0$t1_z,cumsum(lambda_01_z0*S_B_10_v))

dsemi_part_AB_10= dsemi_part_A_10* semi_part_B_10_v  
Ssemi_part_AB_10=pLA(P.time,Ans_1$t3_z,cumsum(dsemi_part_AB_10)) *  tpye_CS

###########################################################

#semi_part_A_00_j
lambda_021_z0=Ans_0$dL3_z 
LA3_j_z0=pLA(Ans_0$t3_z,Ans_0$t3_z,cumsum(Ans_0$dL3_z))
dsemi_part_A_00=lambda_021_z0*exp(-LA3_j_z0) 
#semi_part_B_00_v

LA3_v0_z0=pLA(Ans_0$t1_z,Ans_0$t3_z,cumsum(Ans_0$dL3_z))
LA2_v0_z0=pLA(Ans_0$t1_z,Ans_0$t2_z,cumsum(Ans_0$dL2_z)) 
LA1_v0_z0=pLA(Ans_0$t1_z,Ans_0$t1_z,cumsum(Ans_0$dL1_z)) 
S_B_00_v=exp(-(LA1_v0_z0+LA2_v0_z0-LA3_v0_z0))
semi_part_B_00_v=pLA(Ans_0$t3_z,Ans_0$t1_z,cumsum(lambda_01_z0*S_B_00_v))

dsemi_part_AB_00= dsemi_part_A_00* semi_part_B_00_v * tpye_CS 
Ssemi_part_AB_00=pLA(P.time,Ans_0$t3_z,cumsum(dsemi_part_AB_00))  *  tpye_CS

###########################################################
#semi_part_A_11_j
#
dsemi_part_A_11=dsemi_part_A_10

#semi_part_B_11_v

LA3_v1_z1=pLA(Ans_1$t1_z,Ans_1$t3_z,cumsum(Ans_1$dL3_z))
LA2_v1_z1=pLA(Ans_1$t1_z,Ans_1$t2_z,cumsum(Ans_1$dL2_z))
LA1_v1_z1=pLA(Ans_1$t1_z,Ans_1$t1_z,cumsum(Ans_1$dL1_z))
lambda_01_z1=Ans_1$dL1_z
S_B_11_v=exp(-(LA1_v1_z1+LA2_v1_z1-LA3_v1_z1))
semi_part_B_11_v=pLA(Ans_1$t3_z,Ans_1$t1_z,cumsum(lambda_01_z1*S_B_11_v))*  tpye_CS

dsemi_part_AB_11= dsemi_part_A_11* semi_part_B_11_v
Ssemi_part_AB_11=pLA(P.time,Ans_1$t3_z,cumsum(dsemi_part_AB_11))*  tpye_CS



S00=S_competed_part_00+Ssemi_part_AB_00
S10=S_competed_part_10+Ssemi_part_AB_10
S11=S_competed_part_11+Ssemi_part_AB_11

DE=S10-S00
IE=S11-S10




################################################################
################################################################
################################################################


part_S=function(tt,z2,z1){#z1=0;z2=1;tt=0.25
if(z1==0){dLa1=dLa1_z1;t1=t1_z1;d1=d1_z1}
if(z1==1){dLa1=dLa1_z2;t1=t1_z2;d1=d1_z2}
if(z2==0){dLa2=dLa2_z1;dLa3=dLa3_z1;t2=t2_z1;d2=d2_z1;d3=d3_z1}
if(z2==1){dLa2=dLa2_z2;dLa3=dLa3_z2;t2=t2_z2;d2=d2_z2;d3=d3_z2}

La1_j=pLA(t2,t1,cumsum(dLa1))
La2_j=pLA(t2,t2,cumsum(dLa2))
La3_j=pLA(t2,t2,cumsum(dLa3))

La1_v=pLA(t1,t1,cumsum(dLa1))
La2_v=pLA(t1,t2,cumsum(dLa2))
La3_v=pLA(t1,t2,cumsum(dLa3))

t1u=t1[d1>0]

nt2=length(t2)
nt1=length(t1)

nt1u=length(t1u)

t2M=matrix(t2,nt2,nt1u)
t1Mu=matrix(t1u,nt2,nt1u,byrow=TRUE)

CL=(tt>=t2)*1
#################################################################
#partial dLa1
Risk_ujM=ifelse(t1Mu<=t2M,1,0)

#j-version
Au=-apply(Risk_ujM*matrix(dLa2*exp( -(La1_j+La2_j) )*CL,nt2,nt1u),2,sum)

#v-version
t2M_all=matrix(t2,nt2,nt1)
t1M_all=matrix(t1,nt2,nt1,byrow=TRUE)

cv=matrix(dLa1*exp( -(La1_v+La2_v-La3_v) ),nt2,nt1,byrow=TRUE)
Risk_uvM=ifelse(t1M_all<t2M_all,1,0)
Bu=-apply(((t(apply(cv*Risk_uvM,1,Rcumsum)))[,d1>0])*dLa3*exp(-La3_j)*CL,2,sum)

Euv=matrix(exp(-(La1_v+La2_v-La3_v)),nt2,nt1,byrow=TRUE)
Cu=apply((Euv*Risk_uvM)*dLa3*exp(-La3_j)*CL,2,sum)[d1>0]

pla1=Au+Bu+Cu
#####################################################################
#partial dLa2
t2u=t2[d2>0]
CLu=ifelse(t2u<=tt,1,0)
nt2u=length(t2u)
La1_u=pLA(t2u,t1,cumsum(dLa1))
La2_u=pLA(t2u,t2,cumsum(dLa2))

AAu=exp(-(La1_u+La1_u))*CLu

BBU_0=t(apply(cv*Risk_uvM,1,Rcumsum))
BBu_1=BBU_0*dLa3*exp(-La3_j)*CL
BBu=apply(BBu_1,1,sum)[d2>0]

t2Mu=matrix(t2u,nt2,nt2u,byrow=TRUE)
t2M_dLa2=matrix(t2,nt2,nt2u)

Risk_ujM_dLa2=ifelse(t2Mu<=t2M_dLa2,1,0)

CCu=-apply(Risk_ujM_dLa2*matrix(dLa2*exp( -(La1_j+La2_j) )*CL,byrow=TRUE,nt2,nt2u),2,sum)

pla2=AAu+BBu+CCu
#####################################################################
#partial dLa3

AAAu=(apply((cv*Risk_uvM),1,sum)*exp(-La3_j)*CL)[d3>0]

Risk_uvM_3=ifelse(t1M_all<=t2M_all,1,0)
BBBu=-apply((t(apply(cv*Risk_uvM_3,1,Rcumsum)))*dLa3*exp(-La3_j)*CL,1,sum)[d3>0]

pla3=AAAu+BBBu


if(z2==0&z1==0){ans=c(pla1,pla2,pla3,rep(0,sum(c(d1_z2,d2_z2,d3_z2))))}
if(z2==1&z1==1){ans=c(rep(0,sum(c(d1_z1,d2_z1,d3_z1))),pla1,pla2,pla3)}
if(z2==1&z1==0){ans=c(pla1,rep(0,sum(c(d2_z1,d3_z1))),rep(0,sum(d1_z2)),pla2,pla3)}

ans
}

JMz2=diag(ans_Z2$JM)
JMz1=diag(ans_Z1$JM)
JM=diag(c(JMz1,JMz2))

DEsd=IEsd=rep(NA,length(P.time))
m=dim(data)[1]
for(ii in 1:length(P.time)){
part21=part_S(P.time[ii],z2=1,z1=0)
PSSDE=part21-part_S(P.time[ii],z2=0,z1=0)
PSSIE=part_S(P.time[ii],z2=1,z1=1)-part21
DEsd[ii]=(t(PSSDE)%*%JM%*%PSSDE)^0.5
IEsd[ii] =(t(PSSIE) %*%JM%*%PSSIE)^0.5
}



report_ans=list(P.time=P.time,DE=DE,IE=IE,
S10=S10,S00=S00,S11=S11,
Comp_part_IE=S_competed_part_11-S_competed_part_10,
Semi_part_IE=Ssemi_part_AB_11-Ssemi_part_AB_10,
DEsd=DEsd,
IEsd=IEsd
)
return(report_ans)
}



#P.time=seq(0,1.5,by=0.01)
#ans=SeparatedIDE(P.time,data)
#plot(P.time,ans$DE,type="l")
#plot(P.time,ans$DE,type="l")

