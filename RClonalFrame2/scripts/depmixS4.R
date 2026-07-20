#this is a script doing forward-backward in depmixS4,
#just to see how quick it is compared to my own implementation
library(depmixS4)
rho=0.1
delta=500
nu=0.01
response=as.factor(runif(1e6)<1e-4)
len=1e-4
ex=exp(-len*4/3)
pid =1/4+3/4*ex
pnid=1/4-1/4*ex
ex_nu=exp(-4*nu/3)
pid_nu=0.25+0.75*ex_nu
pnid_nu=0.25-0.25*ex_nu
a=1-exp(-rho*len)
if (a==0) b=0.5 else b=1-exp(-rho*len/expm1(rho*len*delta))#Using M/M/inf queue
m=depmix(response~1,nstates=2,family=binomial(link='identity'),ntimes=length(response))
m=setpars(m,c(1-a,a,b,1-b,pid,pnid,pid_nu,pnid_nu))#this is probably not correct
a=forwardbackward(m)#Local decoding
a$logLike#Likelihood
a$gamma#Local decoding state probabilities
#b=viterbi(m)$state #Global decoding
