lambda=0.01#arrival rate
delta=0.01#departure rate
mean_idle=1/lambda
mean_busy=expm1(lambda/delta)/lambda
cat(c(mean_idle,mean_busy))

rho=1
len=1e-3
delta=500
a=1.0-exp(-rho*len)
b=1.0-exp(-1.0/delta)
b2=1.0-exp(-rho*len/expm1(rho*len*delta))
cat(c(a,b,b2))
