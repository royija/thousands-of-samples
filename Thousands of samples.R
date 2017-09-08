#Calling Packages:
#=================
library("psych", lib.loc="~/R/win-library/3.4")
#=====================
#Functions Definition
#=====================
TKL<-function(x,h)
{
  l<-length(x)
  x1<-t (rbind(x,seq(1 : l )))
  x1<-x1[order(x1[ ,1],x1[,2],decreasing=TRUE), ]
  return (length (intersect(x1[1:h,2],seq(from=1,to=h,by=1)))/h)
}
TKL2<-function(x,y,K)
{
  l<-length(x)
  x1<-t(rbind(x,seq(1:l)))
  x1<-x1[order(x1[,1],x1[,2],decreasing=TRUE), ]
  y1<-t(rbind(y,seq(1:l)))
  y1<-y1[order(y1[,1],y1[,2],decreasing=TRUE), ]
  return(length(intersect(x1[1:K,2],y1[1:K,2]))/K)
}

#===================================================
#Direct v.s. Approximated Approach (Estimated sig-q)
#===================================================

N=20000 #Number of genes
alpha<-0.005 #Proportion of the genes that are correlated with Y
u<-alpha*N
K<-seq(from=600,to=1200,by=200)
l<-length(K)
m<-NA
tkl<-NA
U<-N-u
beta<-c(rep(1,times=u),rep(0,times=U))
c1<-NA
f1<-NA
v2<-NA

#approximated model
draw<-NA
empiric<-NA
p<-NA
sig_q<-NA
tkl2<-NA
c<-NA
f<-NA
r<-NA
v<-NA

for (t in 1:l)
{
  k<-K[t]-1
  for(i in 1:10)
  {
    x<-rnorm(n=N)
    y<-sum(x[1:u])
    c<-NA
    for (j in 1:k)
    {
      z<-rnorm(n=N)
      x<-rbind(x,z)
      y=rbind(y,sum(z[1:u]))
    }
    for (j in 1:N)
    {
      r<-cor(x[,j],y)
      c[j]<-fisherz(r)
      c1[j]<-abs(r)
    }
    p<-sqrt(1/(K[t]-3))
    sig_q<-sqrt(var(c)-p^2)
    draw<-rnorm(n=N,mean=0,sd=sig_q)
    empiric<-abs(draw+rnorm(n=N,mean=0,sd=p))
    draw<-abs(draw)
    f[i]<-TKL2(draw,empiric,u)
    f1[i]<-TKL(c1,u)
  }
  
  tkl[t]<-mean(f1)
  v[t]<-sd(f1)
  tkl2[t]<-mean(f)
  v2[t]<-sd(f)  
}  

plot(tkl~K,type='b',col='blue',xlab='Sample???Size',ylab='Proportion',xlim=c(600,1200), ylim=c(0,0.7))
lines(tkl2~K,type='b',col='red')

#=================
#Histogram+QQ???plot
#=================

#u=100

N=20000 #Number o f genea
alpha=0.005 #Proportion of the genes that are correlated with Y
u<-alpha*N
K<-seq(from=30,to=600,by=30)
l<-length(K)
c<-NA
f<-NA
m<-NA

tkl<-NA
U<-N-u
beta<-c(rep(1,times=u),rep(0,times=U))
t<-2
k<-K[t]-1
x<-rnorm(n=N)
y<-sum(x[1:u])
c<-NA
for(j in 1:k)
{
  z<-rnorm(n=N)
  x<-rbind(x,z)
  y=rbind(y,sum(z[1:u]))
}
for (j in 1:N)
{
  c[j]<-fisherz(cor(x[,j],y))
}
par(mfrow=c(1,2))
hist(c,xlab='Fisher???Z')
qqnorm(c)

