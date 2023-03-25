
rm(list=ls())
library(stable)

n=100
delta=0;gamma=1;beta=0;alpha=1
pr=0 #1 OR 0 parametrization
theta1 <- 2
theta0 <- 5

p <- 2 # number of predictors
h <- floor(n/2) + floor((p+1)/2)+floor(1-(2/alpha));h # Breakdown point

x=runif(n,0,100)
e=rstable(n,alpha,beta,gamma,delta,pr)
y=theta0+(theta1*x)+e
all.data <- data.frame(cbind(y,x))

#================== ols ==========================
initial_OLS <- lm(y~x, data = all.data)

res.ols <- initial_OLS$residuals

cutp= stable.fit.mle.restricted(res.ols, c(alpha,beta,1,0), restriction=c(0,0,1,1), param=pr)

if(round(cutp[2],4) < 1 & round(cutp[2],4)> -1){
  c=ceiling(2/cutp[1]);d=floor(n+1-(2/cutp[1]))
}
if(round(cutp[1],4)>=1 & round(cutp[2],4)==1 |round(cutp[2],4)==-1){
  c=ceiling(2/cutp[1]);d=floor(n+1-(2/cutp[1]))
}
if(round(cutp[1],4)<1 & round(cutp[2],4)==1 |round(cutp[2],4)==-1){
  c=ceiling(2/cutp[1]);d=floor(n+1-(2/cutp[1]))
}
c;d
#c=2
#d=99

ord.data <- all.data[order(e),]
trimmed.data <- ord.data[c:d,]
rownames(trimmed.data)<- array(1:(d-c+1))

#================LTS based on OLS ===============================
rep=500
A <- matrix(ncol=p+1,nrow=rep)
B <- matrix(ncol=p,nrow=rep)
C <- c()
D <- matrix(ncol=h,nrow=rep)
for (J in 1:rep) {
  m.p <- sample(1:n,p+1)
  p.trim.1 <- trimmed.data[m.p,]
  m1.p <- lm(y ~ x, data=p.trim.1)
  e1.p <- trimmed.data$y - predict(m1.p, newdata=trimmed.data)
  e1.p.s <- sort(abs(e1.p))
  H0 <- as.numeric(names(e1.p.s[1:h]))
  
  d.trim.1 <- trimmed.data[H0,]
  m1.h <- lm(y ~ x, data=d.trim.1)
  e1 <- trimmed.data$y - predict(m1.h, newdata=trimmed.data)
  e1.s <- sort(abs(e1))
  
  repeat{
    s1=sum(e1[H0]^2);s1
    H0 <- as.numeric(names(e1.s[1:h]))
    h.trim.1 <- trimmed.data[H0,]
    m1.h <- lm(y ~ x, data=h.trim.1)
    e1 <- trimmed.data$y - predict(m1.h, newdata=trimmed.data)
    e1.s <- sort(abs(e1))
    s2=sum(e1[H0]^2);s2
    if(s1==s2){
      break
    }
  }
  
  A[J,]= as.numeric(m.p)
  B[J,]=m1.h$coefficients
  C[J]=s2
  D[J,]=H0
}
result <- as.data.frame(cbind(A,B,C,D))
result.lts <- result[order(result[,6]),]
lts.ols=as.numeric(result.lts[1,4:5]) 
lts.ols

# 4.846134 2.008244
#===============================================
lts.ols.res <- trimmed.data$y - lts.ols[1]- (lts.ols[2]*trimmed.data$x)
names(lts.ols.res) <- array(1:(d-c+1))
lts.ols.res.s <- sort(abs(lts.ols.res))
hi.res.LTS.ols <-   lts.ols.res.s[1:h]

H1 <- as.numeric(names(hi.res.LTS.ols))
init.lts=stable.fit.mle.restricted(lts.ols.res, c(alpha,beta,1,0), restriction=c(0,0,1,1), param=pr)

MLE.LTS=function(par){
  y=ord.data[,1]
  x=ord.data[,2]
  bet1=par[1]
  bet2=par[2]
  i=H1[!is.na(H1)]
  f3=function(i,beta1,beta2){log(factorial(d-c+1)/(factorial(i-1)*factorial(d-c+1-i)))+log(dstable(y[i]-beta1-(beta2*x[i]),init.lts[1],init.lts[2],1,0,pr))+(i-1)*log(pstable(y[i]-beta1-(beta2*x[i]),init.lts[1],init.lts[2],1,0,pr))+(n-i)*log(1-pstable(y[i]-beta1-(beta2*x[i]),init.lts[1],init.lts[2],1,0,pr))}
  g=-(sum(sapply(i,f3,beta1=bet1,beta2=bet2)))
}

theta.lts <- nlminb(lts.ols,MLE.LTS)
theta.lts$par

#[1] 5.081401 2.000729

