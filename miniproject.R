########Q3
####(b)
a <- c(3.68,1.683,1.402,1.385,1.149,1.102,0.9324,0.7592,0.7482,0.6201)
a <- a/100
n <- 8046
N <- c(176,191,119,134,149,48,50,75,40,60)
l_temp <- c()
for (i in 1:10){
  temp <- ((a[i]*sum(N)/n/sum(a))/(N[i]/n))^N[i]
  l_temp <- c(l_temp,temp)
}
l_x <- prod(l_temp)
T <- -2*log(l_x)
pchisq(T,df=9)


####(c)
set.seed(6008)
p <- a*sum(N)/n/sum(a)
p <- p/sum(p)
len_sample <- sum(N)
B <- 10000
resample_array <- array(0,dim = c(B,len_sample))
for (i in 1:B){
  resample_array[i,] <- sample(1:10,len_sample,replace = TRUE, prob = p)
}
bootstrap_T <- c()

for(i in 1:B){
  bootstrap_N <- c()
  for (j in 1:10){
    temp_N <- length(which(resample_array[i,]==j))
    bootstrap_N <- c(bootstrap_N,temp_N)
  }
  bootstrap_T <- c(bootstrap_T,-2*sum(bootstrap_N * (log(a*sum(bootstrap_N))-log(bootstrap_N*sum(a)))))
}

quantile(bootstrap_T,prob = 0.95)
quantile(bootstrap_T,prob = 0.9999999999999)

####(d)
bootstrap_cdf <- ecdf(bootstrap_T)
plot(bootstrap_cdf,main='cdf')
x <- seq(0,35, length=10000) 
y <- pchisq(x,df=9)
lines(x,y,col=3)
legend('bottomright',pch=c(15,15),legend = c('chisq','bootstrap'),col=c(1,3))







########Q4
####(b)
otherN_length <- n-sum(N)  #Remaining length except the first 10 N_i
ff <- function(zz){        #Define f(Z)
  p <- c()
  for (i in 1:10){
    p_temp <- N[i]*zz/otherN_length
    pi <- max(p_temp,a[i])
    p <- c(p,pi)
  }
  f <- sum(p) + zz - 1
  return(f)
}
root <- uniroot(ff,c(0.5,1),tol=0.01)
z <- root$root            #Z = 0.85
p_H0 <- c()
for (i in 1:10){
  pi <- max(N[i]*z/otherN_length,a[i])
  p_H0 <- c(p_H0,pi)
}
p_H1 <- N/n
ratio_first_10 <- 1       #likelihood ratio of the first 10 term
temp <- p_H0/p_H1
for (j in 1:10){
  ratio_first_10 <- ratio_first_10*(temp[j]^N[j])
}
T <- -2*log(ratio_first_10*(n*z/otherN_length)^otherN_length)  #102.25

####(c)
p_H0 == a


####(d)
set.seed(6008)
p_resample <- c(p_H0,1-sum(p_H0))                            #probability of resample
len_sample <- n                                              
B <- 10000
resample_array <- array(0,dim = c(B,len_sample))             #define an array to store resamples
for (i in 1:B){
  resample_array[i,] <- sample(1:11,len_sample,replace = TRUE, prob = p_resample)
}
bootstrap_T <- c()                                           #define an variable to store results of bootstrap statistic
#The following ff1 defines f(Z). There's an undefined variable bootstrap_N in ff1. This variable represents the length of each N_i in each resample.
#I'll define bootstrap_N later
ff1 <- function(zz1){                                        
  p <- c()
  for (i in 1:10){
    p_temp <- bootstrap_N[i]*zz1/bootstrap_N[11]
    pi <- max(p_temp,a[i])
    p <- c(p,pi)
  }
  f <- sum(p) + zz1 - 1
  return(f)
}

for(k in 1:B){
  bootstrap_N <- c()                                       #define boostrap_N to store the results of each resample
  for (j in 1:11){
    temp_N <- length(which(resample_array[k,]==j))
    bootstrap_N <- c(bootstrap_N,temp_N)
  }
  i <- 0
  z <- uniroot(ff1,c(0,1),tol = 0.0001)$root                 #compute Z for B'th resample
  p_H0 <- c()                                              #the following procedure is similar to (b)
  for (l in 1:10){
    pi <- max(bootstrap_N[l]*z/bootstrap_N[11],a[l])
    p_H0 <- c(p_H0,pi)
  }
  p_H1 <- bootstrap_N/n
  ratio_first_10 <- 1
  temp <- p_H0/p_H1[1:10]
  for (m in 1:10){
    ratio_first_10 <- ratio_first_10*temp[m]^bootstrap_N[m]
  }
  T <- -2*log(ratio_first_10*(n*z/bootstrap_N[11])^bootstrap_N[11])
  bootstrap_T <- c(bootstrap_T,T)
}

bootstrap_cdf <- ecdf(bootstrap_T)
plot(bootstrap_cdf,main='cdf')
quantile(bootstrap_T,prob = 0.95)







