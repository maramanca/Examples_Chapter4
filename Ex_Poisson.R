# Comparison 4 parameters Poisson
library(readxl)
library(ggplot2)
library(ContourFunctions)
library(latex2exp)
library(rootSolve)

# Example 1
s_1 = 34
s_2 = 4
s_3 = 6
s_4 = 5
n_1 = 154
n_2 = 19
n_3 = 20
n_4 = 11

# Example 2 
s_1 = 4
s_2 = 1
s_3 = 7
s_4 = 9
n_1 = 1259
n_2 = 2082
n_3 = 1417
n_4 = 1647

# zeta function : profile function
posterior_lambda <- function(lambda){
  post = dgamma(x = lambda, shape = 0.5+s_1, rate = n_1)*
    dgamma(x = lambda, shape = 0.5+s_2, rate = n_2)*
    dgamma(x = lambda, shape = 0.5+s_3, rate = n_3)*
    dgamma(x = lambda, shape = 0.5+s_4, rate = n_4)
  return(post)
}

# Solution with the max
lambda = seq(0.001,4,0.00001)
x11()
plot(lambda,posterior_lambda(lambda), type = "l")
L = max(posterior_lambda(lambda))
L
P_comp = lambda[which.max(posterior_lambda(lambda))]
P_comp 

# Solution with the derivative
A = s_1+s_2+s_3+s_4-2
A
P_comp = A/(n_1+n_2+n_3+n_4)
P_comp

P_star = rep(P_comp, 4)
P_star

# h distribution
posterior_TOT <- function(lambda){
  post = dgamma(x = lambda[1], shape = 0.5+s_1, rate = n_1)*
         dgamma(x = lambda[2], shape = 0.5+s_2, rate = n_2)*
         dgamma(x = lambda[3], shape = 0.5+s_3, rate = n_3)*
         dgamma(x = lambda[4], shape = 0.5+s_4, rate = n_4)
  return(post)
}

# Optimal vector
c_star = gradient(f = posterior_TOT, x = P_star)
c_star = c_star/sqrt(sum(c_star^2))
c_star
sum(c_star)
round(c_star,3)

psi = c_star[1]*(rgamma(n = 1000, shape = 0.5+s_1, rate = n_1))+
  c_star[2]*(rgamma(n = 1000, shape = 0.5+s_2, rate = n_2))+
  c_star[3]*(rgamma(n = 1000, shape = 0.5+s_3, rate = n_3))+
  c_star[4]*(rgamma(n = 1000, shape = 0.5+s_4, rate = n_4))

int <- sum(psi>0)/length(psi)
int

# Discrepancy measure
d_H <- function(Int){
  
  if(Int < 0.5){
    out <- 1 - 2*Int
  }
  else {
    out <- 1 - 2*(1 - Int)
  }
}

delta_H <- d_H(int)
delta_H



