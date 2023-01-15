# Comparison 5 populations Binomial
library(readxl)
library(ggplot2)
library(ContourFunctions)
library(latex2exp)
library(rootSolve)

set.seed(1)
k_1 <- 12
k_2 <- 10
k_3 <- 36
k_4 <- 36
k_5 <- 37
n_1 <- 60
n_2 <- 40
n_3 <- 80
n_4 <- 60
n_5 <- 50

# zeta function : profile function
posterior_theta <- function(theta){
  post = 
    dbeta(x = theta, shape1 = k_1 + 0.5, 
          shape2 = n_1 - k_1 + 0.5)*
    dbeta(x = theta, shape1 = k_2 + 0.5, 
          shape2 = n_2 - k_2 + 0.5)*
    dbeta(x = theta, shape1 = k_3 + 0.5, 
          shape2 = n_3 - k_3 + 0.5)*
    dbeta(x = theta, shape1 = k_4 + 0.5, 
          shape2 = n_4 - k_4 + 0.5)*
    dbeta(x = theta, shape1 = k_5 + 0.5, 
          shape2 = n_5 - k_5 + 0.5)
  return(post)
}

# max
theta = seq(0.001,1,0.00001)
x11()
plot(theta,posterior_theta(theta), type = "l")
L = max(posterior_theta(theta))
L
P_comp = theta[which.max(posterior_theta(theta))]
P_comp 

P_star <- rep(P_comp, 5)

# h distribution
posterior_TOT <- function(theta){
  post = dbeta(x = theta[1], shape1 = k_1 + 0.5, 
               shape2 = n_1 - k_1 + 0.5)*
         dbeta(x = theta[2], shape1 = k_2 + 0.5, 
               shape2 = n_2 - k_2 + 0.5)*
         dbeta(x = theta[3], shape1 = k_3 + 0.5, 
               shape2 = n_3 - k_3 + 0.5)*
    dbeta(x = theta[4], shape1 = k_4 + 0.5, 
          shape2 = n_4 - k_4 + 0.5)*
    dbeta(x = theta[5], shape1 = k_5 + 0.5, 
          shape2 = n_5 - k_5 + 0.5)
  
  return(post)
}

# Optimal vector
c_star = gradient(f = posterior_TOT, x = P_star)
c_star = c_star/sqrt(sum(c_star^2))
c_star
sum(c_star)
round(c_star,3)

psi = c_star[1]*(rbeta(n = 1000, shape1 = k_1 + 0.5, 
                               shape2 = n_1 - k_1 + 0.5))+
      c_star[2]*(rbeta(n = 1000, shape1 = k_2 + 0.5, 
                               shape2 = n_2 - k_2 + 0.5))+
      c_star[3]*(rbeta(n = 1000, shape1 = k_3 + 0.5, 
                               shape2 = n_3 - k_3 + 0.5))+
      c_star[4]*(rbeta(n = 1000, shape1 = k_4 + 0.5, 
                               shape2 = n_4 - k_4 + 0.5))+
      c_star[5]*(rbeta(n = 1000, shape1 = k_5 + 0.5, 
                               shape2 = n_5 - k_5 + 0.5))

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
# 1
