# Confronto parametri 4 Gaussiane
library(wiqid)
library(rootSolve)

# Function that gives the hyperparameters of the posterior distribution
hyperparam <- function(n, mean, s_square){
  eta <- mean          
  nu <- n
  alpha <- (n - 1)/2
  beta <- (n*s_square)/2
  list(eta = eta, nu = nu, alpha = alpha, beta = beta)
} 

# Function that computes the BDM 
d_H <- function(Int){
  if(Int < 0.5){
    out <- 1 - 2*Int
  }
  else {
    out <- 1 - 2*(1 - Int)
  }
}

data0 <- read.table("ciclisti.txt", header = TRUE)
data <- matrix(data0$time, ncol=4, byrow=T,
               dimnames = list(subj=1:9, 
                               cond=c("0","5","9","13")))
# [1]
x_1 <- data[,1]
n_1 <- length(x_1)
x1_mean <- mean(x_1)
s_square1 <- var(x_1)

# [2]
x_2 <- data[,2]
n_2 <- length(x_2)
x2_mean <- mean(x_2)
s_square2 <- var(x_2)

# [3]
x_3 <- data[,3]
n_3 <- length(x_3)
x3_mean <- mean(x_3)
s_square3 <- var(x_3)

# [4]
x_4 <- data[,4]
n_4 <- length(x_4)
x4_mean <- mean(x_4)
s_square4 <- var(x_4)

# Hyperparameters of the 4 populations
h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)
h3 <- hyperparam(n = n_3, mean = x3_mean, s_square = s_square3)
h4 <- hyperparam(n = n_4, mean = x4_mean, s_square = s_square4)

# Sampling of the parameter vector theta = (phi,phi)
set.seed(1)
posterior_phi <- function(phi){
  post = dgamma(x = phi,
                shape = h1$alpha,
                scale = h1$beta)*
         dgamma(x = phi,
                shape = h2$alpha,
                scale = h2$beta)*
         dgamma(x = phi,
                shape = h3$alpha,
                scale = h3$beta)*
         dgamma(x = phi,
                shape = h4$alpha,
                scale = h4$beta)
  return(post)
}

# max
phi = seq(2700,2850,0.00001) 
L = max(posterior_phi(phi))
L
P_comp = phi[which.max(posterior_phi(phi))]
P_comp 
P_star = rep(P_comp, 4)
P_star

# Plot zeta function
x11()
plot(phi, posterior_phi(phi), type = "l")

# h distribution
posterior_TOT <- function(phi){
  post = dgamma(x = phi[1],
                shape = h1$alpha,
                scale = h1$beta)*
         dgamma(x = phi[2],
                shape = h2$alpha,
                scale = h2$beta)*
         dgamma(x = phi[3],
                shape = h3$alpha,
                scale = h3$beta)*
         dgamma(x = phi[4],
                shape = h4$alpha,
                scale = h4$beta)
  return(post)
}

# Optimal vector
c_star = gradient(f = posterior_TOT, 
                  x = P_star)
c_star = c_star/sqrt(sum(c_star^2))
c_star
sum(c_star)
round(c_star,3)

psi = c_star[1]*(rgamma(n = 10000,
                 shape = h1$alpha,
                 scale = h1$beta))+
      c_star[2]*(rgamma(n = 10000,
                 shape = h2$alpha,
                 scale = h2$beta))+
      c_star[3]*(rgamma(n = 10000,
                 shape = h3$alpha,
                 scale = h3$beta))+
      c_star[4]*(rgamma(n = 10000,
                 shape = h4$alpha,
                 scale = h4$beta)) 

int <- sum(psi>0)/length(psi)
int

delta_H <- d_H(int)
delta_H
# delta_H = 0.5812