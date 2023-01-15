# Comparison 4 means Gaussians
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

# # Sampling of the parameter vector theta = (phi,mu)
# # for the 4 populations
set.seed(1)
posterior_mu <- function(mu){
  post = dt2(x = mu,
             location = h1$eta,
             scale =  sqrt(h1$beta/(h1$nu*h1$alpha)),
             df = 2*h1$alpha)*
    dt2(x = mu,
        location = h2$eta,
        scale =  sqrt(h2$beta/(h2$nu*h2$alpha)),
        df = 2*h2$alpha)*
    dt2(x = mu,
        location = h3$eta,
        scale =  sqrt(h3$beta/(h3$nu*h3$alpha)),
        df = 2*h3$alpha)*
    dt2(x = mu,
        location = h4$eta,
        scale =  sqrt(h4$beta/(h4$nu*h4$alpha)),
        df = 2*h4$alpha)
  
  return(post)
}

# max
mu = seq(40,70,0.000001)
L = max(posterior_mu(mu))
L
P_comp = mu[which.max(posterior_mu(mu))]
P_comp 
P_star = rep(P_comp, 4)
P_star

# Plot zeta function
x11()
plot(mu, posterior_mu(mu), type = "l")

# h distribution
posterior_TOT <- function(mu){
  post = dt2(x = mu[1],
             location = h1$eta,
             scale =  sqrt(h1$beta/(h1$nu*h1$alpha)),
             df = 2*h1$alpha)*
    dt2(x = mu[2],
        location = h2$eta,
        scale =  sqrt(h2$beta/(h2$nu*h2$alpha)),
        df = 2*h2$alpha)*
    dt2(x = mu[3],
        location = h3$eta,
        scale =  sqrt(h3$beta/(h3$nu*h3$alpha)),
        df = 2*h3$alpha)*
    dt2(x = mu[4],
        location = h4$eta,
        scale =  sqrt(h4$beta/(h4$nu*h4$alpha)),
        df = 2*h4$alpha)
 
  return(post)
}

# Optimal vector
c_star = gradient(f = posterior_TOT, 
                  x = P_star, 
                  pert = 1e-11)
c_star
c_star = c_star/sqrt(sum(c_star^2))
c_star
sum(c_star)
round(c_star,3)

psi = c_star[1]*(rt2(n = 10000,
                    location = h1$eta,
                    scale =  sqrt(h1$beta/(h1$nu*h1$alpha)),
                    df = 2*h1$alpha))+
      c_star[2]*(rt2(n = 10000,
                    location = h2$eta,
                    scale =  sqrt(h2$beta/(h2$nu*h2$alpha)),
                    df = 2*h2$alpha))+
      c_star[3]*(rt2(n = 10000,
                    location = h3$eta,
                    scale =  sqrt(h3$beta/(h3$nu*h3$alpha)),
                    df = 2*h3$alpha))+
      c_star[4]*(rt2(n = 10000,
                    location = h4$eta,
                    scale =  sqrt(h4$beta/(h4$nu*h4$alpha)),
                    df = 2*h4$alpha))

int <- sum(psi>0)/length(psi)
int

delta_H <- d_H(int)
delta_H

# ciclisti
# p-value = 0.00359
# delta_H = 0.9392 


