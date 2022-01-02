#Confounders, common cause and mediators

#import modules

library(dagitty)
library(car)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
} 

#Imagine you want to study what variables influence the Z variable, having as well
#X and Y as covariates
#You may want to consider the DAG:
comm.cause.DAG <- dagitty("dag {
X -> Y
X -> Z
Y -> Z
e_y -> Y
e_z -> Z
}")

coordinates(comm.cause.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                    y = c(Y = 1, X = 3, Z = 1, e_y = 0.75, e_z = 0.75))
drawdag(comm.cause.DAG)

#As well as the DAG:
X.is.cause.DAG <- dagitty("dag {
X -> Y
X -> Z
e_y -> Y
e_z -> Z
}")

coordinates(X.is.cause.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                    y = c(Y = 1, X = 3, Z = 1, e_y = 0.75, e_z = 0.75))
drawdag(X.is.cause.DAG)

#How can you know which is the case? How do you know if you want to consider X, Y or both?
#Let's analyze all the possibilities!

#define the numerical relationship between variables, and size of the dataset (N)
N <- 500
b_xy <- 10
b_xz <- 3
b_yz <- 2
e_x <- 1
e_y <- 1
e_z <- 1

#create dataset in which X enhances Z, but Y decreases it
set.seed(13)
comm.cause.df1 <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
comm.cause.df1$Y <- comm.cause.df1$X * b_xy + rnorm(N, sd = e_y)
comm.cause.df1$Z <- comm.cause.df1$X * b_xz - comm.cause.df1$Y * b_yz + rnorm(N, sd = e_z)

#there are three ways to adjust this model
summary(lm(Z~X+Y, data = comm.cause.df1)) #both significant, only direct effect of Z
summary(lm(Z~X, data = comm.cause.df1)) #significant, total effect of X
summary(lm(Z~Y, data = comm.cause.df1)) #significant, direct effect of Y why 0.3 less even if you change it ?

b_yz <- 0

comm.cause.df1 <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
comm.cause.df1$Y <- comm.cause.df1$X * b_xy + rnorm(N, sd = e_y)
comm.cause.df1$Z <- comm.cause.df1$X * b_xz - comm.cause.df1$Y * b_yz + rnorm(N, sd = e_z)

#there are three possibilities to analyze the causal relations
summary(lm(Z~X+Y, data = comm.cause.df1)) #only X is significant
summary(lm(Z~X, data = comm.cause.df1)) #significant, total effect of X
summary(lm(Z~Y, data = comm.cause.df1)) #significant, 0.3













#create a dataset in which Z is enhanced by both X and Y
comm.cause.df2 <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
comm.cause.df2$Y <- comm.cause.df2$X * b_xy + rnorm(N, sd = e_y)
comm.cause.df2$Z <- comm.cause.df2$X * b_xz + comm.cause.df2$Y * b_yz + rnorm(N, sd = e_z)

#how does each of the scenarios behave in the adjustment of Y?
summary(lm(Z~X+Y, data = comm.cause.df2))
summary(lm(Z~X, data = comm.cause.df2))
summary(lm(Z~Y, data =comm.cause.df2))


