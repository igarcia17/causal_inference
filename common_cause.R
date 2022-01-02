#Confounder case: common cause

#import modules

library(dagitty)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
} 

#Consider the DAG:
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

#define the numerical relationship between variables, and size (N)
N <- 500
b_xy <- 10
b_xz <- 3
b_yz <- 2
e_x <- 1
e_y <- 1
e_z <- 1

#create a dataset in which Z is enhanced by both X and Y
set.seed(13)
comm.cause.df1 <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
comm.cause.df1$Y <- comm.cause.df1$X * b_xy + rnorm(N, sd = e_y)
comm.cause.df1$Z <- comm.cause.df1$X * b_xz + comm.cause.df1$Y * b_yz + rnorm(N, sd = e_z)

#create dataset in which X enhances Z, but Y decreases it
comm.cause.df2 <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
comm.cause.df2$Y <- comm.cause.df2$X * b_xy + rnorm(N, sd = e_y)
comm.cause.df2$Z <- comm.cause.df2$X * b_xz - comm.cause.df2$Y * b_yz + rnorm(N, sd = e_z)

#how does each of the scenarios behave in the adjustment of Y?

summary(lm(Z~X+Y, data = comm.cause.df1))
summary(lm(Z~X, data = comm.cause.df1))

summary(lm(Z~X+Y, data = comm.cause.df2))
summary(lm(Z~X, data = comm.cause.df2))


