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

#create dataset in which both X and Y positively increase Z
set.seed(13)
rows <- 500
comm.cause.df1 <- data.frame(X = runif(rows, 1, 100))
comm.cause.df1$Y <- comm.cause.df1$X * 10 + rnorm(rows)
comm.cause.df1$Z <- comm.cause.df1$X * 3 + comm.cause.df1$Y * 2 + rnorm(rows)

#create dataset in which X enhances Z, but Y decreases it
comm.cause.df2 <- data.frame(X = runif(rows, 1, 100))
comm.cause.df2$Y <- comm.cause.df2$X * 10 + rnorm(rows)
comm.cause.df2$Z <- comm.cause.df2$X * 3 - comm.cause.df2$Y * 2 + rnorm(rows)


