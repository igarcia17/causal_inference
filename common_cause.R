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
XYaffect.DAG <- dagitty("dag {
X -> Y
X -> Z
Y -> Z
e_y -> Y
e_z -> Z
}")

coordinates(XYaffect.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                    y = c(Y = 3, X = 1, Z = 3, e_y = 2.75, e_z = 2.75))
drawdag(XYaffect.DAG)

#As well as the DAG:
X.is.cause.DAG <- dagitty("dag {
X -> Y
X -> Z
e_y -> Y
e_z -> Z
}")

coordinates(X.is.cause.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                  y = c(Y = 3, X = 1, Z = 3, e_y = 2.75, e_z = 2.75))

drawdag(X.is.cause.DAG)

#In both cases X is common cause of Y and Z. In the first case, Y has a causal effect
#on Z, while in the second case they are independent.

#How can you know which is the case? How do you know if you want to consider X, Y or both?
#Let's analyze all the possibilities!

#Let's create a function that creates the different datasets that we need.
create.dataset <- function(b_yz, N = 500, b_xy = 3, b_xz = 3,
                           e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz - name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}

set.seed(13)
Ynoinfluences <- create.dataset(0)
Yinfluences <- create.dataset(2)
non.influences <- create.dataset(0, b_xz = 0)
#Yistruecause <- create.dataset(5, b_xz = 0)

summary(lm(Z~X+Y, data = Yinfluences)) #both significant, only direct effect of Z
summary(lm(Z~X, data = Yinfluences)) #significant, total effect of X
summary(lm(Z~Y, data = Yinfluences)) #significant, direct effect of Y why 0.3 less even if you change it ?

#there are three possibilities to analyze the causal relations
summary(lm(Z~X+Y, data = Ynoinfluences)) #only X is significant
summary(lm(Z~X, data = Ynoinfluences)) #significant, total effect of X
summary(lm(Z~Y, data = Ynoinfluences)) #significant, 0.3

#We create a function that checks which is the scenario

#The argument must be a dataframe, which column names are Z (the variable of study)
#Y and X. X should be the common cause.

Y_check <- function (dataset, conflevel = 0.01) {
  model_with_Y <- lm(Z~X+Y, data = dataset)
  p.v.X <-(summary(model_with_Y)$coefficients['X','Pr(>|t|)'])
  p.v.Y <- (summary(model_with_Y)$coefficients['Y', 'Pr(>|t|)'])
  if ((p.v.X <= conflevel)&(p.v.Y <= conflevel))
  {cat("The provided outcome of the dataset is influenced by both X and Y\n")
    cat('See plot')
    drawdag(XYaffect.DAG)}
  if ((p.v.X <= conflevel)&(p.v.Y > conflevel))
    {cat("The provided outcome of the dataset is not influenced by Y\n")
    cat('See plot')
    drawdag(X.is.cause.DAG)}
  if ((p.v.X > conflevel)&(p.v.Y > conflevel))
  {cat("It seems that neither X or Y affect Z\n You may want to review your experimental model")}
  }

Y_check(non.influences)
Y_check(Yinfluences)
Y_check(Ynoinfluences)

#Let's focus first on the case in which Y doesn't has a causal relationship with Z
#That is, the dataset Ynoinfluences

impliedConditionalIndependencies(X.is.cause.DAG)
Y.noin.condboth <- lm(Z~X + Y, data = Ynoinfluences)
Y.noin.condY <- lm(Z~Y, data = Ynoinfluences)
Y.noin.condX <- lm(Z~X, data = Ynoinfluences)
summary(Y.noin.condboth)
summary(Y.noin.condY)
summary(Y.noin.condX)

