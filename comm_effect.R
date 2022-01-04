# Collider and selection bias, cases:

# Import modules
library(tidyverse)
library(dagitty)
library(car)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}
# library(rethinking) ## for drawdag
# Installing rethinking can be complicated just for a few graphs
# So have a fallback if rethinking not available



# We are going to try two specific scenarios, one where we adjust for Z and
# another where we don't.
# Lets first generate our DAG according to a collider structure where X and Y
# are independent:
comm.effect.DAG <- dagitty("dag {
X -> Z
Y -> Z
e_z -> Z
}")

coordinates(comm.effect.DAG) <- list(x = c(X = 1, Y = 3, Z = 2, 
                                           e_z = 1.25),
                                     y = c(X = 1, Y = 1, Z = 3, 
                                           e_z = 3))
drawdag(comm.effect.DAG)


# Let us generate the data according to the DAG:

N <- 500 # Our sample size will be 500
b_xz <- 3 # 
b_yz <- 2
sd_z <- 5

set.seed(11)

X <- runif(N, 1, 5)
Y <- runif(N, 2, 4)
Z <- b_xz * X + b_yz * Y + rnorm(N, 0, sd = sd_z)

df_collider_ind <- data.frame(X, Y, Z)
df_collider_ind


# Let's take a look at our data distribution:

reg_line_ZX <- lm(Z ~ X)
plot(Z ~ X)
abline(reg_line_ZX)

reg_line_ZY <- lm(Z ~ Y)
plot(Z ~ Y)
abline(reg_line_ZY)

reg_line_XY <- lm(X ~ Y)
plot(X ~ Y)
abline(reg_line_XY)

# As we can see, there is a positive correlation between Z and X and also
# Z and Y. On the other hand, there is no visible correlation between X and Y 
# (as should). We can check it by calculating the estimates and significance of
# each pair:

summary(lm(Z ~ X)) #we can see a strong correlation between Z and X
summary(lm(Z ~ Y)) #we can see a strong correlation between Z and Y
summary(lm(X ~ Y)) #we can't any correlation between X and Y



summary(lm(X ~ Y + Z))


# As we can see, two independent variables (X and Y) when we adjust by Z 
# raise a correlation that wasn't supposed to be there.


#_______________________________________________________________________________
# We could also have cases in which our X and Y variables have some sort of
# relation but the estimate changes when we condition on Z.

comm.effect.DAG_2 <- dagitty("dag {
X2 -> Z2
Y2 -> Z2
X2 -> Y2
e_z -> Z2
}")

coordinates(comm.effect.DAG_2) <- list(x = c(X2 = 1, Y2 = 3, Z2 = 2, 
                                             e_z = 1.25),
                                       y = c(X2 = 1, Y2 = 1, Z2 = 3, 
                                             e_z = 3))
drawdag(comm.effect.DAG_2)


N2 <- 500
X2 <- runif(N2, 1, 10)
Y2 <- X2 * 1.7 + rnorm(N, mean = 0, sd = 0.1)
Z2 <- 1.5 * X2 + 2.5 * Y2 + rnorm(N2, mean = 0, sd = 0.1)


##Again, let's take a look at our data distribution:

reg_line_ZX2 <- lm(Z2 ~ X2)
plot(Z2 ~ X2)
abline(reg_line_ZX2)

reg_line_ZY2 <- lm(Z2 ~ Y2)
plot(Z2 ~ Y2)
abline(reg_line_ZY2)

reg_line_XY2 <- lm(X2 ~ Y2)
plot(X2 ~ Y2)
abline(reg_line_XY2)

# As we can see, there is a positive correlation between Z2 and X2; Z2 and Y2; and X2
# and Y2. In this case we are interested in the relationship between X2 and Y2; 
# we can check it: 

summary(lm(X2 ~ Y2)) # Positive significant correlation

# But if we adjust for Z2 the estimate of the X2-Y2 correlation switches to
# a negative correlation. 

summary(lm(X2 ~ Y2 + Z2)) # Negative significant correlation

## As we can see, in this particular case, conditioning on Z changes the sign 
##of the estimate for the Y2 variable, as in Simpson's paradox




dist_x <- function(N = 50, b_xz = 3, b_yx = 2, sd_x = 1, sd_y = 5, B = 2000,
                   z_min = 1, z_max = 5) {
  
  X_noZ <- rep(NA, B)
  X_yesZ <- rep(NA, B)
  
  pv_X_noZ <- rep(NA, B)
  pv_X_yesZ <- rep(NA, B)
  
  
  for(i in 1:B) {
    Z <- runif(N, z_min, z_max)
    X <- b_xz * Z + rnorm(N, 0, sd = sd_x)
    Y <- b_yx * X + rnorm(N, 0, sd = sd_y)
    
    m1 <- lm(Y ~ X)
    m2 <- lm(Y ~ X + Z)
    
    X_noZ[i] <- coefficients(m1)["X"]
    X_yesZ[i] <- coefficients(m2)["X"]
    
    pv_X_noZ[i] <- summary(m1)$coefficients["X", "Pr(>|t|)"]
    pv_X_yesZ[i] <- summary(m2)$coefficients["X", "Pr(>|t|)"]
    
    rm(Z, X, Y)
  }
  cat("\n Summary X without Z\n")
  print(summary(X_noZ))
  cat("\n s.d. estimate = ", sd(X_noZ))
  cat("\n\n Summary X with Z\n")
  print(summary(X_yesZ))
  cat("\n s.d. estimate = ", sd(X_yesZ), "\n")
  
  op <- par(mfrow = c(2, 2))
  hist(X_noZ, main = "X without Z in the model", xlab = "Estimate")
  abline(v = b_yx, lty = 2)
  hist(X_yesZ, main = "X with Z in the model", xlab = "Estimate")
  abline(v = b_yx, lty = 2)
  
  hist(pv_X_noZ, main = "X without Z in the model", xlab = "p-value")
  hist(pv_X_yesZ, main = "X with Z in the model", xlab = "p-value")
  
  par(op)
  
}


dist_x(b_yx = -2, N = 20)

## Increasing sd_x minimizes problem
dist_x(b_yx = -2, sd_x = 5)




drawdag(dag1)
drawdag(dag1, xlim = c(0.5, 3.5), ylim = c(-2, 1))


impliedConditionalIndependencies(dag1)

