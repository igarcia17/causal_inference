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
b_xz <- 3 # HDHRSHSTRHJTRSJTSRJTDJDJTDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
b_yz <- 2
sd_z <- 5

set.seed(11)

X <- runif(N, 1, 5)
Y <- runif(N, 2, 4)
Z <- b_xz * X + b_yz * Y + rnorm(N, 0, sd = sd_z)



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
X -> Z
Y -> Z
X -> Y
e_z -> Z
}")

coordinates(comm.effect.DAG_2) <- list(x = c(X = 1, Y = 3, Z = 2, 
                                             e_z = 1.75),
                                       y = c(X = 1, Y = 1, Z = 3, 
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
## of the estimate for the Y2 variable, as in Simpson's paradox


# if the estimate of both 



#_______________________________________________________________________________

# Collider with ancestor and descendant:

# Imagine we have now a variable (Z) that has a descendant (variable W),
# does conditioning on any of those two affect the relation between
# variable X and Y?

# ANCESTOR

comm.effect.DAG_3 <- dagitty("dag {
X -> Z
Y -> Z
W -> Z
e_z -> Z
}")

coordinates(comm.effect.DAG_3) <- list(x = c(X = 1, Y = 3, Z = 2, 
                                             e_z = 1.75, W = 2),
                                       y = c(X = 1, Y = 1, Z = 2, 
                                             e_z = 2, W = 3))
drawdag(comm.effect.DAG_3)

N <- 500 
b_xz <- 3 
b_yz <- 2
b_wz <- 2.5
sd_z <- 5

set.seed(11)

X3 <- runif(N, 1, 5)
Y3 <- runif(N, 2, 4)
W3 <- runif(N, 1.5, 3)
Z3 <- b_xz * X3 + b_yz * Y3 + W3 * b_wz + rnorm(N, 0, sd = sd_z)

summary(lm(X3 ~ W3)) # No significant correlation found
summary(lm(X3 ~ W3 + Z3)) # We find a negative correlation between X3 and W3
summary(lm(X3 ~ Y3 + Z3)) # We find a negative correlation between X3 and Y3

# Therefore, conditioning on Z3 modifies the independence between X3, Y3 and W3.

summary(lm(X3 ~ Y3)) # No significant correlation found
summary(lm(X3 ~ Y3 + W3)) # No significant correlation found

# But conditioning on W3 does not change the independence between X3 and Y3.




#DESCENDANT

comm.effect.DAG_4 <- dagitty("dag {
X -> Z
Y -> Z
Z -> W
e_z -> Z
}")

coordinates(comm.effect.DAG_4) <- list(x = c(X = 1, Y = 3, Z = 2, 
                                             e_z = 1.75, W = 2),
                                       y = c(X = 1, Y = 1, Z = 2, 
                                             e_z = 2, W = 3))
drawdag(comm.effect.DAG_4)


N <- 500 
b_xz <- 3 
b_yz <- 2
b_zw <- 2.5
sd_z <- 5
sd_w <- 2

set.seed(11)

X4 <- runif(N, 1, 5)
Y4 <- runif(N, 2, 4)
Z4 <- b_xz * X3 + b_yz * Y3 +  rnorm(N, 0, sd = sd_z)
W4 <- b_zw * Z4 + rnorm(N, 0, sd = sd_w)

summary(lm(X4 ~ Y4)) # No significant correlation found
summary(lm(X4 ~ Y4 + W4)) # We find a negative correlation between X3 and W3
summary(lm(X4 ~ Y4 + Z4)) # We find a negative correlation between X3 and W3

# In this case adjusting by Z4 and W4 (its descendant) resulted in a correlation
# between the variables X4 and Y4. 
