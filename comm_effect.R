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
b_xz <- 3 # The value of the relation between X and Z is 3
b_yz <- 2 # The value of the relation between Y and Z is 2
sd_z <- 5 # The standard deviation for Z will be 5

# We will not be repeating this information description as we already know what
# each variable represents.

set.seed(11)

X <- runif(N, 1, 10)
Y <- runif(N, 1.5, 12)
Z <- b_xz * X + b_yz * Y + rnorm(N, 0, sd = sd_z)


# Let's take a look at our data distribution:

reg_line_ZX <- lm(Z ~ X)
plot(Z ~ X, main = 'Var distribution X-Z')
abline(reg_line_ZX)


reg_line_ZY <- lm(Z ~ Y)
plot(Z ~ Y, main = 'Var distribution Y-Z')
abline(reg_line_ZY)

reg_line_XY <- lm(X ~ Y)
plot(X ~ Y, main = 'Var distribution X-Y')
abline(reg_line_XY)

# As we can see, there is, apparently, a positive correlation between Z and X and also
# Z and Y. On the other hand, there is no visible correlation between X and Y 
# (as should). We can check it by calculating the estimates and significance of
# each pair:

P.V.XZ <- summary(lm(X ~ Z))$coefficients['Z','Pr(>|t|)']
E.XZ <- summary(lm(X ~ Z))$coefficients['Z','Estimate']

cat('Your estimate value is', E.XZ, 'and your p value is', P.V.XZ)


P.V.YZ <- summary(lm(Y ~ Z))$coefficients['Z','Pr(>|t|)']
E.YZ <- summary(lm(Y ~ Z))$coefficients['Z','Estimate']

cat('Your estimate value is', E.YZ, 'and your p value is', P.V.YZ)



P.V.XY <- summary(lm(X ~ Y))$coefficients['Y','Pr(>|t|)']
E.XY <- summary(lm(X ~ Y))$coefficients['Y','Estimate']

cat('Your estimate value is', E.XY, 'and your p value is', P.V.XY)



P.V.XY.Z <- summary(lm(X ~ Y + Z))$coefficients['Y','Pr(>|t|)']
E.XY.Z <- summary(lm(X ~ Y + Z))$coefficients['Y','Estimate']

cat('Your estimate value is', E.XY.Z, 'and your p value is', P.V.XY.Z)


P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
assign('X_Y', P.V.XY, P.VAL.DICT)
assign('X_Y_Z', P.V.XY.Z, P.VAL.DICT)
assign('X_Z', P.V.XZ, P.VAL.DICT)
assign('Y_Z', P.V.YZ, P.VAL.DICT)

E.DICT <- new.env(hash = T, parent = emptyenv())
assign('X_Y', E.XY, E.DICT)
assign('X_Y_Z', E.XY.Z, E.DICT)
assign('X_Z', E.XZ, E.DICT)
assign('Y_Z', E.YZ, E.DICT)

for (i in ls(P.VAL.DICT)) {
  if (P.VAL.DICT[[i]] > 0.05)
    cat('The variables', i, 'do not show a correlation with a p value of',  
        P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
  else
    cat('The variables', i, 'show a correlation with a p value of',  
        P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
}


    
summary(lm(Z ~ X)) #we can see a strong correlation between Z and X
summary(lm(Z ~ Y)) #we can see a strong correlation between Z and Y
summary(lm(X ~ Y)) #we can't any correlation between X and Y
summary(lm(X ~ Y + Z))
# But the two independent variables (X and Y) when we adjust by Z 
# raise a correlation that wasn't supposed to be there.





#______________________________________________________________________________


# We can exemplify this with a more realistic dataset: 
# Here we have a set of variables: cigarettes smoked in a day, COVID19 virus
# load and % of lung capacity.


DAG.Lung <- dagitty("dag {
cigarettes_day -> lung_capacity
vir_load_COV19 -> lung_capacity
}")

coordinates(DAG.Lung) <- list(x = c(cigarettes_day = 1, vir_load_COV19 = 3, 
                                    lung_capacity = 2),
                                     y = c(cigarettes_day = 1, vir_load_COV19 = 1,
                                           lung_capacity = 3))
drawdag(DAG.Lung)

set.seed(11)

cigarettes_day <- floor(runif(N, min = 0, max = 45))
vir_load_COV19 <- floor(runif(N, min = 0, max = 40)) #Ct
lung_capacity <- 100 - 0.7 * cigarettes_day - 0.8 * vir_load_COV19 + 
  rnorm(N, 0, sd = 0.01) # % of lung capacity

data_frame(cigarettes_day, vir_load_COV19, lung_capacity)

summary(lm(cigarettes_day ~ vir_load_COV19)) # No correlation between nº of
# cigarettes smoked a day and the viral load of COV19 infection
summary(lm(cigarettes_day ~ vir_load_COV19 + lung_capacity)) # A correlation
# between these two variables arises when ajusting for the collider 
# (lung_capacity)


# Here we have the plots showing the distribution of the variables cigarettes_day
# and vir_load_COV19, where apparently there is no correlation.
reg_line_c_v <- lm(cigarettes_day ~ vir_load_COV19)
plot(cigarettes_day ~ vir_load_COV19)
abline(reg_line_c_v)

# A plot for comparison, cigarettes and lung capacity, a negative correlation:
reg_line_c_l <- lm(cigarettes_day ~ lung_capacity)
plot(cigarettes_day ~ lung_capacity)
abline(reg_line_c_l)



#_____________________________________________________________________________

# Let's take a look at a bit more complicated example: imagine we are studying 
# the performance of a professional runner. We have a few variables to take into 
# consideration, but we will focus on the following ones: leg length,
# metabolism (which affects weight, meaning it also affects performance), heart
# disease (MEASURED THROUGH ???????????????????????????????????????????),
# cholesterol (as a mesure of overall health); this last variable does not
# affect directly the performance of the runner (speed), but it might through
# its health. 
# Let's take a look at the DAG:

DAG.Runner <- dagitty("dag {
cholesterol -> heart
metabolism -> heart
metabolism -> speed
leg_lenght -> speed
}")

coordinates(DAG.Runner) <- list(x = c(cholesterol = 1, heart = 3, 
                                              metabolism = 1, speed = 3, 
                                              leg_lenght = 1),
                                       y = c(cholesterol = 5, heart = 4, 
                                             metabolism = 3, speed = 2, 
                                             leg_lenght = 1))
drawdag(DAG.Runner)


# To represent the effect of conditioning on a collider in this case we will
# be using data for every variable, but we will consider as if we could not 
# measure the metabolism in order to control the counfounder. 


# CORREGIR ALGUNOS DATOS PARA QUE TENGAN SENTIDO LOS VALORES DE HEART; AHORA
# DEPENDIENDO DE QUE MIDAMOS NO TIENE SENTIDO (POR EJEMPLO SI NOS BASAMOS
# EN LAS PULSACIONES EN REPOSO HAY UN COLEGA QUE TIENE 17, ESTÁ A PUNTITO
# DE PALMARLA ME DA A MI.)
set.seed(11)
N <- 500
leg_lenght <- runif(N, max = 49.75, min = 42.09)
metabolism <- runif(N, 1, 100) #corregir estos valores.
cholesterol <- runif(N2, 125, 200) #mg/dL
speed <- leg_lenght * 1.2 + metabolism * 1 + rnorm(N, mean = 0, sd = 5) #mph
heart <- cholesterol * 0.5 + metabolism * (-0.5) + rnorm(N2, mean = 0, sd = 0.1)

data.frame(leg_lenght, metabolism, cholesterol, speed, heart)


summary(lm(leg_lenght ~ metabolism)) # no correlation
summary(lm(leg_lenght ~ cholesterol)) # no correlation
summary(lm(leg_lenght ~ heart)) # no correlation
summary(lm(leg_lenght ~ heart + speed)) # we find a correlation between heart 
# and leg length that should not be there, there is no correlation between 
# how long your legs are and cardiac disease.

#Let's look at our data and see whether this correlation we find when adjusting
# by our collider (speed) is actually present:

reg_line_Runner_ll_h <- lm(leg_lenght ~ heart)
plot(leg_lenght ~ heart)
abline(reg_line_Runner_ll_h)

# The following plot shows the correlation that actually exists between
# the variables cholesterol and heart, as opposed to the previous ones.
reg_line_Runner_ch_h <- lm(cholesterol ~ heart)
plot(cholesterol ~ heart)
abline(reg_line_Runner_ch_h)

# SOLO LOS PONGO PARA VER MAS O MENOS LA DISTRIBUCIÓN DE LOS DATOS

reg_line_Runner_s_m <- lm(speed ~ metabolism)
plot(speed ~ metabolism)
abline(reg_line_Runner_s_m)

reg_line_Runner_h_m <- lm(heart ~ metabolism)
plot(heart ~ metabolism)
abline(reg_line_Runner_h_m)
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

# As we can see, there is a positive correlation between Z2 and X2; Z2 and Y2; 
# and X2 and Y2. In this case we are interested in the relationship between X2 
# and Y2; we can check it: 

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

