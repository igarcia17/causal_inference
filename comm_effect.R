# Collider and selection bias, cases:

# Import modules
library(dagitty)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}



# We are going to try two specific scenarios, one where we adjust for Z and
# another where we don't. The purpose of these two scenarios is to check whether
# adjusting for a collider will induce a bias in our estimates and how exactly
# this bias will change these estimates. 


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

PLOT.REG <- function(reg_line, plot.var) {
  Regres.line <- lm(reg_line)
  plot(reg_line, main = plot.var)
  abline(Regres.line)
}

PLOT.REG(Z ~ X, 'Var distribution X-Z')
PLOT.REG(Z ~ Y, 'Var distribution Y-Z')
PLOT.REG(X ~ Y, 'Var distribution X-Y')


# As we can see, there is, apparently, a positive correlation between Z and X and also
# Z and Y. On the other hand, there is no visible correlation between X and Y 
# (as should). We can check it by calculating the estimates and significance of
# each pair:

P.V.XZ <- summary(lm(X ~ Z))$coefficients['Z','Pr(>|t|)']
E.XZ <- summary(lm(X ~ Z))$coefficients['Z','Estimate']

P.V.YZ <- summary(lm(Y ~ Z))$coefficients['Z','Pr(>|t|)']
E.YZ <- summary(lm(Y ~ Z))$coefficients['Z','Estimate']

P.V.XY <- summary(lm(X ~ Y))$coefficients['Y','Pr(>|t|)']
E.XY <- summary(lm(X ~ Y))$coefficients['Y','Estimate']

P.V.XY.Z <- summary(lm(X ~ Y + Z))$coefficients['Y','Pr(>|t|)']
E.XY.Z <- summary(lm(X ~ Y + Z))$coefficients['Y','Estimate']

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

#we can see a strong correlation between Z and X we can see a strong correlation
# between Z and Y we can't any correlation between X and Y

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

b_cd <- (-0.7) # Relation between lung capacity and cigarettes
b_vlc <- (-0.8) # Relation between lung capacity and virus load

cigarettes_day <- floor(runif(N, min = 0, max = 45))
vir_load_COV19 <- floor(runif(N, min = 0, max = 40)) #Ct
lung_capacity <- 100 - b_cd * cigarettes_day - b_vlc * vir_load_COV19 + 
  rnorm(N, 0, sd = 0.01) # % of lung capacity

data.frame(cigarettes_day, vir_load_COV19, lung_capacity)

P.V.C.L <- summary(lm(cigarettes_day ~ lung_capacity)
                   )$coefficients['lung_capacity','Pr(>|t|)']
E.C.L <- summary(lm(cigarettes_day ~ lung_capacity)
                 )$coefficients['lung_capacity','Estimate']

P.V.V.L <- summary(lm(vir_load_COV19 ~ lung_capacity)
                   )$coefficients['lung_capacity','Pr(>|t|)']
E.V.L <- summary(lm(vir_load_COV19 ~ lung_capacity)
                 )$coefficients['lung_capacity','Estimate']

P.V.C.V <- summary(lm(cigarettes_day ~ vir_load_COV19)
                  )$coefficients['vir_load_COV19','Pr(>|t|)']
E.C.V <- summary(lm(cigarettes_day ~ vir_load_COV19)
                 )$coefficients['vir_load_COV19','Estimate']



P.VAL.DICT.Lung <- new.env(hash = T, parent = emptyenv())
assign('cig_lung', P.V.C.L, P.VAL.DICT.Lung)
assign('vir_lung', P.V.V.L, P.VAL.DICT.Lung)
assign('cig_vir', P.V.C.V, P.VAL.DICT.Lung)

E.DICT.Lung <- new.env(hash = T, parent = emptyenv())
assign('cig_lung', E.C.L, E.DICT.Lung)
assign('vir_lung', E.V.L, E.DICT.Lung)
assign('cig_vir', E.C.V, E.DICT.Lung)



for (i in ls(P.VAL.DICT.Lung)) {
  if (P.VAL.DICT.Lung[[i]] > 0.05)
    cat('The variables', i, 'do not show a correlation with a p value of',  
        P.VAL.DICT.Lung[[i]], 'and an estimate of', E.DICT.Lung[[i]], '\n')
  else
    cat('The variables', i, 'show a correlation with a p value of',  
        P.VAL.DICT.Lung[[i]], 'and an estimate of', E.DICT.Lung[[i]], '\n')
}


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
# disease (measured through the values of high-sensitivity cardiac troponin (hs-cTn)),
# cholesterol (as a measure of overall health); this last variable does not
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


set.seed(11)

b_ll_s <- 1.2
b_m_s <- 1
b_ch_h <- 0.05
b_m_h <- (-0.3)

leg_lenght <- runif(N, max = 49.75, min = 42.09)
metabolism <- runif(N, 1, 20) # Theoretically this variable is unmeasured
cholesterol <- runif(N, 125, 200) #mg/dL
speed <- leg_lenght * b_ll_s + metabolism * b_m_s + rnorm(N, mean = 0, sd = 5) #mph
heart <- cholesterol * b_ch_h + metabolism * b_m_h + rnorm(N2, mean = 0, sd = 0.1) #ng/L

data.frame(leg_lenght, metabolism, cholesterol, speed, heart)


P.V.L.M <- summary(lm(leg_lenght ~ metabolism))$coefficients['metabolism','Pr(>|t|)']
E.L.M <- summary(lm(leg_lenght ~ metabolism))$coefficients['metabolism','Estimate']

P.V.L.Ch <- summary(lm(leg_lenght ~ cholesterol))$coefficients['cholesterol','Pr(>|t|)']
E.L.Ch <- summary(lm(leg_lenght ~ cholesterol))$coefficients['cholesterol','Estimate']

P.V.L.H <- summary(lm(leg_lenght ~ heart))$coefficients['heart','Pr(>|t|)']
E.L.H <- summary(lm(leg_lenght ~ heart))$coefficients['heart','Estimate']

P.V.L.H_Sp <- summary(lm(leg_lenght ~ heart + speed))$coefficients['heart','Pr(>|t|)']
E.L.H_Sp <- summary(lm(leg_lenght ~ heart + speed))$coefficients['heart','Estimate']


P.VAL.DICT.Speed <- new.env(hash = T, parent = emptyenv())
assign('len_metab', P.V.L.M, P.VAL.DICT.Speed)
assign('len_chol', P.V.L.Ch, P.VAL.DICT.Speed)
assign('len_heart', P.V.L.H, P.VAL.DICT.Speed)
assign('len_heart+Speed', P.V.L.H_Sp, P.VAL.DICT.Speed)



E.DICT.Speed <- new.env(hash = T, parent = emptyenv())
assign('len_metab', E.L.M, E.DICT.Speed)
assign('len_chol', E.L.Ch, E.DICT.Speed)
assign('len_heart', E.L.H, E.DICT.Speed)
assign('len_heart+Speed', E.L.H_Sp, E.DICT.Speed)



for (i in ls(P.VAL.DICT.Speed)) {
  if (P.VAL.DICT.Speed[[i]] > 0.05)
    cat('The variables', i, 'do not show a correlation with a p value of',  
        P.VAL.DICT.Speed[[i]], 'and an estimate of', E.DICT.Speed[[i]], '\n')
  else
    cat('The variables', i, 'show a correlation with a p value of',  
        P.VAL.DICT.Speed[[i]], 'and an estimate of', E.DICT.Speed[[i]], '\n')
}

summary(lm(leg_lenght ~ metabolism)) # no correlation
summary(lm(leg_lenght ~ cholesterol)) # no correlation
summary(lm(leg_lenght ~ heart)) # no correlation
summary(lm(leg_lenght ~ heart + speed)) # we find a correlation between heart 
# and leg length that should not be there, there is no correlation between 
# how long your legs are and cardiac disease.

#Let's look at our data and see whether this correlation we find when adjusting
# by our collider (speed) is actually present:
PLOT.REG(leg_lenght ~ heart, "LEG - HEART")

# The following plot shows the correlation that actually exists between
# the variables cholesterol and heart, as opposed to the previous ones.
PLOT.REG(cholesterol ~ heart, "CHOLESTEROL - HEART")


# SOLO LOS PONGO PARA VER MAS O MENOS LA DISTRIBUCIÓN DE LOS DATOS
PLOT.REG(speed ~ metabolism, "SPEED - METABOLISM")
PLOT.REG(heart ~ metabolism, "HEART - METABOLISM")


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

PLOT.REG(Z2 ~ X2, "Var distribution Z2 - X2")
PLOT.REG(Z2 ~ Y2, "Var distribution Z2 - Y2")
PLOT.REG(X2 ~ Y2, "Var distribution X2 - Y2")


P.V.XY2 <- summary(lm(X2 ~ Y2))$coefficients['Y2','Pr(>|t|)']
E.L.XY2 <- summary(lm(X2 ~ Y2))$coefficients['Y2','Estimate']

P.V.XY2_Z2 <- summary(lm(X2 ~ Y2 + Z2))$coefficients['Y2','Pr(>|t|)']
E.L.XY2_Z2 <- summary(lm(X2 ~ Y2 + Z2))$coefficients['Y2','Estimate']


P.VAL.DICT.Scenario2 <- new.env(hash = T, parent = emptyenv())
assign('X2_Y2', P.V.XY2, P.VAL.DICT.Scenario2)
assign('X2_Y2+Z2', P.V.XY2_Z2, P.VAL.DICT.Scenario2)


E.DICT.Scenario2 <- new.env(hash = T, parent = emptyenv())
assign('X2_Y2', E.L.XY2, E.DICT.Scenario2)
assign('X2_Y2+Z2', E.L.XY2_Z2, E.DICT.Scenario2)


# AQUI ALGO NO VA BIEN PORQUE ME SACA UN P VALUE DE 0 EN X2_Y2 PERO COMO AHORA
# MISMO NO LO VEO LO DEJO PARA MAS ADELANTE.
for (i in ls(P.VAL.DICT.Scenario2)) {
  if (P.VAL.DICT.Scenario2[[i]] > 0.05)
    cat('The variables', i, 'do not show a correlation with a p value of',  
        P.VAL.DICT.Scenario2[[i]], 'and an estimate of', E.DICT.Scenario2[[i]], '\n')
  else
    cat('The variables', i, 'show a correlation with a p value of',  
        P.VAL.DICT.Scenario2[[i]], 'and an estimate of', E.DICT.Scenario2[[i]], '\n')
}
# As we can see, there is a positive correlation between Z2 and X2; Z2 and Y2; 
# and X2 and Y2. In this case we are interested in the relationship between X2 
# and Y2; we can check it: 

summary(lm(X2 ~ Y2)) # Positive significant correlation

# But if we adjust for Z2 the estimate of the X2-Y2 correlation switches to
# a negative correlation. 

summary(lm(X2 ~ Y2 + Z2)) # Negative significant correlation

## As we can see, in this particular case, conditioning on Z changes the sign 
## of the estimate for the Y2 variable, as in Simpson's paradox


# Practical example:




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

