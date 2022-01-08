#_____________________COLLIDER AND SELECTION BIAS, CASES_______________________

# Import modules
library(dagitty)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}

set.seed(11)
# We are going to try two specific scenarios, one where we adjust for Z and
# another where we don't. The purpose of these two scenarios is to check whether
# adjusting for a collider will induce a bias in our estimates and how exactly
# this bias will change these estimates. 


#_________________________________SITUATION 1___________________________________

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

N <- 500 # Our sample size throughout the whole code will be 500
b_xz <- 3 # The value of the relation between X and Z is 3
b_yz <- 2 # The value of the relation between Y and Z is 2
sd_z <- 1 # The standard deviation for Z will be 5

# We will not be repeating this information description as we already know what
# each variable represents.

set.seed(11)

X <- runif(N, 1, 10)
Y <- runif(N, 1.5, 12)
Z <- b_xz * X + b_yz * Y + rnorm(N, 0, sd = sd_z)


# Let's take a look at our data distribution, for that we will create a 
# function that plots our variables with a regression line:
PLOT.REG <- function(reg_line, plot.var) {
  Regres.line <- lm(reg_line)
  plot(reg_line, main = plot.var)
  abline(Regres.line)
}

op <- par(mfrow= c(3, 1))
PLOT.REG(Z ~ X, 'Var distribution X-Z')
PLOT.REG(Z ~ Y, 'Var distribution Y-Z')
PLOT.REG(X ~ Y, 'Var distribution X-Y')
par(op)

# As we can see, there is, apparently, a positive correlation between Z and X and also
# Z and Y. On the other hand, there is no visible correlation between X and Y 
# (as should). We can check it by calculating the estimates and significance of
# each pair. For that we will create a function that accepts the following 
# arguments: variables to check, dependent variable (to get the estimate and
# p value), and the title of the variables; the last two ones should be text
# input in order to make it work.

SUM.2VAR <- function(variable1, v_dep1, title_var1, variable2, 
                     v_dep2, title_var2, ...){
  PV1 <- summary(lm(variable1))$coefficients[v_dep1, 'Pr(>|t|)']
  E1 <- summary(lm(variable1))$coefficients[v_dep1, 'Estimate']
  
  PV2 <- summary(lm(variable2))$coefficients[v_dep2, 'Pr(>|t|)']
  E2 <- summary(lm(variable2))$coefficients[v_dep2, 'Estimate']
  
  P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, PV1, P.VAL.DICT)
  assign(title_var2, PV2, P.VAL.DICT)
  
  E.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, E1, E.DICT)
  assign(title_var2, E2, E.DICT)
  
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0.05)
      cat('The variables', i, 'do not show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
} # SUMMARY for 2 variables

SUM.3VAR <- function(variable1, v_dep1, title_var1, variable2, 
                     v_dep2, title_var2, variable3, v_dep3, title_var3){
  PV1 <- summary(lm(variable1))$coefficients[v_dep1, 'Pr(>|t|)']
  E1 <- summary(lm(variable1))$coefficients[v_dep1, 'Estimate']
  
  PV2 <- summary(lm(variable2))$coefficients[v_dep2, 'Pr(>|t|)']
  E2 <- summary(lm(variable2))$coefficients[v_dep2, 'Estimate']
  
  PV3 <- summary(lm(variable3))$coefficients[v_dep3, 'Pr(>|t|)']
  E3 <- summary(lm(variable3))$coefficients[v_dep3, 'Estimate']
  
  P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, PV1, P.VAL.DICT)
  assign(title_var2, PV2, P.VAL.DICT)
  assign(title_var3, PV3, P.VAL.DICT)
  
  E.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, E1, E.DICT)
  assign(title_var2, E2, E.DICT)
  assign(title_var3, E3, E.DICT)
  
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0.05)
      cat('The variables', i, 'do not show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
}
                      # SUMMARY for 3 variables

SUM.4VAR <- function(variable1, v_dep1, title_var1, variable2, 
                     v_dep2, title_var2, variable3, v_dep3, title_var3, 
                     variable4, v_dep4, title_var4){
  PV1 <- summary(lm(variable1))$coefficients[v_dep1, 'Pr(>|t|)']
  E1 <- summary(lm(variable1))$coefficients[v_dep1, 'Estimate']
  
  PV2 <- summary(lm(variable2))$coefficients[v_dep2, 'Pr(>|t|)']
  E2 <- summary(lm(variable2))$coefficients[v_dep2, 'Estimate']
  
  PV3 <- summary(lm(variable3))$coefficients[v_dep3, 'Pr(>|t|)']
  E3 <- summary(lm(variable3))$coefficients[v_dep3, 'Estimate']
  
  PV4 <- summary(lm(variable4))$coefficients[v_dep4, 'Pr(>|t|)']
  E4 <- summary(lm(variable4))$coefficients[v_dep4, 'Estimate']
  
  P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, PV1, P.VAL.DICT)
  assign(title_var2, PV2, P.VAL.DICT)
  assign(title_var3, PV3, P.VAL.DICT)
  assign(title_var4, PV4, P.VAL.DICT)
  
  
  E.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, E1, E.DICT)
  assign(title_var2, E2, E.DICT)
  assign(title_var3, E3, E.DICT)
  assign(title_var4, E4, E.DICT)
  
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0.05)
      cat('The variables', i, 'do not show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
}
                      #SUMMARY for 4 variables


#Let's check our variables:
SUM.4VAR(X ~ Z, "Z", "X_Z", Y ~ Z, "Z", "Y_Z", X ~ Y, "Y", "X_Y", X ~ Y + Z, 
         "Y", "X_Y_Z")


# We can see a strong correlation between Z and X. There is also a positive
# relation between Z and Y, but none between X and Y. But the two independent 
# variables (X and Y) when we adjust by Z raise a correlation that wasn't 
# supposed to be there.

# When controling a collider, you can induce a bias in the estimate of the
# parent variables when there is, actually, none.




#_______________________REALISTIC EXAMPLE SITUATION 1___________________________


# We can exemplify this with a more realistic dataset: 
# Here we have a set of variables: cigarettes smoked in a day, COVID19 virus
# load and % of lung capacity. Both cigarette consumption and virus load should
# affect the lung capacity, but there should not be a relation between the 
# number of cigarettes consumed daily and the COV19 virus load.

#Let's take a look at the DAG:

DAG.Lung <- dagitty("dag {
cigarettes_day -> lung_capacity
vir_load_COV19 -> lung_capacity
}")

coordinates(DAG.Lung) <- list(x = c(cigarettes_day = 1, vir_load_COV19 = 3, 
                                    lung_capacity = 2),
                              y = c(cigarettes_day = 1, vir_load_COV19 = 1,
                                     lung_capacity = 3))
drawdag(DAG.Lung)


b_cd <- (-0.7) # Relation between cigarettes and lung capacity
b_vlc <- (-0.8) # Relation between virus load and lung capacity 

cigarettes_day <- floor(runif(N, min = 0, max = 45)) #number of cig a day
vir_load_COV19 <- floor(runif(N, min = 0, max = 40)) #Ct
lung_capacity <- 100 + b_cd * cigarettes_day + b_vlc * vir_load_COV19 + 
  rnorm(N, 0, sd = 0.01) # % of lung capacity


#Checking the data: 
data.frame(cigarettes_day, vir_load_COV19, lung_capacity)

#Let's see the relations between our variables.
SUM.4VAR(cigarettes_day ~ lung_capacity, "lung_capacity", "cig_lung",
         vir_load_COV19 ~ lung_capacity, "lung_capacity", "vir_lung",
         cigarettes_day ~ vir_load_COV19, "vir_load_COV19", "cig_vir",
          cigarettes_day ~ vir_load_COV19 + lung_capacity, "vir_load_COV19",
         "cig_vir_lung")


# No correlation between nº of cigarettes smoked a day and the viral load of 
# COV19 infection but a correlation between these two variables arises when 
# adjusting for the collider (lung_capacity)


# Here we have the plots showing the distribution of the variables cigarettes_day
# and vir_load_COV19, where apparently there is no correlation.
op <- par(mfrow= c(2, 1))
PLOT.REG(cigarettes_day ~ vir_load_COV19, "Cigarettes - Virus_load")

# A plot for comparison, cigarettes and lung capacity, a negative correlation:
PLOT.REG(cigarettes_day ~ lung_capacity, "Cigarettes - Lung_capacity")

par(op)


#___________________________MORE COMPLEX EXAMPLE________________________________

# Let's take a look at a bit more complicated example: imagine we are studying 
# the performance of a professional runner. We have a few variables to take into 
# consideration, but we will focus on the following ones: leg length,
# metabolism (which affects weight, meaning it also affects performance), heart
# disease (measured through the values of high-sensitivity cardiac troponin 
# (hs-cTn)), cholesterol (as a measure of overall health); this last variable 
# does not affect directly the performance of the runner (speed), but it might 
# through its health. 

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


b_ll_s <- 1.2
b_m_s <- 1
b_ch_h <- 0.05
b_m_h <- (-0.3)

leg_lenght <- runif(N, max = 49.75, min = 42.09) #cm
metabolism <- runif(N, 1, 20) # Theoretically this variable is unmeasured
cholesterol <- runif(N, 125, 200) #mg/dL
speed <- leg_lenght * b_ll_s + metabolism * b_m_s + rnorm(N, mean = 0, 
                                                          sd = 5) #mph
heart <- cholesterol * b_ch_h + metabolism * b_m_h + rnorm(N2, mean = 0, 
                                                           sd = 0.1) #ng/L

#Checking the data:
data.frame(leg_lenght, metabolism, cholesterol, speed, heart)

#Let's, one more time, check the relations:
SUM.4VAR(leg_lenght ~ metabolism, "metabolism", "len_metab", leg_lenght ~ cholesterol,
         "cholesterol", "len_chol", leg_lenght ~ heart, "heart", "len_heart", 
         leg_lenght ~ heart + speed, "heart", "len_heart_speed")

# There is no correlation between leg lenght and metabolism, leg lenght and 
# cholesterol and leg lenght and heart disease; but when we adjust for the 
# collider speed, we open a path between leg lenght and heart disease. We 
# do find  a correlation between heart when there is no correlation between 
# how long your legs are and cardiac disease.

#Let's look at our data and see whether this correlation we find when adjusting
# by our collider (speed) is actually present:
op <- par(mfrow= c(2, 1))
PLOT.REG(leg_lenght ~ heart, "LEG - HEART")

# The following plot shows the correlation that actually exists between
# the variables cholesterol and heart, as opposed to the previous ones.
PLOT.REG(cholesterol ~ heart, "CHOLESTEROL - HEART")

par(op)





#________________________________SITUATION 2___________________________________
# We could also have a case in which our X and Y variables have some sort of
# relation but the estimate (its correlation) changes from positive to negative
# when we condition on Z.

#Let's take a look at the DAG:

DAG_Situation2 <- dagitty("dag {
X -> Z
Y -> Z
X -> Y
e_z -> Z
}")

coordinates(DAG_Situation2) <- list(x = c(X = 1, Y = 3, Z = 2, 
                                             e_z = 1.75),
                                       y = c(X = 1, Y = 1, Z = 3, 
                                             e_z = 3))
drawdag(DAG_Situation2)

b_xy2 <- 1.7 
b_xz2 <- 1.5
b_yz2 <- 2.5

X2 <- runif(N, 1, 10)
Y2 <- X2 * b_xz2 + rnorm(N, mean = 0, sd = 1)
Z2 <- X2 * b_xz2 + Y2 * b_yz2 + rnorm(N, mean = 0, sd = 1)


# Again, let's take a look at our data distribution:
op <- par(mfrow= c(3, 1))

PLOT.REG(Z2 ~ X2, "Var distribution Z2 - X2")
PLOT.REG(Z2 ~ Y2, "Var distribution Z2 - Y2")
PLOT.REG(X2 ~ Y2, "Var distribution X2 - Y2")

par(op)

# As we can see, there is, apparently a positive correlation between Z2 and X2;  
# Z2 and Y2; and X2 and Y2. In this case we are interested in the relationship  
# between X2 and Y2; we can check it: 

SUM.2VAR(X2 ~ Y2, "Y2", "X2_Y2", X2 ~ Y2 + Z2, "Y2", "X2_Y2_Z2")

# The positive correlation we find between X2 and Y2 switches to
# a negative correlation when we adjust for Z2

## As we can see, in this particular case, conditioning on Z changes the sign 
## of the estimate for the Y2 variable, as in Simpson's paradox





#____________________________PRACTICAL EXAMPLE SITUATION 2______________________

# We want to study the effect of UV radiation (hours of exposure a day) and 
# INK4a on the mutation of gen p53. 

#First of all we create our DAG:

DAG.p53 <- dagitty("dag {
UV_radiation -> mutated_p53
INK4a -> mutated_p53
UV_radiation -> INK4a
e_z -> mutated_p53
}")

coordinates(DAG.p53) <- list(x = c(UV_radiation = 1, INK4a = 3, 
                                   mutated_p53 = 2,  e_z = 1.5),
                            y = c(UV_radiation = 1, INK4a = 1, 
                                  mutated_p53 = 3, e_z = 3))
drawdag(DAG.p53)


#X = UV.radiation (minutes of exposure/day), Y = INK4a, Z= mutatedp53
# Let's generate our data for this particular case:

b_UV_INK4a <- (0.2)
b_INK4a_mp53 <- 2
b_UV_mp53 <- 0.7


UV_radiation <- runif(N, 0, 9)
INK4a <- UV_radiation * b_UV_INK4a + rnorm(N, mean = 0, sd = 0.5)
mutatedp53 <- UV_radiation * b_UV_mp53 + INK4a * b_INK4a_mp53 + 
  rnorm(N, mean = 0, sd = 0.1)

# Again lets gather our data in a dataframe to check it:
data.frame(UV_radiation, INK4a, mutatedp53)


SUM.2VAR(UV.radiation ~ INK4a, "INK4a", "UV_INK4a", UV.radiation ~ INK4a
         + mutatedp53, "INK4a", "UV_INK4a_Mp53")

#___________________________ANCESTOR SITUATION__________________________________

# Imagine we have now a variable (Z) that has a descendant (variable W),
# does conditioning on any of those two affect the relation between
# variable X and Y?

# ANCESTOR

DAG_Situation3 <- dagitty("dag {
X -> Z
Y -> Z
W -> Z
e_z -> Z
}")

coordinates(DAG_Situation3) <- list(x = c(X = 1, Y = 3, Z = 2, 
                                             e_z = 1.75, W = 2),
                                       y = c(X = 1, Y = 1, Z = 2, 
                                             e_z = 2, W = 3))
drawdag(DAG_Situation3)


b_xz <- 3 
b_yz <- 2
b_wz <- 2.5
sd_z <- 1

set.seed(11)

X3 <- runif(N, 1, 5)
Y3 <- runif(N, 2, 4)
W3 <- runif(N, 1.5, 3)
Z3 <- b_xz * X3 + b_yz * Y3 + W3 * b_wz + rnorm(N, 0, sd = sd_z)

SUM.3VAR(X3 ~ W3, "W3", "X3_W3", X3 ~ W3 + Z3, "W3", "X3_W3_Z3", X3 ~ Y3 + Z3, 
         "Y3", "X3_Y3_Z3")



summary(lm(X3 ~ W3)) # No significant correlation found
summary(lm(X3 ~ W3 + Z3)) # We find a negative correlation between X3 and W3
summary(lm(X3 ~ Y3 + Z3)) # We find a negative correlation between X3 and Y3

# Therefore, conditioning on Z3 modifies the independence between X3, Y3 and W3.

summary(lm(X3 ~ Y3)) # No significant correlation found
summary(lm(X3 ~ Y3 + W3)) # No significant correlation found

# But conditioning on W3 does not change the independence between X3 and Y3.

# Let's take a look at a more realistic example: we want to study the effect
# of cortisol and cholesterol on diabetes, but we also have to take into 
# consideration the variable "sugar consumption in gr". Here is the DAG:

DAG.Diabetes.Ancestor <- dagitty("dag {
cortisol -> diabetes
cholesterol -> diabetes
sugar_consumpt -> diabetes
e_z -> diabetes
}")

coordinates(DAG.Diabetes.Ancestor) <- list(x = c(cortisol = 1, cholesterol = 3,
                                           diabetes = 2, e_z = 1.5, 
                                            sugar_consumpt = 2),
                                       y = c(cortisol = 1, cholesterol = 1, 
                                             diabetes = 2, e_z = 2, 
                                             sugar_consumpt = 3))
drawdag(DAG.Diabetes.Ancestor)

#X = cortisol, Y = colesterol, Z = diabetes, W = sugar consumption in gr

b_co_d <- 0.8
b_ch_d <- 2
b_sc_d <- 2.5
sd_z <- 1
sd_w <- 0.5

set.seed(11)

cortisol <- runif(N, 1, 5)
cholesterol <- runif(N, 2, 4)
diabetes <- cortisol * b_co_d + cholesterol * b_ch_d + sug_consumption * 
  b_sc_d+ rnorm(N, 0, sd = sd_z)
sug_consumption <- runif(N, 0, 150)



#___________________________DESCENDANT SITUATION________________________________

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


SUM.3VAR(X4 ~ Y4, "Y4", "X4_Y4", X4 ~ Y4 + W4, "Y4", "X4_Y4_W4", X4 ~ Y4 + Z4,
         "Y4", "X4_Y4_Z4")


summary(lm(X4 ~ Y4)) # No significant correlation found
summary(lm(X4 ~ Y4 + W4)) # We find a negative correlation between X3 and W3
summary(lm(X4 ~ Y4 + Z4)) # We find a negative correlation between X3 and W3

# In this case adjusting by Z4 and W4 (its descendant) resulted in a correlation
# between the variables X4 and Y4. 





#Practical example: we want to see how does diabetes affect heart disease
# for that we will also measure cortisol and cholesterol, which affect 
# diabetes. Let's see what happens when we condition the collider or its
# descendant in this particular case. First we create the DAG:

DAG.Diabetes.Descendant <- dagitty("dag {
cortisol -> diabetes
cholesterol -> diabetes
diabetes -> heart_disease
e_z -> diabetes
}")

coordinates(DAG.Diabetes.Descendant) <- list(x = c(cortisol = 1, 
                                        cholesterol = 3, diabetes = 2, 
                                        e_z = 1.5, heart_disease = 2),
                                       y = c(cortisol = 1, cholesterol = 1, 
                                             diabetes = 2, e_z = 2, 
                                             heart_disease = 3))
drawdag(DAG.Diabetes.Descendant)

#X = cortisol, Y = colesterol, Z = diabetes, W = heart disease


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


