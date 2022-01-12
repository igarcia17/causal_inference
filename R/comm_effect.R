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
}")

coordinates(comm.effect.DAG) <- list(x = c(X = 1, Y = 3, Z = 2),
                                     y = c(X = 1, Y = 1, Z = 3))
drawdag(comm.effect.DAG)


# Let us generate the data according to the DAG:

set.seed(11)

N <- 500 # Our sample size throughout the whole code will be 500
b_xz <- 3 # The value of the relation between X and Z is 3
b_yz <- 2 # The value of the relation between Y and Z is 2
sd_z <- 1 # The standard deviation for Z will be 5

# We will try not to repeat this information description as we already know what
# each variable represents.

X <- runif(N, 1, 10)
Y <- runif(N, 1.5, 12)
Z <- b_xz * X + b_yz * Y + rnorm(N, 0, sd = sd_z)


# Let's take a look at our data distribution, for that we will use a 
# function that plots our variables with a regression line:
PLOT.REG <- function(reg_line, plot.var = '') {
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
      cat('The variables', i, 'do not show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show an association with a p value of',  
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
      cat('The variables', i, 'do not show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
} # SUMMARY for 3 variables

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
      cat('The variables', i, 'do not show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
} #SUMMARY for 4 variables


#Let's check our variables:
SUM.4VAR(X ~ Z, "Z", "X and Z", Y ~ Z, "Z", "Y and Z", X ~ Y, "Y", "X and Y", X ~ Y + Z, 
         "Y", "X and Y (when Z)")


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


# Here we have the plots showing the distribution of the variables cigarettes_day
# and vir_load_COV19, where apparently there is no correlation.
op <- par(mfrow= c(2, 1))
PLOT.REG(cigarettes_day ~ vir_load_COV19, "Cigarettes ~ Virus_load")

# A plot for comparison, cigarettes and lung capacity, a negative correlation:
PLOT.REG(cigarettes_day ~ lung_capacity, "Cigarettes ~ Lung_capacity")

par(op)



#Let's see the relations between our variables.
SUM.4VAR(cigarettes_day ~ lung_capacity, "lung_capacity", "cigarettes and lung capacity",
         vir_load_COV19 ~ lung_capacity, "lung_capacity", "vir_load_COV19 and lung capacity",
         cigarettes_day ~ vir_load_COV19, "vir_load_COV19", "cigarettes and vir_load_COV19",
          cigarettes_day ~ vir_load_COV19 + lung_capacity, "vir_load_COV19",
         "cigarettes and vir_load_COV19 (when lung capacity)")


# No correlation between nº of cigarettes smoked a day and the viral load of 
# COV19 infection but a correlation between these two variables arises when 
# adjusting for the collider (lung_capacity)



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
heart <- cholesterol * b_ch_h + metabolism * b_m_h + rnorm(N, mean = 0, 
                                                           sd = 0.1) #ng/L

#Checking the data:
data.frame(leg_lenght, metabolism, cholesterol, speed, heart)

PLOT.REG(leg_lenght ~ heart, "LEG_LEN ~ HEART")

#Let's, one more time, check the relations:
SUM.4VAR(leg_lenght ~ metabolism, "metabolism", "leg_len and metabolism", leg_lenght ~ cholesterol,
         "cholesterol", "leg_len and cholesterol", leg_lenght ~ heart, "heart", "leg_len and heart", 
         leg_lenght ~ heart + speed, "heart", "leg_len and heart (when speed)")

# There is no correlation between leg lenght and metabolism, leg lenght and 
# cholesterol and leg lenght and heart disease; but when we adjust for the 
# collider speed, we open a path between leg lenght and heart disease. We 
# do find  a correlation between heart when there is no correlation between 
# how long your legs are and cardiac disease.

#Let's look at our data and see whether this correlation we find when adjusting
# by our collider (speed) is actually present:
op <- par(mfrow= c(2, 1))
PLOT.REG(leg_lenght ~ heart, "LEG ~ HEART")

# The following plot shows the correlation that actually exists between
# the variables cholesterol and heart, as opposed to the previous ones.
PLOT.REG(cholesterol ~ heart, "CHOLESTEROL ~ HEART")

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
}")

coordinates(DAG_Situation2) <- list(x = c(X = 1, Y = 3, Z = 2),
                                       y = c(X = 1, Y = 1, Z = 3))
drawdag(DAG_Situation2)

b_xy2 <- 1.7 
b_xz2 <- 1.5
b_yz2 <- 2.5

X2 <- runif(N, 1, 10)
Y2 <- X2 * b_xz2 + rnorm(N, mean = 0, sd = 1)
Z2 <- X2 * b_xz2 + Y2 * b_yz2 + rnorm(N, mean = 0, sd = 1)


# Again, let's take a look at our data distribution:
op <- par(mfrow= c(3, 1))

PLOT.REG(Z2 ~ X2, "Var distribution Z2 ~ X2")
PLOT.REG(Z2 ~ Y2, "Var distribution Z2 ~ Y2")
PLOT.REG(X2 ~ Y2, "Var distribution X2 ~ Y2")

par(op)

# As we can see, there is, apparently a positive correlation between Z2 and X2;  
# Z2 and Y2; and X2 and Y2. In this case we are interested in the relationship  
# between X2 and Y2; we can check it: 

SUM.2VAR(X2 ~ Y2, "Y2", "X2 and Y2", X2 ~ Y2 + Z2, "Y2", "X2 and Y2 (when Z2)"){
if (P.VAL.DICT[[i]] > 0)
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0)
      cat('The variables', i, 'have a positive correlation',
          'with an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'have a negative correlation',
          'with an estimate of', E.DICT[[i]], '\n')}}


# The positive correlation we find between X2 and Y2 switches to
# a negative correlation when we adjust for Z2

# As we can see, in this particular case, conditioning on Z changes the sign 
# of the estimate for the Y2 variable, as in Simpson's paradox




#____________________________PRACTICAL EXAMPLE SITUATION 2______________________

# We want to study the effect of UV radiation (hours of exposure a day) and 
# INK4a on the mutation of gen p53. 

#First of all we create our DAG:

DAG.p53 <- dagitty("dag {
UV_radiation -> mutated_p53
INK4a -> mutated_p53
UV_radiation -> INK4a
}")

coordinates(DAG.p53) <- list(x = c(UV_radiation = 1, INK4a = 3, 
                                   mutated_p53 = 2),
                            y = c(UV_radiation = 1, INK4a = 1, 
                                  mutated_p53 = 3))
drawdag(DAG.p53)


# Let's generate our data for this particular case:

b_UV_INK4a <- 0.2
b_INK4a_mp53 <- 2
b_UV_mp53 <- 0.7


UV_radiation <- runif(N, 0, 9)
INK4a <- UV_radiation * b_UV_INK4a + rnorm(N, mean = 0, sd = 0.5)
mutatedp53 <- UV_radiation * b_UV_mp53 + INK4a * b_INK4a_mp53 + 
  rnorm(N, mean = 0, sd = 0.1)

# Again lets gather our data in a dataframe to check it:
data.frame(UV_radiation, INK4a, mutatedp53)

PLOT.REG(UV_radiation ~ INK4a, 'UV ~ INK4a')


SUM.2VAR(UV_radiation ~ INK4a, "INK4a", "UV and INK4a", UV_radiation ~ INK4a
         + mutatedp53, "INK4a", "UV and INK4a (when Mp53)")

# As shown in the Summary function, the correlation between the variables UV 
# radiation is maintained after adjusting for the collider but the relation
# that was once positive becomes negative when adjusting for Mp53 (collider)




#___________________________ANCESTOR SITUATION__________________________________

# Imagine we have now a variable Z (a collider) that has an ancestor 
# (variable W), does conditioning on any of those two affect the relation 
# between variable X and Y?

# Let's draw our DAG: 

DAG_Ancestor <- dagitty("dag {
X -> Z
Y -> Z
W -> Z
}")

coordinates(DAG_Ancestor) <- list(x = c(X = 1, Y = 3, Z = 2, W = 2),
                                       y = c(X = 1, Y = 1, Z = 2, W = 3))
drawdag(DAG_Ancestor)


b_xz <- 3 
b_yz <- 2
b_wz <- 2.5
sd_z <- 1


X3 <- runif(N, 1, 5)
Y3 <- runif(N, 2, 4)
W3 <- runif(N, 1.5, 3)
Z3 <- b_xz * X3 + b_yz * Y3 + W3 * b_wz + rnorm(N, 0, sd = sd_z)

SUM.4VAR(X3 ~ W3, "W3", "X3 and W3", X3 ~ W3 + Z3, "W3", "X3 and W3 (when Z3)", X3 ~ Y3 + Z3, 
         "Y3", "X3 and Y3 (when Z3)", X3 ~ Y3, "Y3", "X3 and Y3")

# There is no significant correlation between X3 and W3, but when conditioned
# for Z3 a negative correlation arises. Same thing hapopens between the 
# variables X3 and Y3.
# Therefore, conditioning on Z3 modifies the independence between X3, Y3 and W3

SUM.2VAR(X3 ~ Y3, "Y3", "X3 and Y3", X3 ~ Y3 + W3, "Y3", "X3 and Y3 (when W3)")

# No significant correlation found between X3 and Y3, nor can we find a 
# correlation when adjusting for the ancestor (W3). 
# So, conditioning on W3 does not change the independence between X3 and Y3.




#__________________________PRACTICAL EXAMPLE ANCESTOR_________________________

# Let's take a look at a more realistic example: we want to study the effect
# of cortisol and cholesterol on diabetes, but we also have to take into 
# consideration the variable "sugar consumption in gr". Here is the DAG:

DAG.Diabetes.Ancestor <- dagitty("dag {
cortisol -> diabetes
cholesterol -> diabetes
sugar_consumpt -> diabetes
}")

coordinates(DAG.Diabetes.Ancestor) <- list(x = c(cortisol = 1, cholesterol = 3,
                                           diabetes = 2, 
                                            sugar_consumpt = 2),
                                       y = c(cortisol = 1, cholesterol = 1, 
                                             diabetes = 2, 
                                             sugar_consumpt = 3))
drawdag(DAG.Diabetes.Ancestor)

#Once more, we create our data:

b_co_d <- 0.9
b_ch_d <- 0.95
b_sc_d <- 1.2
sd_d <- 1


cortisol <- runif(N, 1, 50)
cholesterol <- runif(N, 10, 250)
diabetes <- cortisol * b_co_d + cholesterol * b_ch_d + sug_consumption * 
  b_sc_d+ rnorm(N, 0, sd = sd_d)
sug_consumption <- runif(N, 0, 150)

PLOT.REG(cortisol ~ cholesterol, 'cortisol ~ cholesterol')

SUM.2VAR(cortisol ~ cholesterol, "cholesterol", "cortisol and cholesterol", cortisol 
         ~ cholesterol + sug_consumption, "cholesterol", "colesterol and cortisol (when sugar)")

# Again we find no correlation between cortisol and cholesterol, even if we
# adjust for the ancestor (sugar_consumption).

SUM.2VAR(cortisol ~ cholesterol, "cholesterol", "cortisol and cholesterol", cortisol 
         ~ cholesterol + diabetes, "cholesterol", "colesterol and cortisol (when diabetes)")

# But we do find a correlation between those two variables when adjusting
# for diabetes (the collider)




#___________________________DESCENDANT SITUATION________________________________

# Let's take a look at our second possible scenario with the ancestor and 
# descendant situations.
# We now have a variable Z (a collider) that has a descendant (variable W),
# does conditioning on any of those two affect the relation between
# variable X and Y?

# Starting with the DAG:

DAG.Descendant <- dagitty("dag {
X -> Z
Y -> Z
Z -> W
}")

coordinates(DAG.Descendant) <- list(x = c(X = 1, Y = 3, Z = 2, W = 2),
                                       y = c(X = 1, Y = 1, Z = 2, W = 3))
drawdag(DAG.Descendant)


b_xz <- 1.1
b_yz <- 0.7
b_zw <- 2
sd_z <- 0.5
sd_w <- 0.5

X4 <- runif(N, 1, 10)
Y4 <- runif(N, 2, 20)
Z4 <- b_xz * X4 + b_yz * Y4 +  rnorm(N, 0, sd = sd_z)
W4 <- b_zw * Z4 + rnorm(N, 0, sd = sd_w)


SUM.3VAR(X4 ~ Y4, "Y4", "X4 and Y4", X4 ~ Y4 + W4, "Y4", "X4 and Y4 (when W4)", X4 ~ Y4 + Z4,
         "Y4", "X4 and Y4 (when Z4)")

# We find no correlation between X4 and Y4, but when we condition on the 
# collider (Z) or its descendant (W4) we find a significantly negative 
# correlation between the two variables. 
# In this case adjusting by Z4 and W4 (its descendant) resulted in a correlation
# between the variables X4 and Y4. 




#________________________PRACTICAL EXAMPLE DESCENDANT___________________________

#Practical example: we want to see how does diabetes affect heart disease
# for that we will also measure cortisol and cholesterol, which affect 
# diabetes. Let's see what happens when we condition the collider or its
# descendant in this particular case. First we create the DAG:

DAG.Diabetes.Descendant <- dagitty("dag {
cortisol -> diabetes
cholesterol -> diabetes
diabetes -> heart_disease
}")

coordinates(DAG.Diabetes.Descendant) <- list(x = c(cortisol = 1, 
                                        cholesterol = 3, diabetes = 2,
                                        heart_disease = 2),
                                       y = c(cortisol = 1, cholesterol = 1, 
                                             diabetes = 2,  
                                             heart_disease = 3))
drawdag(DAG.Diabetes.Descendant)


# As the example is very simmilar to the previous one, we well be using the
# same values for the data, excepting the heart disease variable.
 
b_co_d_d <- 1.5
b_ch_d_d <- 0.95
b_d_hd <- 2.5
sd_d_d <- 5
sd_hd <- 2


cortisol_desc <- runif(N, 1, 50)
cholesterol_desc <- runif(N, 10, 250)
diabetes_desc <- cortisol_desc * b_co_d_d + cholesterol_desc * b_ch_d_d + 
  rnorm(N, 0, sd = sd_d_d)
heart_disease <- diabetes_desc * b_d_hd + rnorm(N, 0, sd = sd_hd)


data.frame(cortisol_desc, cholesterol_desc, diabetes_desc, heart_disease)

SUM.3VAR(cortisol_desc ~ cholesterol_desc, "cholesterol_desc", "cortisol and cholesterol", 
         cortisol_desc ~ cholesterol_desc + diabetes_desc, "cholesterol_desc",
         "cortisol and cholesterol (when diabetes)", cortisol_desc ~ cholesterol_desc + heart_disease,
         "cholesterol_desc", "cortisol and cholesterol (when heart disease)")

# Just like in our previous simple example, in this case there should not be
# any relation between cortisol and cholesterol (as we did intend with our data)
# and when we see our summary we can confirm there is no relation. But when we 
# condition on our collider or its descendant (heart disease) a correlation 
# between the variables cortisol and cholesterol arises. 

op <- par(mfrow= c(2, 1))

PLOT.REG(cortisol_desc ~ cholesterol_desc, 'cortisol ~ cholesterol')

# To better visualize this non existent correlation between our two variables
# we can check the plot for the regression of both of them and it's pretty clear
# that there is no apparent positive or negative correlation.

PLOT.REG(cortisol_desc ~ heart_disease, 'cortisol ~ heart_disease')

# To reinforce this data, we can also check and compare the plot for the 
# variables cortisol and heart disease, where we see a positive correlation.
