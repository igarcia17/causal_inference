#Confounders, common cause and mediators

#import modules
#library(tidyverse)
library(dagitty)
library(car)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
} 

#Imagine you want to study what variables influence the Z variable, having as well
#X and Y as covariates
#You may want to consider the DAG:

scenario1.DAG <- dagitty("dag {
X -> Y
X -> Z
e_y -> Y
e_z -> Z
}")

coordinates(scenario1.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                    y = c(Y = 3, X = 1, Z = 3, e_y = 2.75, e_z = 2.75))

drawdag(scenario1.DAG)

#As well as the DAG:

scenario2.DAG <- dagitty("dag {
X -> Y
X -> Z
Y -> Z
e_y -> Y
e_z -> Z
}")

coordinates(scenario2.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                  y = c(Y = 3, X = 1, Z = 3, e_y = 2.75, e_z = 2.75))
drawdag(scenario2.DAG)
#In both cases X is common cause of Y and Z. In the first case, Y has a causal effect
#on Z, while in the second case they are independent.

#How can you know which is the case? How do you know if you want to consider X, Y or both?
#Let's analyze all the possibilities!

#Let's create a function that creates the different datasets that we need.
create.dataset <- function(b_yz, N = 500, b_xy = 3, b_xz = 3,
                           e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}


#This function will be handy later
create.datasetv2 <- function(b_ax, b_yz=0, N = 500, b_xy = 3, b_xz = 3,
                             e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(A = runif(N, 1, 100) + rnorm(N))
  name_df$X <- name_df$A * b_ax + rnorm(N, sd = e_x)
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}

#We create a function that checks which is the scenario. We assume that X is in any
#case a common cause of both Y and Z

#The argument must be a data frame, which column names are Z (the variable of study)
#Y and X, the common cause.

Y_check <- function (dataset, conflevel = 0.01) {
  
  model_with_Y <- lm(Z~X+Y, data = dataset)
  p.v.X <-(summary(model_with_Y)$coefficients['X','Pr(>|t|)'])
  p.v.Y <- (summary(model_with_Y)$coefficients['Y', 'Pr(>|t|)'])
  
  if ((p.v.X <= conflevel)&(p.v.Y > conflevel))
  {cat("The variable of analysis is not influenced by Y\n")
    cat('See plot\n')
    drawdag(scenario1.DAG)
    return(invisible(1))
    }
  
  if ((p.v.X <= conflevel)&(p.v.Y <= conflevel))
  {cat("The variable of analysis is influenced by both X and Y\n")
    cat('See plot\n')
    drawdag(scenario2.DAG)
    return(invisible(2))}

  if ((p.v.X > conflevel)&(p.v.Y > conflevel))
  {cat("It seems that neither X or Y affect Z\nYou may want to review your working model\n")
  return(invisible(0))}
  
  if ((p.v.Y <= conflevel)&(p.v.X > conflevel))
  {cat('It looks like Y is related to Z, but not Z\nYou may want to revisit the hypothesis \'X = common cause of Y and Z\'')
  return(invisible(0))}
  }

a <- Y_check(non.influences)
b <- Y_check(Yinfluences)
c <- Y_check(Ynoinfluences)
#d <- Y_check(wrongcommoncause)

#Another sensitive question that may arise before the analysis is if we have identified
#the common cause correctly: the objective of the analysis is to 

#We will now create a toy example to illustrate the problems of a bad modeling of 
#scenario 1.
#We will study the expression of gene INK4a, key for melanoma development. It will be
#the outcome variable of Z.
#It is directly affected by UV radiation. UV radiation can come from sunbathing, which
#increases the appetite for ice cream consumption (measured in ml of consumed ice cream)
#You collect data of potentially cancerous tissue from 100 people, from which you know the
#hours they have spent in the sun the last year, the amount of consumed ice cream and
#the expression of INK4a.

b_xy_i_uv_i <- 5
b_xz_i_uv_i <- 10
b_yz_i_uv_i <- 0
samplesize <- 100

i_uv_i.DAG <- dagitty("dag {
UV.radiation -> Ice.cream.consumption
UV.radiation -> INK4a
}")
#coordinates(i_uv_i.DAG) <- list(x = c(Y = 1, X = 2, Z = 3),
                                   #y = c(Y = 3, X = 1, Z = 3))
drawdag(i_uv_i.DAG)

sc1.comm <- function(b_yz, N, b_xz, b_xy, reps = 100, ...){
  onlyY_pv <- rep(NA, reps)
  both_pv <- rep(NA,reps)
  
  onlyX_coefX <- rep(NA,reps)
  both_coefX <- rep(NA, reps)
  
  for (i in 1:reps) {
    
    dataset <- create.dataset(b_yz, N = N, b_xz = b_xz, b_xy = b_xy, ...)
    
    both <- lm(Z~X+Y, data = dataset)
    onlyY <- lm(Z~Y, data = dataset)
    onlyX <- lm(Z~X, data = dataset)
    
    onlyY_pv[i] <- summary(onlyY)$coefficients["Y", "Pr(>|t|)"]
    both_pv[i] <- summary(both)$coefficients["Y", "Pr(>|t|)"]
    
    onlyX_coefX[i] <- summary(onlyX)$coefficients["X", "Estimate"]
    both_coefX[i] <- summary(both)$coefficients["X", "Estimate"]
    
    rm(dataset)
  }
  cat('\n Change in relevance of Y on Z\n')
  cat('\nWhen Z ~ Y: \nThe p value of Y is ', mean(onlyY_pv),'\n')
  cat('\nWhen Z ~Y + X: \nThe p value of Y is ', mean(both_pv), '\n')
  
  cat('\n Change in effect of X over Z')

  cat('\nWhen Z ~ X: \nThe estimate for X is ', mean(onlyX_coefX),'and its s.d. is',
      sd(onlyX_coefX),'\n')
  cat('\nWhen Z ~Y + X: \nThe estimate for X is ', mean(both_coefX), 
      'and its s.d. is', sd(both_coefX),'\n')
  cat('\nBeing input x -> z: ', b_xz)
  #This illustrates how, even if the estimate of the coefficient for X is similar in both cases, the variance is higher in the presence of Y
  op <- par(mfrow= c(2,1), mar = rep(3,4))
  hist(onlyX_coefX, main = 'Z ~ X', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  hist(both_coefX, main = 'Z ~ X + Y', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  par(op)
  
  list <- list('onlyY_pv' = onlyY_pv, 'both_pv' = both_pv, 
               'onlyX_coefX' = onlyX_coefX, 'both_coefX' = both_coefX)
  invisible(list)
  }

sc1.comm(b_yz = b_yz_i_uv_i, N = samplesize, b_xz = b_xz_i_uv_i, b_xy = b_xy_i_uv_i)

#The ice cream consumption is not significant relevant to INK4a expression
#when Uv radiation is present in the model. UV radiation is a confounder.
#This corresponds to:
impliedConditionalIndependencies(i_uv_i.DAG)
#The power of UV radiation over INK4a doesn't vary much with the presence or absence of Y in the model.
#But its variance increases when Y is taken into account.
#Hence, in scenario 1 we should always condition on the common cause X or we could see 
#a fake, but significant, causal relation between Y and Z.

#A variation would be to consider that uv raditaion depends on sun exposure, 
#having the following DAG:
i_uv_i_sun.DAG <- dagitty("dag {
Sun -> UV.radiation
UV.radiation -> Ice.cream.consumption
UV.radiation -> INK4a
}")

#coordinates(i_uv_i_sun.DAG) <- list(x = c(Y = 1, X = 2, Z = 3),
                                #y = c(Y = 3, X = 1, Z = 3))
drawdag(i_uv_i_sun.DAG)

b_ax_i_uv_i2 <- 2
b_xy_i_uv_i2 <- 5
b_xz_i_uv_i2 <- 10
b_yz_i_uv_i2 <- 0
#As the significance of ice cream, Y, was covered in the previous function, 
#it will be skipped in this one. We are interested in knowing if the sun, A, plays a role 
#in the value of INK4a, Z, and how important is it.

sc1.comm.plusancestor <- function(b_yz, N, b_xz, b_xy, b_ax, reps = 30, e_x= 1, ...){
  onlyA_pvA <- rep(NA, reps)
  bothXA_pvA <- rep(NA,reps)
  three_pvA <- rep(NA, reps)
  
  onlyA_coefA <- rep(NA,reps)
  bothXA_coefA <- rep(NA, reps)
  three_coefA <- rep(NA, reps)
  
  onlyX_coefX <- rep(NA,reps)
  bothXA_coefX <- rep(NA, reps)
  three_coefX <- rep(NA, reps)
  
  onlyX_pvX <- rep(NA, reps)
  bothXA_pvX <- rep(NA, reps)
  three_pvX <- rep(NA, reps)
  
  #set.seed(13) #can be uncommented for reproducibility
  for (i in 1:reps) {
    
    dataset <- create.datasetv2(b_yz= b_yz, N = N, b_xz = b_xz, b_ax = b_ax,
                                b_xy = b_xy, e_x =e_x, ...)
    three <- lm(Z~X+A+Y, data = dataset)
    bothXA <- lm(Z~X+A, data = dataset)
    onlyA <- lm(Z~A, data = dataset)
    onlyX <- lm(Z~X, data = dataset)
    
    onlyA_pvA[i] <- summary(onlyA)$coefficients["A", "Pr(>|t|)"]
    bothXA_pvA[i] <- summary(bothXA)$coefficients["A", "Pr(>|t|)"]
    three_pvA[i] <- summary(three)$coefficients["A", "Pr(>|t|)"]
    
    onlyA_coefA[i] <- summary(onlyA)$coefficients["A", 'Estimate']
    bothXA_coefA[i] <- summary(bothXA)$coefficients["A", "Estimate"]
    three_coefA[i] <- summary(three)$coefficients["A", "Estimate"]
    
    onlyX_coefX[i] <- summary(onlyX)$coefficients["X", 'Estimate']
    bothXA_coefX[i] <- summary(bothXA)$coefficients["X", "Estimate"]
    three_coefX[i] <- summary(three)$coefficients["X", "Estimate"]
    
    onlyX_pvX[i] <- summary(onlyX)$coefficients["X", "Pr(>|t|)"]
    bothXA_pvX[i] <- summary(bothXA)$coefficients["X", "Pr(>|t|)"]
    three_pvX[i] <- summary(three)$coefficients["X", "Pr(>|t|)"]

    rm(dataset)
  }
  
  ###Changes in A
  #p valor de A en los modelos, es relevante o no
  cat('\n____Change in p value of A on Z\n')
  cat('\nWhen Z ~ A: \nThe p value of A is ', mean(onlyA_pvA),'\n')
  cat('\nWhen Z ~ X +A: \nThe p value of A is ', mean(bothXA_pvA), '\n')
  cat('\nWhen Z ~ Y + X + A: \nThe p value of A is ', mean(three_pvA), '\n')
  #estimate de A con y sin X
  cat('\n____Effect of A over Z\n')
  cat('Input A -> X: ', b_ax,'\nInput X -> Z:', b_xz,'\nTotal effect A -> Z', b_xz * b_ax, '\n')
  cat('\nWhen Z ~ A: \nCoefficient of A is ', mean(onlyA_coefA),'and its s.d. is',
      sd(onlyA_coefA),'\n')
  cat('\nWhen Z ~ X +A: \nCoefficient of A is ', mean(bothXA_coefA),'and its s.d. is',
      sd(bothXA_coefA),'\n')
  cat('\nWhen Z ~ Y + X + A: \nCoefficient of A is ', mean(three_coefA),'and its s.d. is',
      sd(three_coefA), '\nSee plots:\n')
  op <- par(mfrow= c(2,3), mar = rep(2,4))
  hist(onlyA_pvA, main = 'Z ~ A', xlab = 'p value of A')
  hist(bothXA_pvA, main = 'Z ~ X + A', xlab = 'p value of A')
  hist(three_pvA, main='Z ~X + A + Y', xlab = 'p value of A')
  
  hist(onlyA_coefA, main = 'Z ~ A', xlab = 'Effect A over Z')
  abline(v = b_xz*b_ax, col = 'red')
  hist(bothXA_coefA, main = 'Z ~ X + A', xlab = 'Effect A over Z')
  abline(v = b_xz*b_ax, col = 'red')
  hist(three_coefA, main='Z ~X + A + Y', xlab = 'Effect A over Z')
  abline(v = b_xz*b_ax, col = 'red')
  
  ###Changes in X
  #p value
  cat('\n____Change in p value of X on Z\n')
  cat('\nWhen Z ~ A: \nThe p value of A is ', mean(onlyX_pvX),'\n')
  cat('\nWhen Z ~ X +A: \nThe p value of A is ', mean(bothXA_pvX), '\n')
  cat('\nWhen Z ~ Y + X + A: \nThe p value of A is ', mean(three_pvX), '\n')
  
  #estimate de X en los modelos y error estandar
  cat('\n____Effect of X over Z\n\n')
  cat('Input X -> Z: ', b_xz,'\n')
  cat('\nWhen Z ~ A: \nCoefficient of X is ', mean(onlyX_coefX),'and its s.d. is',
      sd(onlyX_coefX),'\n')
  cat('\nWhen Z ~ X +A: \nCoefficient of X is ', mean(bothXA_coefX),'and its s.d. is',
      sd(bothXA_coefX),'\n')
  cat('\nWhen Z ~ Y + X + A: \nCoefficient of X is ', mean(three_coefX),'and its s.d. is',
      sd(three_coefX),'\nSee plots:\n')
  
  #p valor de X en los modelos, en histograma
  
  hist(onlyX_pvX, main = 'Z ~ X', xlab = 'p value of X')
  hist(bothXA_pvX, main = 'Z ~ X + A', xlab = 'p value of X')
  hist(three_pvX, main='Z ~X + A + Y', xlab = 'p value of X')
  
  hist(onlyX_coefX, main = 'Z ~ X', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  hist(bothXA_coefX, main = 'Z ~ X + A', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  hist(three_coefX, main='Z ~X + A + Y', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  
  par(op)
}

sc1.comm.plusancestor(b_yz = b_yz_i_uv_i2, N = samplesize, b_xz=b_xz_i_uv_i2,
                      b_ax = b_ax_i_uv_i2, b_xy = b_xy_i_uv_i2)
#From this it can be concluded that A is only significant when X is not in the model.
#The total effect of A is only appreciated in this model as well.
#From this it can be concluded that the adjustment of A is required only if we are interested
#on the effect of A over Z. By conditioning by X, A loses its relevance. This is supported by:
impliedConditionalIndependencies(i_uv_i_sun.DAG)
#The p value of X increases with the complexity of the model, but in any case it is significant.
#The estimate of X is around the expected even if complexity is increased, but its standard error gets higher.
#It the cause of study is X, conditioning by A is detrimental as the standard error of its coefficient
#increases, though its p value is never below any sensible significance level by adjusting by other
#covariates.

#As seen in the cause of the cause previous work from Ramón Díaz Uriarte, when the standard
#error of X increases the variance of the estimate doesn't change as much with the 
#presence of A on the model. We are still working on why this happens.
sc1.comm.plusancestor(b_yz = b_yz_i_uv_i2, N = samplesize, b_xz=b_xz_i_uv_i2,
                      b_ax = b_ax_i_uv_i2, b_xy=b_xy_i_uv_i2, e_x =10)

#This function also illustrates a key property of causal inference: same rules for
#simple models (toy example in Z_X_Y_adjust.R) can be applied to more complex models
#(as in this case).

###

#Let's move on to the scenario 2. We will illustrate with another example what to expect
#conditioning on the different possibilities. We have concluded that INK4a over expression
#is caused by UV radiation. A recent study shows that there is also intervention of 
#MATP in the French population in this process. It seems to follow the following
#DAG.
m_uv_i.DAG <- dagitty("dag {
UV.radiation -> MATP
UV.radiation -> INK4a
MATP -> INK4a
}")

#coordinates(i_uv_i_sun.DAG) <- list(x = c(Y = 1, X = 2, Z = 3),
#y = c(Y = 3, X = 1, Z = 3))
drawdag(m_uv_i.DAG)

samplesize <- 100
b_xz_m_uv_i <- 4
b_yz_m_uv_i <- (-3)
b_xy_m_uv_i <- 2


sc2.comm <- function(b_xz, b_yz, b_xy, N, reps = 50, ...) {
  #for simplicity, p value is not checked as in any model the covariates are causally related
  onlyX_coefX <- rep(NA,reps)
  both_coefX <- rep(NA, reps)
  onlyY_coefY <- rep(NA,reps)
  both_coefY <- rep(NA,reps)
  
  for (i in 1:reps){
    dataset <- create.dataset(N=N, b_xz = b_xz, b_yz = b_yz, b_xy = b_xy)
    both <- lm(Z~X + Y, dataset)
    onlyY <- lm(Z~Y, dataset)
    onlyX <- lm(Z~X, dataset)
    
    onlyX_coefX[i] <- summary(onlyX)$coefficients['X', 'Estimate']
    both_coefX[i] <- summary(both)$coefficients['X', 'Estimate']
    onlyY_coefY[i] <- summary(onlyY)$coefficients['Y', 'Estimate']
    both_coefY[i] <- summary(both)$coefficients['Y', 'Estimate']
    
    rm(dataset)
  }
  
  cat('___Changes in X\n')
  cat('When Z ~ X:\nX coefficient:', mean(onlyX_coefX), 's.d:', sd(onlyX_coefX))
  cat('\nWhen Z ~ X + Y:\nX coefficient:', mean(both_coefX), 's.d:', sd(both_coefX))
  
  cat('\n\n___Changes in Y\n')
  cat('When Z ~ Y:\nY coefficient:', mean(onlyY_coefY), 's.d:', sd(onlyY_coefY))
  cat('\nWhen Z ~ Y + X:\nY coefficient:', mean(both_coefY), 's.d:', sd(both_coefY))
  ###histogramas y revisar los efectos directos e indirectos
 }

sc2.comm(b_xz = b_xz_m_uv_i, b_yz = b_yz_m_uv_i, b_xy = b_xy_m_uv_i,
         N = samplesize)
#In this case, Y, thus, MATP, has a significant p value in both models, in presence
#and absence of the common cause X, UV radiation.

#The study goes on and we discovers that both the expression of MATP and INK4a
#is also influenced by another key factor, which is the cortisol level of the 
#individual. This leaves us with the following DAG:

m_uv_i_c.DAG <- dagitty("dag {
UV.radiation -> MATP
UV.radiation -> INK4a
cortisol -> MATP
cortisol -> INK4a
MATP -> INK4a
}")

#coordinates(i_uv_i_sun.DAG) <- list(x = c(Y = 1, X = 2, Z = 3),
#y = c(Y = 3, X = 1, Z = 3))
drawdag(m_uv_i_c.DAG)

samplesize <- 100
b_ay_m_uv_i_c <- 3
b_az_m_uv_i_c <- 2
b_xz_m_uv_i_c <- 4
b_yz_m_uv_i_c <- (-3)
b_xy_m_uv_i_c <- 2



###Backdoor criterion

#_____________________________EXPERIMENTS




Ynoinfluences <- create.dataset(0)
Yinfluences <- create.dataset(-2)
non.influences <- create.dataset(0, b_xz = 0)
#wrongcommoncause <- create.dataset(5, b_xz = 0)

#there are three possibilities to analyze the causal relations
summary(lm(Z~X+Y, data = Ynoinfluences)) #only X is significant
summary(lm(Z~X, data = Ynoinfluences)) #significant, total effect of X
summary(lm(Z~Y, data = Ynoinfluences)) #significant

summary(lm(Z~X+Y, data = Yinfluences)) #both significant, x = 3 e = 0.13, y = -2 e = 0.04
summary(lm(Z~X, data = Yinfluences)) #significant, x = -3 e = 0.003
summary(lm(Z~Y, data = Yinfluences)) #significant, Y = -1 e = 00


set.seed(18)
i_uv_i2 <- create.datasetv2(b_ax = 2, N = samplesize, b_xy = b_xy_i_uv_i2,
                            b_xz = b_xz_i_uv_i2, b_yz = b_yz_i_uv_i2)

onlyX <- lm(Z~X, data = i_uv_i2)
onlyY <- lm(Z~Y, data = i_uv_i2)
onlyA <- lm(Z~A, data = i_uv_i2)
YandX <- lm(Z~Y + X, data = i_uv_i2)
XandA <- lm(Z~X + A, data = i_uv_i2)
YandA <- lm(Z~Y + A, data = i_uv_i2)
YXA <- lm(Z~X +A +Y, data = i_uv_i2)
summary(onlyX) #significant, e = 10, s.e = 0.001
summary(onlyY) #significant, e = 2, s.e = 0.001
summary(onlyA) #significant, e = 20, s.e = 0.03
#summary(YandX) #X is significant, Y is not; Y = 0.07 sd = 0.1; X = 9.6, s.d= 0.54
summary(XandA) #X is significant, A is slightly significant, X= 10 s.d = 0.1, A = -0.05 sd = 0.2
summary(YandA) #Y is strongly significant, A not at all, depends on the size of A, Y = 1.9 sd = 0.04, A = 0.59 sd = 0.04
summary(YXA) #X is significant, A is slightly, X = 9.65 err 0.55, A = -0.05 err 0.2, Y = 0.07 err 0.1

impliedConditionalIndependencies(X.is.cause.DAG)
Y.noin.condboth <- lm(Z~X + Y, data = Ynoinfluences)
Y.noin.condY <- lm(Z~Y, data = Ynoinfluences)
Y.noin.condX <- lm(Z~X, data = Ynoinfluences)
summary(Y.noin.condboth)
summary(Y.noin.condY)
summary(Y.noin.condX)


#Should we condition for the cause of the common cause?
impliedConditionalIndependencies(i_uv_i_sun.DAG)






#Let's focus first on the case in which Y doesn't has a causal relationship with Z
#We will work on a toy example where X, or common cause is 'stress level', Z is
#'cholesterol' and Y is 'Cognition'. The value of all of them is affected by unmeasured covariates.

stress_chol_cog.DAG <- dagitty("dag {
Stress -> Cognition
Stress -> Cholesterol
U_1 -> Stress
U_2 -> Cognition
U_3 -> Cholesterol
}")

#coordinates(stress_chol_cog.DAG) <- list(x = c(Y = 1, X = 2, Z = 3),
                                  #y = c(Y = 3, X = 1, Z = 3))
drawdag(stress_chol_cog.DAG)

N = 100
S_C_H <- data.frame('Stress'= runif(N, 1, 100))
S_C_H$Cognition <- rnorm(N) - S_C_H$Stress * 2
S_C_H$Cholesterol <- S_C_H$Stress * 6 + rnorm(N)

summary(lm(Cholesterol~Cognition + Stress, data = S_C_H))

#el default sera q colesterol y cognicion no estan relacionados (b_yz)
S_C_H <- function(N = 100, reps = 1, signlevel = 0.01, b_yz = 0, ... )
                   {
  pv_st_stCg <- rep(NA, reps)
  pv_st_st <- rep(NA, reps)
  pv_cg_stCg <- rep(NA, reps)
  pv_cg_cg <- rep(NA, reps)
  
  for (i in 1:reps) {
    
    name_df <- create.dataset(b_yz = b_yz, N = N, ...)

    chol.stCh <- lm(Z~X+Y, data = name_df)
    chol.st <- lm(Z~X, data = name_df)
    chol.ch <- lm(Z~Y, data = name_df)
    
    pv_st_stCg[i] <- summary(chol.stCh)$coefficients["X", "Pr(>|t|)"]
    pv_cg_stCg[i] <- summary(chol.stCh)$coefficients["Y", "Pr(>|t|)"]
    pv_st_st[i] <- summary(chol.st)$coefficients["X", "Pr(>|t|)"]
    pv_cg_cg[i] <- summary(chol.ch)$coefficients["Y", "Pr(>|t|)"]
    
    rm(name_df)
  }

  op <- par(mfrow = c(2, 2))
  
  hist(pv_st_stCg, main = "Sig. stress when stress + cognition", xlab = "p-value")
  abline(v = signlevel)
  
  hist(pv_st_st, main = "Sig. stress when NO cognition", xlab = "p-value", xlim = c(0,1))
  abline(v = signlevel)
  
  hist(pv_cg_stCg, main = "Sig. cognition when stress + cognition", xlab = "p-value")
  abline(v = signlevel)
  
  hist(pv_cg_cg, main = "Sig. cognition when NO stress", xlab = "p-value")
  abline(v = signlevel)
  
  par(op)
}
S_C_H()  