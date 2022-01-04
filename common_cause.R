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

Ynoinfluences <- create.dataset(0)
Yinfluences <- create.dataset(-2)
non.influences <- create.dataset(0, b_xz = 0)
#wrongcommoncause <- create.dataset(5, b_xz = 0)

#there are three possibilities to analyze the causal relations
summary(lm(Z~X+Y, data = Ynoinfluences)) #only X is significant
summary(lm(Z~X, data = Ynoinfluences)) #significant, total effect of X
summary(lm(Z~Y, data = Ynoinfluences)) #significant

summary(lm(Z~X+Y, data = Yinfluences)) #both significant, only direct effect of Z
summary(lm(Z~X, data = Yinfluences)) #significant, total effect of X
summary(lm(Z~Y, data = Yinfluences)) #significant

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
#We want to study the expression of gene INK4a, key for melanoma development. It will be
#the outcome variable or Z.
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
  onlyY_pv <- c(NA, reps)
  both_pv <- c(NA,reps)
  
  onlyX_coefX <- c(NA,reps)
  both_coefX <- c(NA, reps)
  
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
  
  cat('\n Change in effect of X over Z\n')
  cat('Input x -> z: ', b_xz,'\n')
  cat('\nWhen Z ~ X: \nThe estimate for X is ', mean(onlyX_coefX),'and its s.d. is',
      sd(onlyX_coefX),'\n')
  cat('\nWhen Z ~Y + X: \nThe estimate for X is ', mean(both_coefX), 
      'and its s.d. is', sd(both_coefX),'\n')
  
  op <- par(mfrow= c(2,1), mar = c(1.5,1.5,1.5,1.5))
  hist(onlyX_coefX, main = 'Z ~ X', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  hist(both_coefX, main = 'Z ~ X + Y', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  par(op)
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

sc1.comm.plusancestor <- function(b_yz, N, b_xz, b_xy, b_ax, reps = 30, ...){
  onlyA_pvA <- c(NA, reps)
  bothXA_pvA <- c(NA,reps)
  three_pvA <- c(NA, reps)
  
  onlyA_coefA <- c(NA,reps)
  bothXA_coefA <- c(NA, reps)
  three_coefA <- c(NA, reps)
  
  onlyX_coefX <- c(NA,reps)
  bothXA_coefX <- c(NA, reps)
  three_coefX <- c(NA, reps)
  
  onlyX_pvX <- c(NA, reps)
  bothXA_pvX <- c(NA, reps)
  three_pvX <- c(NA, reps)
  
  for (i in 1:reps) {
    
    dataset <- create.datasetv2(b_yz= b_yz, N = N, b_xz = b_xz, b_ax = b_ax,
                                b_xy = b_xy, ...)
    three <- lm(Z~X+A+Y)
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
  #p valor de A en los modelos, es relevante o no
  cat('\n Change in relevance of A on Z\n')
  cat('\nWhen Z ~ A: \nThe p value of A is ', mean(onlyA_pvA),'\n')
  cat('\nWhen Z ~ X +A: \nThe p value of A is ', mean(bothXA_pvA), '\n')
  cat('\nWhen Z ~ Y + X + A: \nThe p value of A is ', mean(three_pvA), '\n')
  #estimate de A con y sin X
  cat('\n Effect of A over Z\n')
  cat('Input A -> X: ', b_ax,'\n')
  cat('\nWhen Z ~ A: \nCoefficient of A is ', mean(onlyA_coefA),'\n')
  cat('\nWhen Z ~ X +A: \nCoefficient of A is ', mean(bothXA_coefA), '\n')
  cat('\nWhen Z ~ Y + X + A: \nCoefficient of A is ', mean(three_coefA), '\n')
  #Aqui poner unos histogramas bonitos
  
  #estimate de X en los modelos y error estandar
  cat('\n Effect of A over Z\n')
  cat('Input X -> Z: ', b_xz,'\n')
  cat('\nWhen Z ~ A: \nCoefficient of X is ', mean(onlyX_coefX),'and its s.d. is',
      sd(onlyX_coefX),'\n')
  cat('\nWhen Z ~ X +A: \nCoefficient of X is ', mean(bothXA_coefX),'and its s.d. is',
      sd(bothXA_coefX),'\n')
  cat('\nWhen Z ~ Y + X + A: \nCoefficient of X is ', mean(three_coefX),'and its s.d. is',
      sd(three_coefX),'\n')
  
  #p valor de X en los modelos, en histograma
  
  
  op <- par(mfrow= c(2,1), mar = c(1.5,1.5,1.5,1.5))
  hist(onlyX_coefX, main = 'Z ~ X', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  hist(both_coefX, main = 'Z ~ X + Y', xlab = 'Effect X over Z')
  abline(v = b_xz, col = 'red')
  par(op)
}

sc1.comm.plusancestor(b_yz = b_yz_i_uv_i2, N = samplesize, b_xz=b_xz_i_uv_i2,
                      b_ax = b_ax_i_uv_i2)

#But when the effect of A over X is greater, this changes:

sc1.comm.plusancestor(b_yz = b_yz_i_uv_i2, N = samplesize, b_xz=b_xz_i_uv_i2,
                      b_ax = 100)



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

set.seed(18)
i_uv_i2 <- create.datasetv2(b_ax = 200, N = samplesize, b_xy = b_xy_i_uv_i2,b_xz = b_xz_i_uv_i2, b_yz = b_yz_i_uv_i2)
onlyX <- lm(Z~X, data = i_uv_i2)
onlyY <- lm(Z~Y, data = i_uv_i2)
onlyA <- lm(Z~A, data = i_uv_i2)
YandX <- lm(Z~Y + X, data = i_uv_i2)
XandA <- lm(Z~X + A, data = i_uv_i2)
YandA <- lm(Z~Y + A, data = i_uv_i2)
YXA <- lm(Z~A +Y + X, data = i_uv_i2)
summary(onlyX) #significant, e = 10, s.e = 0.00
summary(onlyY) #significant, e = 2, s.e = 0.00
summary(onlyA) #significant, e = 2000, s.e = 0.03
#summary(YandX) #X is significant, Y is not; Y = 0.07 sd = 0.1; X = 9.6, s.d= 0.54
summary(XandA) #X is significant, A is slightly significant, X= 10 s.d = 0.1, A = -6 sd = 20
summary(YandA) #Y is strongly significant, A not at all, depends on the size of A, Y = 1.9 sd = 0.04, A = 59 sd = 41
summary(YXA) #X is significant, A is slightly, X = 9.65 err 0.55, A = -6 err 20, Y = 0.07 err 0.1


#Should we condition for the cause of the common cause?
impliedConditionalIndependencies(i_uv_i_sun.DAG)


#_____________________________























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
  


impliedConditionalIndependencies(X.is.cause.DAG)
Y.noin.condboth <- lm(Z~X + Y, data = Ynoinfluences)
Y.noin.condY <- lm(Z~Y, data = Ynoinfluences)
Y.noin.condX <- lm(Z~X, data = Ynoinfluences)
summary(Y.noin.condboth)
summary(Y.noin.condY)
summary(Y.noin.condX)