\name{sc1.comm}
\alias{sc1.comm.plusancestor, sc2.comm.extraoverY}
\title{
Description of common cause cases
}
\description{
Given the parameters that affect a set of variables X, Y, Z, (A, B), generate the plots that illustrate a change of the estimates and a terminal output of description.
}
\arguments{
  \item{b_yz}{ Coefficient of Y over Z}
  \item{b_xz}{ Coefficient of X over Z}
  \item{b_xy}{ Coefficient of X over Y}
  \item{b_ax}{ Coefficient of A over X}
  \item{b_by}{ Coefficient of B over ZY}
  \item{N}{Sample size for each time a dataset is done}
  \item{reps}{Number of times models are generated}
}

\value{ Console output of descriptive parameters to pay attention in the analysis.
}

\author{
Ines Garcia Ortiz
}

\examples{
##Not directly executable!!

## The function for simple scenario 1 is currently defined as
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
## The function for complex scenario 1 is currently defined as:
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
  cat('\nWhen Z ~ X + A: \nThe p value of A is ', mean(bothXA_pvA), '\n')
  cat('\nWhen Z ~ Y + X + A: \nThe p value of A is ', mean(three_pvA), '\n')
  #estimate de A con y sin X
  cat('\n____Effect of A over Z\n')
  cat('Input A -> X: ', b_ax,'\nInput X -> Z:', b_xz,'\nTotal effect A -> Z', b_xz * b_ax, '\n')
  cat('\nWhen Z ~ A: \nCoefficient of A is ', mean(onlyA_coefA),'and its s.d. is',
      sd(onlyA_coefA),'\n')
  cat('\nWhen Z ~ X + A: \nCoefficient of A is ', mean(bothXA_coefA),'and its s.d. is',
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
  par(op)
  ###Changes in X
  #p value         ######it is very obvious
  #cat('\n____Change in p value of X on Z\n')
  #cat('\nWhen Z ~ A: \nThe p value of A is ', mean(onlyX_pvX),'\n')
  #cat('\nWhen Z ~ X +A: \nThe p value of A is ', mean(bothXA_pvX), '\n')
  #cat('\nWhen Z ~ Y + X + A: \nThe p value of A is ', mean(three_pvX), '\n')
  
  #estimate de X en los modelos y error estandar
  cat('\n____Effect of X over Z\n\n')
  cat('Input X -> Z: ', b_xz,'\n')
  cat('\nWhen Z ~ A: \nCoefficient of X is ', mean(onlyX_coefX),'and its s.d. is',
      sd(onlyX_coefX),'\n')
  cat('\nWhen Z ~ X +A: \nCoefficient of X is ', mean(bothXA_coefX),'and its s.d. is',
      sd(bothXA_coefX),'\n')
  cat('\nWhen Z ~ Y + X + A: \nCoefficient of X is ', mean(three_coefX),'and its s.d. is',
      sd(three_coefX),'\nSee plots:\n')
  
  op <- par(mfrow= c(2,3), mar = rep(2,4))
  
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
## The function for simple scenario 2 is currently defined as:
sc2.comm <- function(b_xz, b_yz, b_xy, N, reps = 200, ...) {
  onlyY_pvY <- rep(NA, reps)
  both_pvY <- rep(NA, reps)
  onlyX_coefX <- rep(NA,reps)
  both_coefX <- rep(NA, reps)
  onlyY_coefY <- rep(NA,reps)
  both_coefY <- rep(NA,reps)
  
  for (i in 1:reps){
    dataset <- create.dataset(N=N, b_xz = b_xz, b_yz = b_yz, b_xy = b_xy)
    both <- lm(Z~X + Y, dataset)
    onlyY <- lm(Z~Y, dataset)
    onlyX <- lm(Z~X, dataset)
    
    onlyY_pvY[i] <- summary(onlyY)$coefficients['Y', 'Pr(>|t|)']
    both_pvY[i] <- summary(both)$coefficients['Y', 'Pr(>|t|)']
    onlyX_coefX[i] <- summary(onlyX)$coefficients['X', 'Estimate']
    both_coefX[i] <- summary(both)$coefficients['X', 'Estimate']
    onlyY_coefY[i] <- summary(onlyY)$coefficients['Y', 'Estimate']
    both_coefY[i] <- summary(both)$coefficients['Y', 'Estimate']
    
    #rm(dataset)
  }
  
  cat('p value of Y')
  cat('\nWhen Z~Y:', mean(onlyY_pvY))
  cat('\nWhen Z~Y+X', mean(both_pvY),'\n\n')
  
  cat('___Changes in X\n')
  cat('When Z ~ X:\nX coefficient:', mean(onlyX_coefX), 's.d:', sd(onlyX_coefX))
  cat('\nWhen Z ~ X + Y:\nX coefficient:', mean(both_coefX), 's.d:', sd(both_coefX))
  
  cat('\n\n___Changes in Y\n')
  cat('When Z ~ Y:\nY coefficient:', mean(onlyY_coefY), 's.d:', sd(onlyY_coefY))
  cat('\nWhen Z ~ Y + X:\nY coefficient:', mean(both_coefY), 's.d:', sd(both_coefY))
  ##legend: blue, direct effect X, red total effect X, green effect Y
  
  op <- par(mfrow= c(2,2), mar = rep(3,4))
  hist(onlyX_coefX, main = 'Z ~ X', xlab = 'Effect X over Z')
  abline(v = b_xz + b_yz*b_xy, col = 'blue')#total effect, it takes into account both sources of effect
  abline(v = b_xz, col = 'red')#direct effect
  
  hist(both_coefX, main = 'Z ~ X + Y', xlab = 'Effect X over Z')
  abline(v = b_xz + b_yz*b_xy, col = 'blue')#total effect
  abline(v = b_xz, col = 'red')#direct effect
  
  hist(onlyY_coefY, main = 'Z ~ Y', xlab = 'Effect X over Z')
  abline(v = b_yz, col = 'green')
  hist(both_coefY, main = 'Z ~ Y + X', xlab = 'Effect X over Z')
  abline(v = b_yz, col = 'green') 
  par(op)

 }
## The function for complex scenario 2 is currently defined as:
sc2.comm.extraoverY <- function(b_by, b_bz, b_xz, b_xy, b_yz, N, reps = 200, ...) {
  #que variables quiero ver ahora
  onlyB_coefB <- rep(NA,reps)
  bothBX_coefB <- rep(NA, reps)
  three_coefB <- rep(NA,reps)
  
  for (i in 1:reps){
    dataset <- create.datasetv3(N=N, b_xz = b_xz, b_yz = b_yz, b_xy = b_xy,
                                b_by = b_by, b_bz = b_bz, ...)
    onlyB <- lm(Z ~ B, dataset)
    bothBX <- lm(Z ~ X + B, dataset)
    three <- lm(Z ~ Y + X + B, dataset)
    
    onlyB_coefB[i] <- summary(onlyB)$coefficients['B', 'Estimate']
    bothBX_coefB[i] <- summary(bothBX)$coefficients['B', 'Estimate']
    three_coefB[i] <- summary(three)$coefficients['B', 'Estimate']

    rm(dataset)
  }

  cat('___Changes in B\n')
  cat('When Z ~ B:\nB coefficient:', mean(onlyB_coefB), 's.d:', sd(onlyB_coefB))
  cat('\nWhen Z ~ X + B:\nB coefficient:', mean(bothBX_coefB), 's.d:', sd(bothBX_coefB))
  cat('\nWhen Z ~ Y + X + B:\nB coefficient:', mean(three_coefB), 's.d:', sd(three_coefB))
  
 ##legend: blue, total efffect X, red direct effect X, green effect Y
  
  op <- par(mfrow= c(1,3), mar = rep(3,4))
  hist(onlyB_coefB, main = 'Z ~ B', xlab = 'Effect B over Z')
  abline(v = b_bz, col = 'red')#direct effect
  abline(v = b_bz + b_yz*b_by, col = 'blue')#total effect of B
  
  hist(bothBX_coefB, main = 'Z ~ B + X', xlab = 'Effect B over Z')
  abline(v = b_bz, col = 'red')#direct effect
  abline(v = b_bz + b_yz*b_by, col = 'blue')#total effect of B
  
  hist(three_coefB, main = 'Z ~ B + X + Y', xlab = 'Effect B over Z')
  abline(v = b_bz, col = 'red')#direct effect
  abline(v = b_bz + b_yz*b_by, col = 'blue')#total effect of B
  par(op)
}
{
  }
}