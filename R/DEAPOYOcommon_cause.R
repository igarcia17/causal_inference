
#Let's create a function that creates the different datasets that we need.
create.dataset <- function(b_yz, N = 500, b_xy = 3, b_xz = 3,
                           e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}


#These functions will be handy later
create.datasetv2 <- function(b_ax, b_yz=0, N = 500, b_xy = 3, b_xz = 3,
                             e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(A = runif(N, 1, 100) + rnorm(N))
  name_df$X <- name_df$A * b_ax + rnorm(N, sd = e_x)
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}

create.datasetv3 <- function(b_by, b_bz, b_yz=(-3), N = 500, b_xy = 3, b_xz = 3,
                             e_x = 1, e_y = 1, e_z = 1, e_b = 1) {
  name_df <- data.frame(B = runif(N, 1, 100) + rnorm(N, sd = e_b))
  name_df$X <- runif(N, 1, 100) + rnorm(N, sd = e_x)
  name_df$Y <- name_df$X * b_xy + name_df$B * b_by +rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + 
    name_df$B * b_bz+ rnorm(N, sd = e_z)
  return(name_df)
}

Ynoinfluences <- create.dataset(0)
Yinfluences <- create.dataset(-2)
non.influences <- create.dataset(0, b_xz = 0)

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
    scenario1.DAG <- dagitty("dag {
    X -> Y
    X -> Z
    e_y -> Y
    e_z -> Z
    }")
    
    coordinates(scenario1.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                       y = c(Y = 3, X = 1, Z = 3, e_y = 2.75, e_z = 2.75))
    
    drawdag(scenario1.DAG)
    return(invisible(1))
    }
  
  if ((p.v.X <= conflevel)&(p.v.Y <= conflevel))
  {cat("The variable of analysis is influenced by both X and Y\n")
    cat('See plot\n')
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
coordinates(i_uv_i.DAG) <- list(x = c(Ice.cream.consumption = 1, UV.radiation = 2, INK4a = 3),
 y = c(Ice.cream.consumption = 3, UV.radiation = 1, INK4a = 3))

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
b_ax_i_uv_i2 <- 2
b_xy_i_uv_i2 <- 5
b_xz_i_uv_i2 <- 10
b_yz_i_uv_i2 <- 0

i_uv_i_sun.DAG <- dagitty("dag {
Sun -> UV.radiation
UV.radiation -> Ice.cream.consumption
UV.radiation -> INK4a
}")

coordinates(i_uv_i_sun.DAG) <- list(x = c(Ice.cream.consumption = 1, UV.radiation = 2, Sun = 2, INK4a = 3),
                                y = c(Ice.cream.consumption = 3, UV.radiation = 2, Sun = 1, INK4a = 3))
drawdag(i_uv_i_sun.DAG)

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

#As seen in the cause of the cause previous work from Ram?n D?az Uriarte, when the standard
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
samplesize <- 100
b_xz_m_uv_i <- 4
b_yz_m_uv_i <- (-3)
b_xy_m_uv_i <- 2

m_uv_i.DAG <- dagitty("dag {
UV.radiation -> MATP
UV.radiation -> INK4a
MATP -> INK4a
}")

coordinates(m_uv_i.DAG) <- list(x = c(MATP = 1, UV.radiation = 2, INK4a = 3),
                                y = c(MATP = 3, UV.radiation = 1, INK4a = 3))
drawdag(m_uv_i.DAG)

#In this case, we are not just asking the question how UV radiation (X) influences
#INK4a (Z), but we may also be interested in how MATP (Y) affects INK4a.

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

sc2.comm(b_xz = b_xz_m_uv_i, b_yz = b_yz_m_uv_i, b_xy = b_xy_m_uv_i,
         N = samplesize)

#In this case, Y, thus, MATP, has a significant p value in both models, in presence
#and absence of the common cause X, UV radiation. This makes sense and was expected.

#On blue it is shown the total effect of X over Z, as it takes into account the effect
#of X over Y as well. On red, it is shown the direct effect of X over Z. On green, the effect of Y over 
#Z. 
#When the mediator Y is out of the model, the X estimates the total effect over Z,
#including Y contribution. Only when Y is included it is possible to discern what is the
#direct effect of X.
#Depending on the case it would be more interesting to study the total or the direct
#effect of X. For this case, we argue that to know the total effect would be better
#because in the human body MATP expression is unavoidable.
#In any case, when there are two covariates in the model the variance of the estimate increases

#To know the effect of Y over Z, X has to be taken into account. When X is not 
#considered, the estimate for Y is biased; for this reason in this kind of graphs 
#X is called a confounder.

#To condition on Y or not may give unexpected outcomes when looking at X over Z.
#In this case, if MATP is not considered, it seems that the total effect is negative.
when_xz_4 <- create.dataset(b_xz = b_xz_m_uv_i, b_yz = b_yz_m_uv_i, b_xy = b_xy_m_uv_i,
                            N = samplesize)
scatterplot(Z~X, data = when_xz_4, main ='Original case', regLine=TRUE)

#If we input a higher X->Z value
sc2.comm(b_xz = b_xz_m_uv_i*10, b_yz = b_yz_m_uv_i, b_xy = b_xy_m_uv_i,
         N = samplesize)
when_xz_40 <- create.dataset(b_xz = b_xz_m_uv_i*10, b_yz = b_yz_m_uv_i, b_xy = b_xy_m_uv_i,
                             N = samplesize)
scatterplot(Z~X, data = when_xz_40, main = 'If UV radiation effect is stronger', regLine=TRUE)
#This is because the X contribution has a higher impact over Z than Y in this second case.

#The estimate for X in the simpler model isn't negative, it's total effect is lower
#than the direct effect.
#If the Y -> Z value wasn't negative:
sc2.comm(b_xz = b_xz_m_uv_i*10, b_yz = b_yz_m_uv_i*(-1), b_xy = b_xy_m_uv_i,
         N = samplesize)
when_xz_40andnegative <- create.dataset(b_xz = b_xz_m_uv_i*10, b_yz = b_yz_m_uv_i*(-1), b_xy = b_xy_m_uv_i,
                             N = samplesize)
#scatterplot(Z~X, data = when_xz_40andnegative, main = 'If MATP enhances INK4a', regLine=TRUE)

#Then the effect of X, UV radiation, over Z, INK4a, is increased, as it should be obvious.

#The study goes on and we discovers that both the expression of MATP and INK4a
#is also influenced by another key factor, the cortisol level. This leaves us 
#with the following DAG:

m_uv_i_c.DAG <- dagitty("dag {
UV.radiation -> MATP
UV.radiation -> INK4a
Cortisol -> MATP
Cortisol -> INK4a
MATP -> INK4a
}")

coordinates(m_uv_i_c.DAG) <- list(x = c(UV.radiation = 1, MATP = 2, INK4a = 2, Cortisol = 3),
                                y = c(UV.radiation = 1, MATP = 2, INK4a = 3, Cortisol = 1))
drawdag(m_uv_i_c.DAG)

samplesize <- 100
b_by_m_uv_i_c <- 3
b_bz_m_uv_i_c <- 2
b_xz_m_uv_i_c <- 4
b_xy_m_uv_i_c <- 2
b_yz_m_uv_i_c <- (-3)

#From the previous function we have learnt that the condition on the mediator
#variable Y, MATP, would allow us to know the direct effect of UV.radiation. If it is not 
#present in the model, what we can see is the total effect.
#The reasoning behing the UV.radiation-MATP-INK4a set also apply to the Cortisol-MATP-INK4a set.
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

sc2.comm.extraoverY(b_by = b_by_m_uv_i_c, b_bz = b_bz_m_uv_i_c,
                     b_xz = b_xz_m_uv_i_c, b_xy = b_xy_m_uv_i_c,
                     b_yz = b_yz_m_uv_i_c, N = samplesize)
#The total effect of cortisol, B, is well reflected when it is on its own in the model or 
#when UV radition, X, is considered, because they are independent from each other
#as can be seen in
impliedConditionalIndependencies(m_uv_i_c.DAG)
#The variance of the coefficient when it is found along X is smaller than when B (cortisol)
#is checked on its own or when it also considers the collider Y (MATP).
#Therefore in this type of graph it would be prefered to consider both B and X on the model:
#the variance of the estimates is smaller and the toal effect is calculated.
#As in the previous case, condiotioning on Y may be counterproductive.
#The absence of unmeasured confounding for the influence of both the exposure 
#and the intermediate variable on the outcome is required for estimating direct 
#effects. If both of these requirements are not satisfied, no approach can offer
#unbiased estimates of exposure's direct effects.