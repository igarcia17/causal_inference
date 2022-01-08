#Simpson paradox and backdoor criterion
library(simstudy)
library(dagitty)
library(car)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}
#Simpson paradox -> extreme case of confounding (though this is a simplicity of what he
#originally published)
#It is given when a trend appears in multiple data groups but
#disappears or is reversed when the groups are combined.

#Example of Simpson's paradox and unmeasured variables problem.

#Low birth weight paradox: 
LBWsc1.DAG <- dagitty('dag {
Smoking -> LBW
Smoking -> Mortality
U -> LBW
U -> Mortality
}')

coordinates(LBWsc1.DAG) <- list(x = c(LBW = 1, Smoking = 2, U = 2, Mortality = 3),
                                y = c(LBW = 3, Smoking = 2, U = 1, Mortality = 3))

drawdag(LBWsc1.DAG)

LBWsc2.DAG <- dagitty("dag {
Smoking -> LBW
LBW -> Mortality
Smoking -> Mortality
U -> LBW
U -> Mortality
}")

coordinates(LBWsc2.DAG) <- list(x = c(LBW = 1, Smoking = 2, U = 2, Mortality = 3),
                                y = c(LBW = 3, Smoking = 2, U = 1, Mortality = 3))

drawdag(LBWsc2.DAG)

#If U is not taken into consideration, same scenarios as in common cause analysis.
#An unmeasured variable can't be controlled in a model, and may represent for causes that
#we don't even consider in the experiment. In the low birth weight paradox, there are 
#unmeasured causes of low birth weight and mortality, though the only cause that we are accounting for
#is smoking. The questions are: does smoking cause mortality? Does LBW cause mortality?
#Does smoking cause LBW? Let's simulate the data: we will increase the natural proba-
#bilities to maximise the effect. Let's assume that U is an unkown health condition, and
#can be 0 or 1 in the absence or presence.
samplesize <- 10000

var.sc1 <- defData(varname = 'U', dist = 'binary', formula = 0.5)
var.sc1 <- defData(var.sc1, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc1 <- defData(var.sc1, varname = 'LBW', dist = 'binary', formula = '0.5 * U + 0.4 * Smoking', link = 'identity')
var.sc1 <- defData(var.sc1, varname = 'Mortality', dist = 'binary', formula = '0.05 * Smoking + 0.95 * U', link = 'identity')

set.seed(13)
LBWsc1.df <- genData(samplesize, var.sc1)

var.sc2 <- defData(varname = 'U', dist = 'binary', formula = 0.5)
var.sc2 <- defData(var.sc2, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc2 <- defData(var.sc2, varname = 'LBW', dist = 'binary', formula = '0.5 * U + 0.4 * Smoking', link = 'identity')
var.sc2 <- defData(var.sc2, varname = 'Mortality', dist = 'binary', formula = '0.01 * Smoking + 0.7 * U + 0.1 * LBW', link = 'identity')

set.seed(13)
LBWsc2.df <- genData(samplesize, var.sc2)


LBW.fun <- function(dataset, reps = 100){
  onlyS_coefS <- rep(NA, reps)
  both_coefS <- rep(NA, reps)
  
  for (i in 1:reps){
    onlyS <- glm(Mortality~Smoking, data = dataset, family = 'binomial')
    both <- glm(Mortality~LBW+Smoking, data = dataset, family = 'binomial')
    
    onlyS_coefS[i] <- summary(onlyS)$coefficients['Smoking','Estimate']
    both_coefS[i] <- summary(both)$coefficients['Smoking','Estimate']
    
  }

  #This shows the estimate of Smoking
  cat('\nThe estimate for Smoking is:\nWhen Mortality ~ Smoking:', mean(onlyS_coefS),
      '\nWhen Mortality ~ Smoking + LBW:', mean(both_coefS), '\n')
  ##poner scatterplots en una sola imagen, y que la leyenda del seguno sea mas descriptiva o no aparezca
  op <- par(mfrow= c(2,1), mar = rep(3,4))
  scatterplot(Mortality ~ Smoking, data = dataset, main = 'Mortality ~ Smoking', smooth = FALSE, regline = FALSE)
  scatterplot(Mortality ~ Smoking + LBW, data = dataset, main = 'Mortality ~ Smoking + LBW', smooth = FALSE, regline = FALSE)
  par(op)
  }
  
LBW.fun(LBWsc1.df)
LBW.fun(LBWsc2.df)
#Because low birth weight can be caused by other variables rather than smoking,
#that can be more dangerous to the child's health, when conditioning on LBW it is seen
#that Smoking has a negative correlation with mortality, thus, protecting the baby
#of dying. As in the previous cases, coditioning on Y, LBW, is a bad idea.

#If U could be measured as birth defects, it would give the following DAG:
LBWsc3.DAG <- dagitty("dag {
Smoking -> LBW
LBW -> Mortality
Smoking -> Mortality
Birth.defects -> LBW
Birth.defects -> Mortality
}")

coordinates(LBWsc3.DAG) <- list(x = c(LBW = 1.5, Smoking = 1, Birth.defects = 1, Mortality = 3),
                                y = c(LBW = 2, Smoking = 3, Birth.defects = 1, Mortality = 2))

drawdag(LBWsc3.DAG)
#It looks very alike the example with cortisol, UV radiation, INK4a and MATP.
#To analyse a more challenging DAG, we need to set a criteria to avoid being 
#fooled by paradoxes like this.

#The backdoor criteria tries to block all the backdoor paths between the outcome
#variable (X j) and the 'treatment' variable (X i). A backdoor path is a
#a non causal path between X and Y that has an arrow pointing to X i. We will 
#have to adjust for the set of variables Z when:
#(1) All back door paths between X and Y are blocked after conditioning on Z -> that is,
#as seen before, when Z is a mediator or a confounder or not a collider
#(2) No variables in Z are descendants of X.
#Given the DAG:
complex.DAG <- dagitty('dag{
X_1 -> X_3
X_1 -> X_4
X_2 -> X_4
X_2 -> X_5
X_3 -> X_i
X_4 -> X_i
X_4 -> X_j
X_5 -> X_j
X_i -> X_6
X_6 -> X_j
}')

coordinates(complex.DAG) <- list(x = c(X_1 =1, X_3 = 1, X_i =1,
                                       X_4 = 2, X_6 = 2,
                                       X_2 = 3, X_5 = 3, X_j =3),
                                 y = c(X_1=1, X_2 = 1,
                                       X_3 =2, X_4 = 2, X_5 = 2,
                                       X_i = 3, X_6 = 3, X_j =3))
drawdag(complex.DAG)

X_1 <- runif(100, 1, 100)
X_2 <- runif(100, 3, 30)
X_3 <- X_1 * 3
X_4 <- X_1 * 2 + X_2 * 5
X_5 <- 4 * X_2
X_i <- 0.2 * X_3 + 1.5 * X_4
X_6 <- 2 * X_i
X_j <- X_6 * 2 + X_4 * 3 + X_5 * 5

#First, it is necessary to identify all the non-causal paths between X i and X j.
paths(complex.DAG, from = 'X_i', to = 'X_j')
i -> 3 -> 1 -> 4 -> j
i -> 3 -> 1 -> 4 -> 2 -> 5 -> j
i -> 4 -> 2 -> 5 -> j
i -> 4 -> j
#Only X_4 is present in all the paths: conditioning on it will block all the backdoor
#paths. But we have to take into account that X_4 is a collider. It will be 
#necessary to condition on one of their ascendants, or a descendant of the ascendant.
#It gives 4 options: X_4 with X_1, X_2, or X_5, X_3.

identical(summary(lm(X_j ~ X_i + X_4))$coefficients['X_i'],summary(lm(X_j ~ X_i + X_4 + X_1))$coefficients['X_i'])
#The results corresponds to:
adjustmentSets(complex.DAG, "X_i", "X_j")


summary(lm(X_j ~ X_i))
summary(lm(X_j ~ X_i + X_4))
summary(lm(X_j ~ X_i + X_4 + X_1))
summary(lm(X_j ~ X_i + X_4 + X_2))
summary(lm(X_j ~ X_i + X_4 + X_3))
summary(lm(X_j ~ X_i + X_4 + X_5))
#As the data analyst Motoharu Dei says,
#When doing data analysis, you have to know the causal structure of the subject 
#and use it properly. Otherwise, you may end up deriving a wrong insight. You 
#can't forget the causal context and reduce all the variables to statistical terms








##__________EXPERIMENTS
summary(glm(Mortality~Smoking, data = LBWsc1.df, family = 'binomial'))
summary(glm(Mortality~Smoking+LBW, data = LBWsc1.df, family = 'binomial'))
##sin repeticiones
LBW.fun22 <- function(dataset){
  
  onlyS <- glm(Mortality~Smoking, data = dataset, family = 'binomial')
  onlyL <- glm(Mortality~LBW, data = dataset, family = 'binomial')
  both <- glm(Mortality~LBW+Smoking, data = dataset, family = 'binomial')
  
  onlyL_pvL <- summary(onlyL)$coefficients['LBW','Pr(>|z|)']
  both_pvL <- summary(both)$coefficients['LBW','Pr(>|z|)']
  
  onlyS_pvS <- summary(onlyS)$coefficients['Smoking','Pr(>|z|)']
  both_pvS <- summary(both)$coefficients['Smoking','Pr(>|z|)']
  onlyS_coefS <- summary(onlyS)$coefficients['Smoking','Estimate']
  both_coefS <- summary(both)$coefficients['Smoking','Estimate']
  
  #This checks is LBW is related or not
  cat('p value of LBW\nWhen Mortality ~ LBW:', onlyL_pvL, '\nWhen Mortality ~ LBW + Smoking', both_pvL,'\n')
  #Is Smoking relevant?
  cat('\np value of Smoking\nWhen Mortality ~ Smoking:', onlyS_pvS, '\nWhen Mortality ~ LBW + Smoking', both_pvS,'\n')
  
  #This shows the estimate of Smoking
  cat('\nThe estimate for Smoking is:\nWhen Mortality ~ Smoking', onlyS_coefS,
      '\nWhen Mortality ~ Smoking + LBW', both_coefS, '\n')
  scatterplot(Mortality ~ Smoking, data = dataset, main = 'Mortality ~ Smoking')
  scatterplot(Mortality ~ Smoking + LBW, data = dataset, main = 'Mortality ~ Smoking + LBW')
}
LBW.fun22(LBWsc1.df)
LBW.fun22(LBWsc2.df)

var.sc3 <- defData(varname = 'U', dist = 'binary', formula = 0.5)
var.sc3 <- defData(var.sc3, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc3 <- defData(var.sc3, varname = 'LBW', dist = 'binary', formula = '0.75 * U + 0.25 * Smoking', link = 'identity')
var.sc3 <- defData(var.sc3, varname = 'Mortality', dist = 'binary', formula = '0.05 * Smoking + 0.95 * U', link = 'identity')

set.seed(13)
LBWsc3.df <- genData(samplesize, var.sc3)
LBW.fun22(LBWsc3.df)

var.sc3 <- defData(varname = 'Birth_defects', dist = 'binary', formula = 0.5)
var.sc3 <- defData(var.sc3, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc3 <- defData(var.sc3, varname = 'LBW', dist = 'binary', formula = '0.5 * Birth_defects + 0.4 * Smoking', link = 'identity')
var.sc3 <- defData(var.sc3, varname = 'Mortality', dist = 'binary', formula = '0.1 * Smoking + 0.7 * Birth_defects + 0.1 * LBW', link = 'identity')

set.seed(13)
LBWsc3.df <- genData(samplesize, var.sc3)

summary(glm(Mortality ~ Smoking, data = LBWsc3.df, family = 'binomial'))
summary(glm(Mortality ~ Birth_defects, data = LBWsc3.df, family = 'binomial'))
summary(glm(Mortality ~ Smoking + Birth_defects, data = LBWsc3.df, family = 'binomial'))
summary(glm(Mortality ~ Birth_defects + Smoking + LBW, data = LBWsc3.df, family = 'binomial'))

