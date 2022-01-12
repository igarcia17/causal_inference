
library(simstudy)
library(dagitty)
#library(car)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}


#________________________Berkson's paradox: 


# This paradox is a particular kind of selection bias, caused by systematically 
# observing some events more than other.

# Let's take a look at the DAG for our particular example, where we observe 
# hospital patients with diabetes and hospital patients with cholecystitis.


DAG.Berkson <- dagitty("dag {
diabetes -> hosp_patient
cholycistitis -> hosp_patient
}")

coordinates(DAG.Berkson) <- list(x = c(diabetes = 1,  cholycistitis = 3, 
                                       hosp_patient = 2),
                                 y = c(diabetes = 1, cholycistitis = 1, 
                                       hosp_patient = 3))
drawdag(DAG.Berkson)

# Let us generate the data according to the DAG:

set.seed(11)


N=100

diabetes <- rbinom(n=100, size=1, prob=0.7)
cholecystitis <- rbinom(n=100, size=1, prob=0.3) 
hosp_patients <- diabetes*0.6 + cholecystitis*0.4 + rnorm(N, mean = 0, sd = 0.2)

# Let's take a look at our data estimates, for that we will use a 
# function that creates a dataframe containing the p values and estimates for 
# the variables diabetes and cholecystitis and also diabetes and cholecystitis 
# adjusting by the hospital patients. It also returns a string checking the 
# significance of the pvalue for diabetes and cholecystitis and another one
# for the same value but when adjusting by the hospital patients.

df_pval_estimates <- function(var1, var2, var3) {
  
  dataset <- data.frame(var1,var2,var3)
  
  est_v1_v2 <- summary(glm(var1~var2, data = dataset, family = 'binomial')
  )$coefficients['var2', 'Estimate']
  pval_v1_v2 <- summary(glm(var1~var2, data = dataset, family = 'binomial')
  )$coefficients['var2', 'Pr(>|z|)']
  
  
  est_v1_v2_v3 <- summary(glm(var1~var2+var3, data=dataset, 
                              family = 'binomial'))$coefficients['var2', 'Estimate']
  pval_v1_v2_v3 <-summary(glm(var1~var2+var3, data=dataset, 
                              family = 'binomial'))$coefficients['var2', 'Pr(>|z|)']
  
  df_final <- data.frame('values_of' = c('Estimates', 'P_values'), 'diab_chol' = 
                           c(est_v1_v2, pval_v1_v2), 'diab_chol_HP' = 
                           c(est_v1_v2_v3, pval_v1_v2_v3))
  
  print(df_final)
  cat('\n')
  
  {if (pval_v1_v2 > 0.05)
    cat('The variables diabetes and cholecystitis do not show an association. The p value is',
        pval_v1_v2, 'and the estimate is', est_v1_v2, 'for this model', '\n')
    else
      cat('The variables diabetes and cholecystitis do show an association. The p value is',
          pval_v1_v2, 'and the estimate is', est_v1_v2, 'for this model', '\n')}
  
  {if (pval_v1_v2_v3 > 0.05)
    cat('The variables diabetes and cholycistitis do not show an association when adjusting by hospital patients. The p value is',
        pval_v1_v2_v3, 'and the estimate is',est_v1_v2_v3, 'for this model', '\n')
    else
      cat('The variables diabetes and cholycistitis do show an association when adjusting by hospital patients. The p value is',
          pval_v1_v2_v3, 'and the estimate is',est_v1_v2_v3, 'for this model', '\n')}
}

df_pval_estimates(diabetes, cholecystitis, hosp_patients)

#_______________________________Simpson paradox and backdoor criterion
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
samplesize <- 100

var.sc1 <- defData(varname = 'U', dist = 'binary', formula = 0.5)
var.sc1 <- defData(var.sc1, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc1 <- defData(var.sc1, varname = 'LBW', dist = 'binary', formula = '0.5 * U + 0.4 * Smoking', link = 'identity')
var.sc1 <- defData(var.sc1, varname = 'Mortality', dist = 'binary', formula = '0.1 * Smoking + 0.7 * U', link = 'identity')
#two seed to have the data on both dataframes as similar as possible
set.seed(13)
LBWsc1.df <- genData(samplesize, var.sc1)

var.sc2 <- defData(varname = 'U', dist = 'binary', formula = 0.5)
var.sc2 <- defData(var.sc2, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc2 <- defData(var.sc2, varname = 'LBW', dist = 'binary', formula = '0.5 * U + 0.4 * Smoking', link = 'identity')
var.sc2 <- defData(var.sc2, varname = 'Mortality', dist = 'binary', formula = '0.1 * Smoking + 0.7 * U + 0.1 * LBW', link = 'identity')

set.seed(13)
LBWsc2.df <- genData(samplesize, var.sc2)


LBW.fun <- function(dataset){
  #it doesn't check the p value of LBW because it is not the main issue of this type of
  #paradox; we want to focus on how the effect of Smoking changes in presence or absence of LBW

  onlyS <- glm(Mortality~Smoking, data = dataset, family = 'binomial')
  both <- glm(Mortality~LBW+Smoking, data = dataset, family = 'binomial')

  #This shows the estimate of Smoking
  cat('\nThe estimate for Smoking is:\nWhen Mortality ~ Smoking:', 
      summary(onlyS)$coefficients['Smoking','Estimate'],
      '\nWhen Mortality ~ Smoking + LBW:', 
      summary(both)$coefficients['Smoking','Estimate'], '\n')
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

var.sc3 <- defData(varname = 'Birth_defects', dist = 'binary', formula = 0.5)
var.sc3 <- defData(var.sc3, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc3 <- defData(var.sc3, varname = 'LBW', dist = 'binary', formula = '0.5 * Birth_defects + 0.4 * Smoking', link = 'identity')
var.sc3 <- defData(var.sc3, varname = 'Mortality', dist = 'binary', formula = '0.1 * Smoking + 0.7 * Birth_defects + 0.1 * LBW', link = 'identity')

set.seed(13)
LBWsc3.df <- genData(samplesize, var.sc3)
#It looks very alike the example with cortisol, UV radiation, INK4a and MATP.
#To analyse it, it is basically the same as function sc2.comm.extraoverY
#as it has the same causal structure
LBW.fun_noU <- function(dataset){
  onlyS <- glm(Mortality~Smoking, data = dataset, family = 'binomial')
  onlyBD <- glm(Mortality~Birth_defects, data = dataset, family = 'binomial')
  both <- glm(Mortality~Smoking+Birth_defects, data = dataset, family = 'binomial')
  
  cat('\nThe estimate for Smoking is:\nWhen Mortality ~ Smoking:', 
      summary(onlyS)$coefficients['Smoking','Estimate'],'\ns.d:', 
      summary(onlyS)$coefficients['Smoking','Std. Error'],
      '\nWhen Mortality ~ Smoking + Birth defects:', 
      summary(both)$coefficients['Smoking','Estimate'], '\ns.d:', 
      summary(both)$coefficients['Smoking','Std. Error'],'\n')
  cat('\nThe estimate for Birth defects is:\nWhen Mortality ~ Smoking:', 
      summary(onlyBD)$coefficients['Birth_defects','Estimate'],'\ns.d:', 
      summary(onlyBD)$coefficients['Birth_defects','Std. Error'],
      '\nWhen Mortality ~ Smoking + Birth defects:', 
      summary(both)$coefficients['Birth_defects','Estimate'], '\ns.d:', 
      summary(both)$coefficients['Birth_defects','Std. Error'],'\n')
}

LBW.fun_noU(LBWsc3.df)
#estimate for smoking varies slightly more: more test would be required to check why


##__________________________Backdoor criteria

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
X_3 <- X_1 * 5
X_4 <- X_1 * 2 + X_2 * 3
X_5 <- 11 * X_2
X_i <- 7 * X_3 + 13 * X_4
X_6 <- 17 * X_i
X_j <- X_6 * 23 + X_4 * 31 + X_5 * 19

#First, it is necessary to identify all the non-causal paths between X i and X j.
#By hand:
# i -> 3 -> 1 -> 4 -> j
# i -> 3 -> 1 -> 4 -> 2 -> 5 -> j
# i -> 4 -> 2 -> 5 -> j
# i -> 4 -> j
#By command:
paths(complex.DAG, from = 'X_i', to = 'X_j') #It is also included the causal path

#Only X_4 is present in all the paths, except the causal one: conditioning on it 
#will block all the backdoor
#paths. But we have to take into account that X_4 is a collider. It will be 
#necessary to condition on one of their ascendants, or a descendant of the ascendant.
#It gives 4 options: X_4 with X_1, X_2, or X_5, X_3.

identical(summary(lm(X_j ~ X_i + X_4))$coefficients['X_i'],
          summary(lm(X_j ~ X_i + X_4 + X_1))$coefficients['X_i'])
identical(summary(lm(X_j ~ X_i + X_4))$coefficients['X_i'],
          summary(lm(X_j ~ X_i + X_1))$coefficients['X_i'])

#This corresponds to:
adjustmentSets(complex.DAG, "X_i", "X_j")

#summary(lm(X_j ~ X_i))
#summary(lm(X_j ~ X_i + X_4))
#summary(lm(X_j ~ X_i + X_1))
#summary(lm(X_j ~ X_i + X_4 + X_1))
#summary(lm(X_j ~ X_i + X_4 + X_2))
#summary(lm(X_j ~ X_i + X_4 + X_3))
#summary(lm(X_j ~ X_i + X_4 + X_5))

#As the data analyst Motoharu Dei says,
#When doing data analysis, you have to know the causal structure of the subject 
#and use it properly. Otherwise, you may end up deriving a wrong insight. You 
#can't forget the causal context and reduce all the variables to statistical terms.
#Backdoor criteria cannot apply when subyacent DAG is not known.







