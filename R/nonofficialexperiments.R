library(simstudy)
library(dagitty)
#library(car)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}
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


