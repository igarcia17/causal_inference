
#First, it is necessary to identify all the non-causal paths between X i and X j.

# i -> 3 -> 1 -> 4 -> j
# i -> 3 -> 1 -> 4 -> 2 -> 5 -> j
# i -> 4 -> 2 -> 5 -> j
# i -> 4 -> j
#Only X_4 is present in all the paths: conditioning on it will block all the backdoor
#paths. But we have to take into account that X_4 is a collider. It will be 
#necessary to condition on one of their ascendants, or a descendant of the ascendant.
#It gives 4 options: X_4 with X_1, X_2, or X_5, X_3.

#This corresponds to:
adjustmentSets(complex.DAG, "X_i", "X_j")

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
#It looks very alike the example with cortisol, UV radiation, INK4a and MATP
#(common cause, scenario 2, complex).
#To analyse it, it is basically the same as function sc2.comm.extraoverY
#as it has the same causal structure. The backdoor criteria tells us to condition
#on the same covariates as it was seen previously, in an experimental way, that the
#best estimates with less variance are obtained.
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
#estimate for smoking doesn't change its sign.

