X_1 <- runif(100, 1, 100)
X_2 <- runif(100, 3, 30)
X_3 <- X_1 * 5
X_4 <- X_1 * 2 + X_2 * 3
X_5 <- 11 * X_2
X_i <- 7 * X_3 + 13 * X_4
X_6 <- 17 * X_i
X_j <- X_6 * 23 + X_4 * 31 + X_5 * 19

#First, it is necessary to identify all the non-causal paths between X i and X j.
paths(complex.DAG, from = 'X_i', to = 'X_j')
# i -> 3 -> 1 -> 4 -> j
# i -> 3 -> 1 -> 4 -> 2 -> 5 -> j
# i -> 4 -> 2 -> 5 -> j
# i -> 4 -> j
#Only X_4 is present in all the paths: conditioning on it will block all the backdoor
#paths. But we have to take into account that X_4 is a collider. It will be 
#necessary to condition on one of their ascendants, or a descendant of the ascendant.
#It gives 4 options: X_4 with X_1, X_2, or X_5, X_3.

identical(summary(lm(X_j ~ X_i + X_4))$coefficients['X_i'],summary(lm(X_j ~ X_i + X_4 + X_1))$coefficients['X_i'])
identical(summary(lm(X_j ~ X_i + X_4))$coefficients['X_i'],summary(lm(X_j ~ X_i + X_1))$coefficients['X_i'])

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
#It looks very alike the example with cortisol, UV radiation, INK4a and MATP.
#To analyse a more challenging DAG, we need to set a criteria to avoid being 
#fooled by paradoxes like this.
