'Qué pasa si no se ajusta para la causa comun?'
library(car)
N <- 1e4

#primera columna Age, 1000 datos con media 30 y sd 5 que siga distribucion estandar
common_cause <- data.frame(Age = rnorm(N,30,5))
common_cause$Exercise <- 2 * common_cause$Age + rnorm(N)
common_cause$Cholesterol <- 3* common_cause$Age - common_cause$Exercise + rnorm(N)

m_common_cause_adjust <- lm(Cholesterol ~Exercise + Age, data = common_cause)
m_common_cause_no_adjust <- lm(Cholesterol~Exercise, data = common_cause)

S(m_common_cause_adjust)
S(m_common_cause_no_adjust)


#prueba, welch test o anova?
prueba <- data.frame(y = c(rnorm(1000, mean = 0.1, sd = 0.9), rnorm(5000, sd = 0.1)),
                     strain = factor(c(rep('mutado', times= 1000 ), rep('wt', times=5000))
                     ))
summary(prueba)
t.test(y~strain, data = prueba, var.equal = TRUE, alternative = 'two.sided')
t.test(y~strain, data = prueba, var.equal = FALSE, alternative = 'two.sided')

anova(lm(y~strain, data = prueba))
