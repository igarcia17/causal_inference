'Es buena idea ajustar para la causa de una causa?'
## Consider this DAG
## Z -> X  -> Y

## Plots, with dagitty
library(dagitty)

library(rethinking) ## for drawdag
## Installing rethinking can be complicated just for a few graphs
## So have a fallback if rethinking not available
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
    drawdag <- plot
} 

## Be explicit about the error terms
dag1 <- dagitty("dag {
Z -> X
X -> Y
e_x -> X
e_y -> Y
}")

coordinates(dag1) <- list(x = c(Z = 1, X = 2, Y = 3, e_x = 1.75, e_y = 2.75),
                          y = c(Z = 0, X = 0, Y = 0, e_x = -.2, e_y = -.2))

#con plot no se puede modificar el grafo para que quede bonito
## plot(dag1)
#drawdag aunq al principio quede feo luego se pueden poner limites para que se
#ajusten mejor las flechas de error
drawdag(dag1)
drawdag(dag1, xlim = c(0.5, 3.5), ylim = c(-2, 1))

#con esto se ve qué terminos son independientes entre si en el grafo
impliedConditionalIndependencies(dag1)


## This shows that adjusting for Z is actually a bad idea
## as it can increase the variance of the estimate of the effect of X

N <- 50
b_xz <- 3
b_yx <- 2
sd_x <- 1
sd_y <- 5

#haz un vector con tantos numeros como muestra haya, del 1 al 5
Z <- runif(N, 1, 5)
#X dependera de Z segun la relacion b_xz, y se le añade el error estandar 
#especificado antes
X <- b_xz * Z + rnorm(N, 0, sd = sd_x)
#Y dependera de X
Y <- b_yx * X + rnorm(N, 0, sd = sd_y)

summary(lm(Y ~ X))
summary(lm(Y ~ Z))
summary(lm(Y ~ X + Z))


dist_x <- function(N = 50, b_xz = 3, b_yx = 2, sd_x = 1, sd_y = 5, B = 2000,
                   z_min = 1, z_max = 5) {
    #inicializar vectores para almacenar los valores de X cuando se ajusta o no por Z
    X_noZ <- rep(NA, B)
    X_yesZ <- rep(NA, B)
    #inicializar vectores para almacenar p valores cuando se austa o no por Z
    pv_X_noZ <- rep(NA, B)
    pv_X_yesZ <- rep(NA, B)

    #por cada valor en la muestra, es decir, para cada caso   
    for(i in 1:B) {
      #X, Y y Z tendran un valor determinado
        Z <- runif(N, z_min, z_max)
        X <- b_xz * Z + rnorm(N, 0, sd = sd_x)
        Y <- b_yx * X + rnorm(N, 0, sd = sd_y)
        #se crean los modelos lineales teniendo en cuenta la causa
        m1 <- lm(Y ~ X)
        #y la causa de la causa
        m2 <- lm(Y ~ X + Z)
        #se almacena en los vectores los coeficientes
        X_noZ[i] <- coefficients(m1)["X"]
        X_yesZ[i] <- coefficients(m2)["X"]
        #se almacena en los vectores los p valor
        pv_X_noZ[i] <- summary(m1)$coefficients["X", "Pr(>|t|)"]
        pv_X_yesZ[i] <- summary(m2)$coefficients["X", "Pr(>|t|)"]
        #se borra Z, X y Y
        rm(Z, X, Y)
    }
    #se presenta un resumen de los coeficientes del modelo sin ajustar por la 
    #causa de una causa
    cat("\n Summary X without Z\n")
    print(summary(X_noZ))
    #se presenta el error estandar de los coeficientes
    cat("\n s.d. estimate = ", sd(X_noZ))
    #se presenta el resumen de los coeficientes cuando se ajusta por Z
    cat("\n\n Summary X with Z\n")
    print(summary(X_yesZ))
    #se presenta su error estandar
    cat("\n s.d. estimate = ", sd(X_yesZ), "\n")
    
    #se abre un plot device de 2 * 2
    op <- par(mfrow = c(2, 2), mar =c(1.5,1.5,1.5,1.5))
    #se hace un histograma de los coeficientes sin y con Z en el modelo
    hist(X_noZ, main = "X without Z in the model", xlab = "Estimate")
    abline(v = b_yx, lty = 2, col = 'red')
    hist(X_yesZ, main = "X with Z in the model", xlab = "Estimate")
    abline(v = b_yx, lty = 2, col = 'red')
    #se hace un histograma de los p valores con y sin Z en el modelo
    hist(pv_X_noZ, main = "X without Z in the model", xlab = "p-value")
    hist(pv_X_yesZ, main = "X with Z in the model", xlab = "p-value")
    #se restaura el plot device
    par(op)
    
}


dist_x(b_yx = -2, N = 20)

## Increasing sd_x minimizes problem
dist_x(b_yx = -2, sd_x = 5)

