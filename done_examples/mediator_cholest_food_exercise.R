'Efecto de un mediador'

## Plots, with dagitty

library(dagitty)
library(car) ## For S function

library(rethinking) ## for drawdag
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
    drawdag <- plot
} 

#tenemos un grafo en el que el ejercicio es causa comun de comida y colesterol
#y colesterol es un collider de ejercicio y comida
mediator <- dagitty("dag {
Exercise -> Food
Food -> Cholesterol
Exercise -> Cholesterol
}")

drawdag(mediator)


N <- 1e4 ## huge N
#se crean mil datos para ejercicio, que vayan del 1 al 100
mediator <- data.frame(Exercise = runif(N, 1, 100))
#comida sera el fruto de un numero aleatorio mas cierta aportacion del ejercicio
#si haces más ejercicio, comes mas
#la parte aleatoria contribuye a otros factores que afectan a comida
mediator$Food <- mediator$Exercise * 3 + rnorm(N)
## Cholesterol decreases with Exercise (-1) but increases with Food
mediator$Cholesterol <- mediator$Food * 2 - mediator$Exercise + rnorm(N)

'Que pasa si ajustas por comida frente a si no ajustas'
m_mediator_adjust <- lm(Cholesterol ~ Exercise + Food, data = mediator)
m_mediator_no_adjust <- lm(Cholesterol ~ Exercise, data = mediator)

'no entiendo muy bien qué estoy viendo: por qué al NO ajustar por comida vemos
como estimador 5, que es lo que sí tiene en cuenta tambien el efecto
de la comida?'
S(m_mediator_adjust)
S(m_mediator_no_adjust)
## The estimate from the last model is what we would expect for the total effects:
## 3 * 2 - 1, where the 3 * 2 is the usual rule of product along paths.
