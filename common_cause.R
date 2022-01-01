#Confounder case: common cause

#import modules

library(dagitty)
library(rethinking)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
} 

#Consider the DAG:
comm.cause.DAG <- dagitty("dag {
X -> Y
X -> Z
Y -> Z
e_y -> Y
e_z -> Z
}")

coordinates(comm.cause.DAG) <- list(x = c(Y = 1, X = 2, Z = 3, e_y = 0.75, e_z = 2.75),
                                    y = c(Y = 1, X = 3, Z = 1, e_y = 0.75, e_z = 0.75))
drawdag(comm.cause.DAG)


