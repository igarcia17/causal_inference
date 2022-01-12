\name{Ycheck}
\alias{Ycheck, checking for descendant of common cause}
\title{
Y check: What is the relation between Y, descendant of X, and Z, variable of study
}
\description{
Receives a dataset an input and analyses the relevance of variable Y, that is a descendant of X, with Z, which is also a descendant of X, and the variable of the study. Returns the DAG of causal relationship among the three of them. Doesn't take into account other possible confounders.
}

\usage{
Check common cause scenario of dataset(x)
}

\arguments{
  \item{dataset}{
      Data frame that must include as column names X, Y and Z. TODO: improve this.   %%~~Describe \code{x} here~~
}
  \item{conflevel}{
     Confident level to which decide in which scenario are we in. Default = 0.01.
     %%~~Describe \code{x} here~~
}
}

\value{ --Nera: falta rellenar
  Returns the type of relationship between Y and variable of study Z and corresponding DAG~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
  Ines Garcia Ortiz
}


\examples{
## NOT directly executable
## The function is currently defined as
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
