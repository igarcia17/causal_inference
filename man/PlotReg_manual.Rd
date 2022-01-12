\name{PLOT.REG}
\alias{Plot Regression description function}
\title{Plot Regressions}
\description{Function to illustrate visually the distribution of the regression of two
variables with a scatterplot including a regression line.}

\arguments{
  \item{reg_line}
  Insert variables for regression, example: X ~ Y.
   \item{plot.var}
  Insert title for the graph.
}

\value{Returns a plot with the distribution of the regression model for the variables.}

\author{
Nermina Logo Lendo
}

\examples{
##---- Should NOT be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

N <- 500 
b_xz <- 3 
b_yz <- 2 
sd_z <- 1 


X <- runif(N, 1, 10)
Y <- runif(N, 1.5, 12)
Z <- b_xz * X + b_yz * Y + rnorm(N, 0, sd = sd_z)

PLOT.REG(Z ~ X, 'Var distribution X-Z')

## The function is currently defined as
PLOT.REG <- function(reg_line, plot.var = '') {
  Regres.line <- lm(reg_line) #generates the regression values
  plot(reg_line, main = plot.var) #plots this regression between the variables
  abline(Regres.line) #inserts the regression line in the plot
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