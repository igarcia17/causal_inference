\name{LBW}
\alias{LBW description function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low Birth Weight Paradox estimates description%%  ~~function to do ... ~~
}
\description{
Function to illustrate the problem that arises in the Low Birth Weight paradox case, being it a Simpson's paradox example. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  Must contain variables: Smoking, LBW, Mortality, U.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ines Garcia Ortiz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(simstudy)
samplesize <- 10000
var.sc1 <- defData(varname = 'U', dist = 'binary', formula = 0.5)
var.sc1 <- defData(var.sc1, varname = 'Smoking', dist = 'binary', formula = 0.5)
var.sc1 <- defData(var.sc1, varname = 'LBW', dist = 'binary', formula = '0.5 * U + 0.4 * Smoking', link = 'identity')
var.sc1 <- defData(var.sc1, varname = 'Mortality', dist = 'binary', formula = '0.05 * Smoking + 0.95 * U', link = 'identity')
LBWsc1.df <- genData(samplesize, var.sc1)
## The function is currently defined as
LBW.fun <- function(dataset, reps = 100){
  onlyS_coefS <- rep(NA, reps)
  both_coefS <- rep(NA, reps)
  
  for (i in 1:reps){
    onlyS <- glm(Mortality~Smoking, data = dataset, family = 'binomial')
    both <- glm(Mortality~LBW+Smoking, data = dataset, family = 'binomial')
    
    onlyS_coefS[i] <- summary(onlyS)$coefficients['Smoking','Estimate']
    both_coefS[i] <- summary(both)$coefficients['Smoking','Estimate']
    
  }

  #This shows the estimate of Smoking
  cat('\nThe estimate for Smoking is:\nWhen Mortality ~ Smoking:', mean(onlyS_coefS),
      '\nWhen Mortality ~ Smoking + LBW:', mean(both_coefS), '\n')
  op <- par(mfrow= c(2,1), mar = rep(3,4))
  scatterplot(Mortality ~ Smoking, data = dataset, main = 'Mortality ~ Smoking', smooth = FALSE, regline = FALSE)
  scatterplot(Mortality ~ Smoking + LBW, data = dataset, main = 'Mortality ~ Smoking + LBW', smooth = FALSE, regline = FALSE)
  par(op)
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
