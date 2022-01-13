\name{SUM.2VAR}
\alias{Summary for 2 variables}
\title{Summary 2 variables}
\description{
This function recieves 2 variables and loops through each one of the p values
of the summary of a linear regression to check wheter there is a significant
association between them.
}
\usage{
SUM.2VAR(x)
}

\arguments{
  \item{variableX}
  We introduce our model of regression.
  \item{v_depX}
  Dependent variable whose estimates and pvalues we want to see.
  \item{title_varX}
  Title for the model (variables for the regression).
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ It returns a string containing the pvalues and estimates and it also
returns text describing whether there is association between the variables.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Nermina Logo Lendo}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should NOT be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.



## The function is currently defined as
SUM.2VAR <- function(variable1, v_dep1, title_var1, variable2, 
                     v_dep2, title_var2){
  PV1 <- summary(lm(variable1))$coefficients[v_dep1, 'Pr(>|t|)']
  E1 <- summary(lm(variable1))$coefficients[v_dep1, 'Estimate']
  
  PV2 <- summary(lm(variable2))$coefficients[v_dep2, 'Pr(>|t|)']
  E2 <- summary(lm(variable2))$coefficients[v_dep2, 'Estimate']
  
  P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, PV1, P.VAL.DICT)
  assign(title_var2, PV2, P.VAL.DICT)
  
  E.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, E1, E.DICT)
  assign(title_var2, E2, E.DICT)
  
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0.05)
      cat('The variables', i, 'do not show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show a correlation with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
          # SUMMARY for 2 variables

SUM.3VAR <- function(variable1, v_dep1, title_var1, variable2, 
                     v_dep2, title_var2, variable3, v_dep3, title_var3){
  PV1 <- summary(lm(variable1))$coefficients[v_dep1, 'Pr(>|t|)']
  E1 <- summary(lm(variable1))$coefficients[v_dep1, 'Estimate']
  
  PV2 <- summary(lm(variable2))$coefficients[v_dep2, 'Pr(>|t|)']
  E2 <- summary(lm(variable2))$coefficients[v_dep2, 'Estimate']
  
  PV3 <- summary(lm(variable3))$coefficients[v_dep3, 'Pr(>|t|)']
  E3 <- summary(lm(variable3))$coefficients[v_dep3, 'Estimate']
  
  P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, PV1, P.VAL.DICT)
  assign(title_var2, PV2, P.VAL.DICT)
  assign(title_var3, PV3, P.VAL.DICT)
  
  E.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, E1, E.DICT)
  assign(title_var2, E2, E.DICT)
  assign(title_var3, E3, E.DICT)
  
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0.05)
      cat('The variables', i, 'do not show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
} # SUMMARY for 3 variables

SUM.4VAR <- function(variable1, v_dep1, title_var1, variable2, 
                     v_dep2, title_var2, variable3, v_dep3, title_var3, 
                     variable4, v_dep4, title_var4){
  PV1 <- summary(lm(variable1))$coefficients[v_dep1, 'Pr(>|t|)']
  E1 <- summary(lm(variable1))$coefficients[v_dep1, 'Estimate']
  
  PV2 <- summary(lm(variable2))$coefficients[v_dep2, 'Pr(>|t|)']
  E2 <- summary(lm(variable2))$coefficients[v_dep2, 'Estimate']
  
  PV3 <- summary(lm(variable3))$coefficients[v_dep3, 'Pr(>|t|)']
  E3 <- summary(lm(variable3))$coefficients[v_dep3, 'Estimate']
  
  PV4 <- summary(lm(variable4))$coefficients[v_dep4, 'Pr(>|t|)']
  E4 <- summary(lm(variable4))$coefficients[v_dep4, 'Estimate']
  
  P.VAL.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, PV1, P.VAL.DICT)
  assign(title_var2, PV2, P.VAL.DICT)
  assign(title_var3, PV3, P.VAL.DICT)
  assign(title_var4, PV4, P.VAL.DICT)
  
  
  E.DICT <- new.env(hash = T, parent = emptyenv())
  assign(title_var1, E1, E.DICT)
  assign(title_var2, E2, E.DICT)
  assign(title_var3, E3, E.DICT)
  assign(title_var4, E4, E.DICT)
  
  for (i in ls(P.VAL.DICT)) {
    if (P.VAL.DICT[[i]] > 0.05)
      cat('The variables', i, 'do not show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')
    else
      cat('The variables', i, 'show an association with a p value of',  
          P.VAL.DICT[[i]], 'and an estimate of', E.DICT[[i]], '\n')}
} #SUMMARY for 4 variables

}
