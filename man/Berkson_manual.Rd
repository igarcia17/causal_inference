\name{df_pval_estimates}
\alias{Df and summary for Berkson's paradox}
\title{Berkson's analysis}

\description{This function recieves 3 variables and creates a dataframe with those variables. It retrieves through a summary function (gml, as the variables are categorical) the pvalues and estimates of the variables when adjusting for the collider and when not adjusting. These values are printed in a dataframe form where it shows the models and the estimates and pvalues of each one. It also returns a string for each model analysis with the association or lack thereof and the values of the estimates and significance.

}
\usage{
df_pval_estimates(var1, var2, var3)
}

\arguments{
  \item{var1}
  % First variable for the regression model.
  \item{var2}
  % Second variable for the regression model.
  \item{var3}
  % Variable to adjust model.
}

\value{Returns a dataframe with pvalues and estimates. Also returns a string describing whether there is an association between the variables of study.}

\author{Nermina Logo Lendo
}

\examples{
##---- Should NOT be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

{df_pval_estimates <- function(var1, var2, var3) {
  
  dataset <- data.frame(var1,var2,var3)
  
  est_v1_v2 <- summary(glm(var1~var2, data = dataset, family = 'binomial')
  )$coefficients['var2', 'Estimate']
  pval_v1_v2 <- summary(glm(var1~var2, data = dataset, family = 'binomial')
  )$coefficients['var2', 'Pr(>|z|)']
  
  
  est_v1_v2_v3 <- summary(glm(var1~var2+var3, data=dataset, 
                              family = 'binomial'))$coefficients['var2', 'Estimate']
  pval_v1_v2_v3 <-summary(glm(var1~var2+var3, data=dataset, 
                              family = 'binomial'))$coefficients['var2', 'Pr(>|z|)']
  
  df_final <- data.frame('values_of' = c('Estimates', 'P_values'), 'diab_chol' = 
                           c(est_v1_v2, pval_v1_v2), 'diab_chol_HP' = 
                           c(est_v1_v2_v3, pval_v1_v2_v3))
  
  print(df_final)
  cat('\n')
  
  {if (pval_v1_v2 > 0.05)
    cat('The variables diabetes and cholecystitis do not show an association. The p value is',
        pval_v1_v2, 'and the estimate is', est_v1_v2, 'for this model', '\n')
    else
      cat('The variables diabetes and cholecystitis do show an association. The p value is',
          pval_v1_v2, 'and the estimate is', est_v1_v2, 'for this model', '\n')}
  
  {if (pval_v1_v2_v3 > 0.05)
    cat('The variables diabetes and cholycistitis do not show an association when adjusting by hospital patients. The p value is',
        pval_v1_v2_v3, 'and the estimate is',est_v1_v2_v3, 'for this model', '\n')
    else
      cat('The variables diabetes and cholycistitis do show an association when adjusting by hospital patients. The p value is',
          pval_v1_v2_v3, 'and the estimate is',est_v1_v2_v3, 'for this model', '\n')}
}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
