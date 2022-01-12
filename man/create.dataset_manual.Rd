\name{create.dataset}
\alias{create.datasetv2, create.datasetv3}
\title{
Creation of dataset for causal inference examples%%  ~~function to do ... ~~
}
\description{
Create a dataset that at least contains variables X, Y and Z with a given numerical relation between them. Necessary for common cause description functions. Version 2 and 3 required to make dataset with variables A and B respectively, for scenario 1 with ancestor of common cause and scenario 2 when Y is a collider.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Dataset creation(x)
}
\arguments{ -- Nera: falta o quitar
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}

\value{ -- Falta
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Ines Garcia Ortiz
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
b_yz = 6
b_ax = 2
b_by = 1
b_bz = 3

## The function is currently defined as
create.dataset <- function(b_yz, N = 500, b_xy = 3, b_xz = 3,
                           e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(X = runif(N, 1, 100) + rnorm(N, sd = e_x))
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}

create.datasetv2 <- function(b_ax, b_yz=0, N = 500, b_xy = 3, b_xz = 3,
                             e_x = 1, e_y = 1, e_z = 1) {
  name_df <- data.frame(A = runif(N, 1, 100) + rnorm(N))
  name_df$X <- name_df$A * b_ax + rnorm(N, sd = e_x)
  name_df$Y <- name_df$X * b_xy + rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + rnorm(N, sd = e_z)
  return(name_df)
}

create.datasetv3 <- function(b_by, b_bz, b_yz=(-3), N = 500, b_xy = 3, b_xz = 3,
                             e_x = 1, e_y = 1, e_z = 1, e_b = 1) {
  name_df <- data.frame(B = runif(N, 1, 100) + rnorm(N, sd = e_b))
  name_df$X <- runif(N, 1, 100) + rnorm(N, sd = e_x)
  name_df$Y <- name_df$X * b_xy + name_df$B * b_by +rnorm(N, sd = e_y)
  name_df$Z <- name_df$X * b_xz + name_df$Y * b_yz + 
    name_df$B * b_bz+ rnorm(N, sd = e_z)
  return(name_df)
}
create.dataset(b_yz)
create.datasetv2(b_ax)
create.datasetv3(b_by, b_bz)
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
