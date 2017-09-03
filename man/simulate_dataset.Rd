\name{simulate_dataset}
\alias{simulate_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate from a data frame of time-independent data.

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~ 
This function takes as argument an existing dataset in the form of a data frame and outputs a randomized version of all its columns. The function accepts the following types: character variables, numeric variables, and ordered and unordered factor variables.
}
\usage{
simulate_dataset(dataset, digits=2, n=NA, 
                             use.levels=TRUE, use.miss=TRUE, 
                             mvt.method="eigen", het.ML=FALSE, 
                             het.suppress=TRUE, stealth.level=1,
                             level3.noise=FALSE, ignore=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
the data frame from which to generate a randomized version
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
the number of digits after the decimal point to include in the new values
}
  \item{n}{
%%     ~~Describe \code{n} here~~
number of rows in the new data frame. Equal to the number of rows in the original if set to NA, the default.
}
  \item{use.levels}{
%%     ~~Describe \code{use.levels} here~~
when set to true, gives the simulated factor variables the same number of levels as the original.
}
  \item{use.miss}{
%%     ~~Describe \code{use.miss} here~~
when set to TRUE, inserts the missing data like is present in the original (i.e. based on the distribution of missingness in the original data).
}
  \item{mvt.method}{
%%     ~~Describe \code{mvt.method} here~~
specifies the matrix decomposition to be used in sampling from the 
multivariate normal.
}
  \item{het.ML}{
%%     ~~Describe \code{het.ML} here~~
as per the hetcor function, if TRUE, compute maximum-likelihood estimates;if FALSE, compute quick two-step estimates in computing the heterogeneous correlation matrix.
}
  \item{het.suppress}{
%%     ~~Describe \code{het.suppress} here~~
when set to TRUE, suppresses stops from the het.corr function.
}
  \item{stealth.level}{
  when set to 1 (deafult), takes into account the covariances between all the unordered factors and the covariances between the numeric and ordered factors. When set to 2, simulates each variable independently. When set to 3, does not take into account any covariances and instead randomly samples from a uniform distribution ranging from the min to the max of the data for each variable. 
}
  \item{level3.noise}{
  when set to TRUE, add Gaussian noise to the min and max parameter for the uniform distribution in stealth.level 3. The noise term has a variance of one fourth of the range of the data for any particular variable.  
}
  \item{ignore}{
  specifies which columns to ignore (i.e. to leave as is instead of simulate). Takes in a list of column names as input. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function does not account for clustered time series data (see 
simulate_dataset_ts).
  
This function randomly samples each each character and factor variable from the population distribution given in the original dataset. It simulates numeric and ordered factors from a multivariate normal distribution. When both numeric and ordered factors are included, a heterogeneous correlation matrix is used, coercing the means of the ordered factor variables to be 0. 
  
The function only accounts for between-column correlations for numeric and ordered factor variables. Each unordered factor and character column is treated as independent.

The order of the columns in the simulated dataset may differ from the order of the original dataset since the function puts the numeric and ordered factor data in the front and the character and unordered factor data afterwards. The column names stay consistent, however.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a data frame with the same number of columns and same type for each.
}
\references{
%% ~put references to the literature/web site here ~
Inspired by the fakeR function originally created by Ryne Estabrook.
}
\author{
%%  ~~who you are~~
Lily Zhang
Dustin Tingley
}

\examples{
# single column of an unordered, string factor
state_df <- data.frame(division=state.division)
# character variable
state_df$division <- as.character(state_df$division)
# numeric variable
state_df$area <- state.area
# factor variable
state_df$region <- state.region
state_sim <- simulate_dataset(state_df)
}

\keyword{ manip }
\keyword{ datagen }
