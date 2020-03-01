\name{simulate_dataset_ts}
\alias{simulate_dataset_ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate a dataframe of time series data
}
\description{
This function simulates clustered numeric time series data from an ARIMA model fit.
}
\usage{
simulate_dataset_ts(dataset, digits=2, n=NA, cluster=NA, time.variable=NA,
                    date.index=FALSE, complete.panel=FALSE,
                    stealth.level=2, level3.noise=FALSE, use.miss=TRUE, ignore=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
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
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
the column names of the time series variables. Argument should be in the form of a list if multiple values.
}
  \item{time.variable}{
%%     ~~Describe \code{time.variable} here~~
the column name(s) of the time variables corresponding to each time series variable. Should be the same length as cluster, even if that means including the same time variable multiple times.
}
  \item{date.index}{
%%     ~~Describe \code{date.index} here~~
whether the time variable is a date and should be treated as a Date object.
}
  \item{complete.panel}{
%%     ~~Describe \code{complete.panel} here~~
when set to TRUE, indicates a preprocessing step needed to complete 
the time series columns. Specifically, inserts all missing dates and 
zero values at each of those time points.
}
  \item{stealth.level}{
  when set to 2 (default), simulates independent time series observations. When set to 3, does not take into account any covariances between time points and instead randomly samples from a uniform distribution ranging from the min to the max of the data for each variable. No option 1. 
}
  \item{level3.noise}{
  when set to TRUE, add Gaussian noise to the min and max parameter for the uniform distribution in stealth.level 3. The noise term has a variance of one fourth of the range of the data for any particular variable.  
}
  \item{use.miss}{
%%     ~~Describe \code{use.miss} here~~
when set to TRUE, inserts the missing data like is present in the original.
}
  \item{ignore}{
  specifies which columns to ignore (i.e. to leave as is instead of simulate). Takes in a list of column names as input. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Note that this function is specific to two types of numeric time series, stationary ones and zero-inflated count ones. For modeling clustered numeric data assuming a multivariate normal distribution, look at simulate_dataset. Note that this function only accepts numeric observation types. 
  
The function assumes each time series process is independent of the
others, and allows for a different time variable to be associated with each series. Thus, there is no stealth level of 1 for simulate_dataset_ts(), as this function does not simulate multivariate time series. Columns not part of time series values or time indices are ignored and not simulated.
  
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame. Columns alternating time variable and cluster variable, with the cluster/time variable pairs in the order inputted into the function arguments.
}

\author{
%%  ~~who you are~~
Lily Zhang
Dustin Tingley
}

\examples{
# An example using the treering dataset from the R datasets package
tree_ring <- data.frame(treering)
tree_ring$year <- c(1: nrow(tree_ring))
sim_tree_ring <- simulate_dataset_ts(tree_ring, 
                                     cluster="treering", 
                                     time.variable="year")
par(mfrow = c(2, 1), mar = c(3, 3, 4, 2), mgp = 0.9 * 2:0)
plot (tree_ring$year, tree_ring$treering, type='l', 
      main=paste("Original","Normalized ring width"),
      ylab="Ring width", xlab="Year index")
plot (tree_ring$year, tree_ring$treering, type='l', 
      main=paste("Simulated","Normalized ring width"),
      ylab="Ring width", xlab="Year index")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ datagen }% __ONLY ONE__ keyword per line
