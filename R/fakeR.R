# helper function to simulate in level 3 stealth mode 
# (i.e. from a uniform distribution)
.level3_sim <- function(dataset, cluster=NA, n=n, level3.noise=FALSE){
  # for time series data, just look at the cluster variables
  if(!any(is.na(cluster))){
    dataset <- dataset[cluster]}
  mins <- as.numeric(apply(dataset,2,min, na.rm=TRUE))
  maxes <- as.numeric(apply(dataset,2,max, na.rm=TRUE))
  vars <- colnames(dataset)
  num_vars <- length(vars)
  if(level3.noise==TRUE){
    for(i in 1:num_vars){
      sd <- (maxes[i]-mins[i])/4
      mins[i] <- mins[i] + stats::rnorm(1, 0, sd)
      maxes[i] <- maxes[i] + stats::rnorm(1, 0, sd)}
    print('added noise')
  }
  for(i in 1:num_vars){
    dataset[vars[i]] <- stats::runif(n, mins[i], maxes[i])}
  return(dataset)
}

#helper function to get dataset properties
.get_dataset_properties <- function(dataset){
  col <- dim(dataset)[2] # number of columns
  ### variable type check
  num <- rep(NA, col)    # see what's numeric
  ord <- rep(NA, col)    # see what's an ordered factor
  integer <- rep(NA, col)
  # which columns are numeric (the others are factors)?
  for (i in 1:col){
    num[i] <- is.numeric(dataset[,i])
    ord[i] <- is.ordered(dataset[,i])
    integer <- isTRUE(all.equal(dataset[,i], 
                                as.integer(dataset[,i]), 
                                check.attributes = FALSE))
  }
  return(list('col'=col,'num'=num,'ord'=ord,'int'=integer))
}

# helper function to process the time series (i.e. add missing dates and 0 values)
.process_ts <- function(dataset, date, val, date.index){
  #only keep relevant columns in the dataset
  ts_new <- dataset[,(names(dataset) %in% c(date, val)),
                    drop=FALSE]
  
  # make sure value is numeric
  ts_new[,val] <- as.numeric(ts_new[,val])
  # make sure dates are date objects
  if(date.index==TRUE){
    ts_new[,date] <- as.Date(ts_new[,date]) }
  
  # get all dates in the time series range 
  alldates <- seq(min(ts_new[,date]), max(ts_new[,date]), 1)
  
  # filter out the dates already in the dataframe
  dates0 <- alldates[!(alldates %in% ts_new[,date]),drop=FALSE]
  # make all of the values 0
  data0 <- data.frame(a=dates0, b=0)
  colnames(data0) <- c(date, val)
  # combine both parts
  ts_final <- rbind(ts_new, data0)
  # order the column by date
  ts_final <- ts_final[order(ts_final[,date]),]
  return(ts_final)
}

##############################################################
simulate_dataset_ts <- function(dataset, digits=2, n=NA, cluster=NA, time.variable=NA,
                                date.index=FALSE, complete.panel=FALSE,
                                stealth.level=2, level3.noise=FALSE, use.miss=TRUE, ignore=NA){
  
  ## checks
  if((is.data.frame(dataset))==0){
    stop("Data must be a data frame or matrix")}
  if(any(is.na(cluster))){
    stop("Must have at least one time series column.")}
  if(any(is.na(time.variable))){
    stop("Must have a date column for every time series values column.")}
  num_ts <- length(cluster)
  num_date <- length(time.variable)
  dataset <- dataset[c(cluster, time.variable)]
  # check that there is a corresponding date column for every time series values column
  if(num_ts!=num_date){
    stop("Must have a date column for every time series values column.")}
  # check the stealth.level
  if(stealth.level==1){
    stop("No option stealth.level==1. See package MTS to simulate multivariate time series.")}
  else if(stealth.level!=2 & stealth.level!=3){
    stop("Must choose a stealth level of 2 (0 covariance structure between time series) or 3 
         (0 covariance structure between and within time series).")}
  
  # ignore certain columns, if specified
  if(!any(is.na(ignore))){
    ignored_cols <- dataset[ignore]
    num_ignore  <- length(ignore)
    for(i in 1:num_ignore){
      dataset[ignore[i]] <- NULL
    }
  }
  
  # get global dimensions for the data frame
  row <- dim(dataset)[1] # number of rows
  if(is.na(n))(n <- row) # sets unspecified sample size to num rows
  del <- is.na(dataset)
  if(n!=row){
    select <- round(stats::runif(n, 0.5, row+.49),0)
    del <- del[select,]}
  
  print("Some clustered time series data...")
  # create new dataframe for the time series subset of dataset
  fake_ts <- data.frame(matrix(NA, nrow=row, ncol=0))
  
  # if indicated, process each time series by filling in zeros for days of no activity
  if(complete.panel==TRUE){
    processed_df <- rep(NA, num_ts)
    #process each column iteratively
    for(i in 1:num_ts){
      processed_df[[i]] <- .process_ts(dataset, time.variable[i],cluster[i],date.index)
      # add each column to fake_ts
      fake_ts <- do.call(cbind, processed_df)}}
  
  # else, directly move all the time series data over to fake_ts
  else{
    fake_ts <- dataset}
  print("Processing done...")
  
  # for a stealth.level of 2, simulate each time series separately and assume independence
  if(stealth.level==2){
    # use ARIMA to fit each time series and simulate from the ARIMA fit
    for(i in 1:num_ts){
      ts_fit <- stats::arima(fake_ts[,cluster[i]])
      fake_ts[,cluster[i]] <- stats::arima.sim(model=as.list(stats::coef(ts_fit)), 
                                                length(fake_ts[,cluster[i]]))}}

  # at a stealth level of three, simulate time series by sampling from a uniform 
  # distribution of the variable min to variable max, plus Gaussian noise with a standard
  # deviation of 1/4 of the range of the data if specified
  
  else if(stealth.level==3){
    dataset[,cluster] <- .level3_sim(dataset, cluster=cluster, 
                                     n=n, level3.noise=level3.noise)}    
  # round the data to the requested digits
  fake_ts <- round(fake_ts, digits)
  # insert the missing data, if so requested
  if(use.miss==TRUE)(fake_ts[del] <- NA)  
  # reinsert the columns to ignore
  if(!any(is.na(ignore))){
    fake_ts <- cbind(ignored_cols, fake_ts)
  }
  return(fake_ts)
  }
##############################################################
# helper function to sample from a contingency table
.sampct <- function(n, table) {
  # sample with replacement from a multivariate distribution
  # defined by a contingency table
  if(!is.null(dfr <- dimnames(table)) && prod(sapply(dfr, length)) == length(table)){
    dfr <- expand.grid(dfr)}
  else{dfr <- expand.grid(lapply(dim(table), seq))}
  dfr[sample(1:nrow(dfr), n, prob = table, replace = T), ]
}

# helper function to simulate the categorical data from the dataset
.simulate_categorical <- function(dataset, location, row){
  fake_categorical <- data.frame(matrix(NA, nrow=row, ncol=0))
  # for level 1, create a multivariate distribution of all categorical data
  # and simulate from it
  # create a contingency table to sample from
  contingency_table <- table(dataset[location])
  # title the columns the same as the original
  column_names <- colnames(dataset[location])
  num_cat <- length(column_names)
  # add this sampled data to a column in the new dataframe
  fake_categorical[column_names] <- .sampct(row, contingency_table)
  
  for(i in 1:num_cat){
    # give the column vector a type factor and/or character if it began as such
    if(is.factor(dataset[,i])){
      fake_categorical[[column_names[i]]] <- as.factor(fake_categorical[[column_names[i]]])}
    if(is.character(dataset[,i])){
      fake_categorical[[column_names[i]]] <- as.character(fake_categorical[[column_names[i]]])}
  }
  return(fake_categorical)
}

# helper function to simulate the numeric and ordered data
.simulate_num_and_ordered <- function(dataset, row, props, het.suppress, het.ML, mvt.method, use.levels){
  mixedMeans <- rep(0, props$col)
  mixedMeans[props$num] <- apply(dataset[,props$num], 2, mean, na.rm=TRUE)
  # estimate a heterogeneous correlation matrix
  if (het.suppress==TRUE){
    suppressWarnings(het <- polycor::hetcor(dataset, ML=het.ML))} 
  else{het <- hetcor(dataset, ML=het.ML)}
  mixedCov <- het$correlations
  # make a diagonal matrix of standard deviations to turn the 
  # correlation matrix into a covariance matrix
  stand <- matrix(0, props$col, props$col)
  diag(stand) <- rep(1, props$col)
  diag(stand)[props$num] <- apply(dataset[,props$num], 2, sd, na.rm=TRUE)
  # pre and post multiply hetero cor matrix by diagonal sd matrix
  mixedCov <- stand %*% mixedCov %*% stand
  # generate the data
  fake <- as.data.frame(mvtnorm::rmvnorm(row, mixedMeans, mixedCov, mvt.method))
  # turn the appropriate variables into factors
  for (i in (1:props$col)[!props$num]){
    # the original data for this column
    old <- dataset[,i]    
    # the new data for this column, omiting NAs
    new <- fake[!is.na(fake[,i]),i]    
    # the levels of the original factor
    lev <- levels(old)    
    # establish cutpoints in new variable from cdf of old factor
    cut <- cumsum(table(old))/(sum(!is.na(old)))   
    # put continuous variable into a matrix, repeating value across columns
    wide <- matrix(new, length(new), length(lev))    
    # put the cutpoints in a matrix, repeating the cut point values across rows
    crit <- matrix(stats::quantile(new, cut), length(new), length(lev), byrow=TRUE)   
    # for each value (row of the wide matrix), 
    # number of cutpoints the value surpasses is the category
    fake[!is.na(fake[,i]),i] <- apply(wide>crit, 1, sum)    
    # make it a factor
    fake[,i] <- factor(fake[,i], ordered=TRUE)    
    # give the new factor the same levels as the old variable
    if(length(levels(fake[,i]))!=length(lev))message(
      paste("Fewer categories in simulated variable", 
            names(fake)[i], "than in input variable", names(dataset)[i]))
    if(use.levels==TRUE&(length(levels(fake[,i]))==length(lev))){
      levels(fake[,i]) <- lev} 
    else{levels(fake[,i]) <- 1:length(lev)}}
  return(fake)
}
###########################################################
simulate_dataset <- function(dataset, digits=2, n=NA, 
                             use.levels=TRUE, use.miss=TRUE, 
                             mvt.method="eigen", het.ML=FALSE, 
                             het.suppress=TRUE, stealth.level=1,
                             level3.noise=FALSE, ignore=NA){
  
  # requires data frame or matrix
  if((is.data.frame(dataset)+is.matrix(dataset))==0){
    stop("Data must be a data frame or matrix")}
  
  # check the stealth.level
  if(stealth.level!=2 & stealth.level!=3 & stealth.level!=1){
    stop("Must choose a stealth level of 2 (0 covariance structure between variables) or 3 
         (0 covariance structure within variables).")}
  
  # ignore certain columns, if specified
  if(!any(is.na(ignore))){
    ignored_cols <- dataset[ignore]
    num_ignore  <- length(ignore)
    for(i in 1:num_ignore){
      dataset[ignore[i]] <- NULL
    }
  }
  props <- .get_dataset_properties(dataset)
  
  ## get dataset properties that remain constant
  row <- dim(dataset)[1] # number of rows
  del <- is.na(dataset)  # records position of NAs in dataset
  if(is.na(n))(n <- row) # sets unspecified sample size to num rows
  # if n is not equal 
  if(n!=row){
    select <- round(stats::runif(n, 0.5, row+.49),0)
    del <- del[select,]}
  
  if(stealth.level==1){
    # check for unordered factors and characters
    location <- !(props$num|props$ord)
    unorder <- sum(location)
    ### if characters and/or unordered factors, start here
    if(unorder>0){
      print("Some unordered factors...")
      # create new dataframe for the categorical subset of dataset
      fake_categorical <- .simulate_categorical(dataset, location, n)}
    # remove categorical variables from dataset
    dataset[which(location)] <- list(NULL)
    # if there are no more columns to simulate
    if(ncol(dataset)==0){return(fake_categorical)}
    # redetermine dataset properties now that categorical variables removed
    props <- .get_dataset_properties(dataset)
    ### from here, if everything is numeric, start here
    if(sum(!props$num)==0){
      print("Numeric variables. No ordered factors...")
      # draw from a multivariate normal distribution that takes into account covariances
      # generate data with rmvnorm
      fake <- data.frame(mvtnorm::rmvnorm(n, apply(dataset, 2, mean, na.rm=TRUE),
                                          stats::cov(dataset, use="pairwise.complete.obs"),
                                          mvt.method))}
    
    ### if there are ordered factors, skip to here
    else{
      print("Some numeric variables and ordered factors...")
      fake <- .simulate_num_and_ordered(dataset, row, props,het.suppress, 
                                        het.ML, mvt.method, use.levels)}
  }
  # for levels 2 and 3, simulate from each factor column separately 
  # (no distinguishing between ordered and unordered)
  else if(stealth.level==2|stealth.level==3){
    location <- !(props$num)
    fake <- data.frame(matrix(NA, nrow=row, ncol=0))
    for(i in 1:props$col){
      # title the column for the fake data the same as the original
      column_name <- colnames(dataset)[i]
      # for categorical data
      if(location[i]==TRUE){
        if(stealth.level==2){
          # create a contingency table to sample from
          contingency_table <- table(dataset[,i])
          # add this sampled data to a column in the new dataframe
          fake[column_name] <- .sampct(row, contingency_table)}
        else if(stealth.level==3){
          # find all unique values in the column
          col_vals <- unique(dataset[,i])
          # add this sampled data to a column in the new dataframe
          fake[column_name] <- sample(col_vals[!is.na(col_vals)], n, replace=TRUE)}}
      # for numeric data
      else{
        if(stealth.level==2){
          fake[column_name] <- stats::rnorm(n, mean(dataset[,i], na.rm=TRUE), stats::sd(dataset[,i], na.rm=TRUE))}
        else if(stealth.level==3){
          fake <- .level3_sim(dataset, n=n, level3.noise=level3.noise)}}
    }
  }    
  # round the data to the requested digits
  fake[,props$num] <- round(fake[,props$num], digits)    
  # give the noncategorical variables names
  names(fake) <- names(dataset)
  # make integers of the variables that started as such
  fake[,props$int] <- round(fake[,props$int]) 
  # insert the missing data, if so requested
  if(use.miss==TRUE & !length(props$del)==0){
    fake[del] <- NA}   
  # append the categorical data if present 
  if(stealth.level==1){if(unorder>0){fake <- data.frame(fake, fake_categorical)}}
  if(!any(is.na(ignore))){
    fake <- cbind(ignored_cols, fake)
  }
  return(fake)
  }
