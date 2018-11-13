run_vector_change <- function(y.mat, covs, x.mat.red, x.mat.full, permute=999){
  #
  # This script will run an analysis to test two-state phenotypic changes across groups.
  # It will estimate the necessary parameters and run a permutation procedure to determine statistical significance
  # Code was modified from 
  # Collyer & Adams (2007). Analysis of two‐state multivariate phenotypic change in ecological studies. Ecology, 88(3), 683-692.
  # Data stored in figshare: https://doi.org/10.6084/m9.figshare.c.3299645.v1
  #
  # Usage:
  #   y.mat = Y matrix to solve regression
  #   covs  = Covariates to control for during the parameter estimation, should contain the intercept as well.
  #           Make sure that all numeric covariates are scaled, so zero is the mean 
  #   x.mat.red  = X matrix to solve regression, without the interaction term
  #   x.mat.full = X matrix to solve regression, with interaction term of interest in the last columns.
  #                For both X matrices, make sure that the first column contains the two state condition, while groups are next
  #   permute = number of permutations to determine statistical significance (default 999)
  
  # 1. Setting data as matrices and setting general variables
  y.mat      <- as.matrix(y.mat)
  covs       <- as.matrix(covs)
  ncovs      <- ncol(covs) #number of covariates
  ngroups    <- ncol(x.mat.red) #number of groups
  x.mat.red  <- cbind(covs, as.matrix(x.mat.red)) # Full X matrix
  x.mat.full <- cbind(covs, as.matrix(x.mat.full)) # Reduced X matrix
  
  # 2. Estimate parameters for full and reduced matrices
  b.mat.full <- solve((t(x.mat.full) %*% x.mat.full)) %*% (t(x.mat.full) %*% y.mat)
  b.mat.red  <- solve((t(x.mat.red) %*% x.mat.red)) %*% (t(x.mat.red) %*% y.mat)
  
  # 3. Estimate means
  #Arbitrarily, a = groups (might be more than two), and b = two state condition (like sex)
  as_list  <- paste("a", 1:ngroups, sep="")
  bs_list  <- paste("b", 1:2, sep="")
  all_comb <- as.vector(outer(as_list, bs_list, function(x, y) paste(x, y, sep = ''))) # setting names of groups by state conditions
  
  base <- c(1, rep(0,ncovs-1)) # Base contains the intercept and zeros equal to the number of covariates 
  base <- t(replicate(length(all_comb), base))
  
  sex_groups <- expand.grid(rep(list(0:1), ngroups)) # Generating all possible combinations of state by group combinations
  keep       <- !(apply(sex_groups[,1:(ncol(sex_groups)-1)], 1, function(x) sum(x)) > 1)
  sex_groups <- sex_groups[keep,]
  sex_groups <- as.matrix(cbind(sex_groups[,ncol(sex_groups)], sex_groups[,-ncol(sex_groups)]))
  
  for(i in 2:ngroups){
    sex_groups <- cbind(sex_groups, sex_groups[,1] * sex_groups[,i]) # Generating the interaction term
  }
  
  x.ls.full  <- cbind(base, sex_groups)
  colnames(x.ls.full) <- NULL
  rownames(x.ls.full) <- all_comb
  
  x.ls.red <- x.ls.full[,1:(ncol(x.ls.full) - (ngroups-1))]
  
  obs.ls.full <- x.ls.full %*% b.mat.full  # Observed ls means (full model)
  obs.ls.red  <- x.ls.red %*% b.mat.red     # Observed ls means (reduced model)
  
  # 4. Vector and statistics calculations 
  #Vectors
  obs.vect <- matrix(0, nrow = ngroups , ncol = ncol(obs.ls.full)) 
  for (i in 1:ngroups){
    obs.vect[i,] <- obs.ls.full[i,] - obs.ls.full[i+ngroups,] # These are the phenotypic change vectors
  }
  
  #Lengths
  obs.d <- matrix(0, nrow = ngroups, ncol = 1)
  for (i in 1:ngroups){
    obs.d[i,] <- sqrt(t(obs.vect[i,]) %*% obs.vect[i,]) # These are lengths of vectors
  }
  
  #Contrast
  obs.contrast <- matrix(0, nrow = ngroups * (ngroups-1) / 2, ncol = 1)
  row = 1
  for (i in 1:ngroups){
    t = i + 1
    while (t <= ngroups){
      obs.contrast[row,1] <- abs(obs.d[i,1] - obs.d[t,1] )
      t = t + 1
      row = row + 1
    }
  }

  #Angles
  obs.angle <- matrix(0, nrow = ngroups * (ngroups-1) / 2, ncol = 1)
  row = 1
  for (i in 1:ngroups){
    t = i + 1
    while (t <= ngroups){
      obs.angle[row,1] <- acos(t((obs.vect[i,]) / obs.d[i,1]) %*% ((obs.vect[t,])/obs.d[t,1]))
      obs.angle[row,1] <- obs.angle[row,1]*180/pi # This step is only necessary to convert radians to degrees
      t = t + 1
      row = row + 1
    }
  }
  
  # 5. Set-up permutation procedure
  y.hat <- x.mat.red %*% b.mat.red     # Predicted values from reduced model
  y.res <- y.mat - y.hat               # Resdiuals of reduced mode (these are the permuted units)
  
  # PERMUTATION PROCEDURE
  # Need to set-up distributions to be generated
  dist.d <- matrix(0, nrow = (permute + 1), ncol = ngroups)
  dist.d[1,] <- obs.d # Observed values are first random values
  
  dist.contrast <- matrix(0, nrow = (permute + 1), ncol = ngroups * (ngroups-1) / 2)
  dist.contrast[1,] <- obs.contrast # Observed values are first random values
  
  dist.angle <- matrix(0, nrow = (permute + 1), ncol = ngroups * (ngroups-1) / 2)
  dist.angle[1,] <- obs.angle # Observed values are first random values
  
  # In addition to saving random values, it is wise to save the outcome of comparisons
  # of observed and random values.
  # This can be done with logical statements (below).
  # Separate distributions are created for these comparisons.
  # The 'p' indicates that these distributions will be used 
  # to calculate empirical probabilities.
  
  pdist.contrast <- matrix(0, nrow = (permute + 1), ncol = ngroups * (ngroups-1) / 2)
  pdist.contrast[1,] <- 1 
  
  pdist.angle <- matrix(0, nrow = (permute + 1), ncol = ngroups * (ngroups-1) / 2)
  pdist.angle[1,] <- 1
  
  # Create an array from 1 to number of object in data set
  # This will be randomized later
  line <- array(1:(length(x.mat.full[,1])),dim=c(length(x.mat.full[,1])))
  
  for (i in 1:permute){
    
    line.rand <- sample(line, replace=FALSE)
    y.res.temp <- cbind(line.rand, y.res)
    z <- (order(line.rand))
    y.res.temp2 <- as.matrix(y.res.temp[z,])
    y.res.rand  <- y.res.temp2[,-1]  # Rows of residuals are now randomized
    
    # Create random values
    y.rand <- y.hat + y.res.rand
    
    # Estimate parameters
    b.mat.rand <- solve((t(x.mat.full)%*%x.mat.full))%*%(t(x.mat.full)%*%y.rand)
    
    # Calculate LS means
    rand.ls.full <- x.ls.full%*%b.mat.rand
    
    # Repeat fourth step for random data!
    rand.vect <- matrix(0, nrow = ngroups , ncol = ncol(rand.ls.full)) 
    for (l in 1:ngroups){
      rand.vect[l,] <- rand.ls.full[l,] - rand.ls.full[l+ngroups,] # These are the phenotypic change vectors
    }
    
    #Lengths
    rand.d <- matrix(0, nrow = ngroups, ncol = 1)
    for (l in 1:ngroups){
      rand.d[l,] <- sqrt(t(rand.vect[l,]) %*% rand.vect[l,]) # These are lengths of vectors
    }
    
    #Contrast
    rand.contrast <- matrix(0, nrow = ngroups * (ngroups-1) / 2, ncol = 1)
    row = 1
    for (l in 1:ngroups){
      t = l + 1
      while (t <= ngroups){
        rand.contrast[row,1] <- abs(rand.d[l,1] - rand.d[t,1] )
        t = t + 1
        row = row + 1
      }
    }
    
    #Angles
    rand.angle <- matrix(0, nrow = ngroups * (ngroups-1) / 2, ncol = 1)
    row = 1
    for (l in 1:ngroups){
      t = l + 1
      while (t <= ngroups){
        rand.angle[row,1] <- acos(t((rand.vect[l,]) / rand.d[l,1]) %*% ((rand.vect[t,])/rand.d[t,1]))
        rand.angle[row,1] <- rand.angle[row,1]*180/pi # This step is only necessary to convert radians to degrees
        t = t + 1
        row = row + 1
      }
    }
    
    # Append distributions
    
    dist.d[i+1, ] <- rand.d
    dist.contrast[i+1,] <- rand.contrast
    dist.angle[i+1,]    <- rand.angle
    
    aa <- ifelse(rand.contrast>=obs.contrast,1,0)
    bb <- ifelse(rand.angle>=obs.angle,1,0)
    
    pdist.contrast[i+1,] <- aa
    pdist.angle[i+1,]    <- bb
    
  }
  
  # Empirical probabilities are calculated as follows
  
  p.contrast <- apply(pdist.contrast, 2, function(x) sum(x)/(permute+1)) 
  p.angle    <- apply(pdist.angle, 2, function(x) sum(x)/(permute+1)) 
  
  # OUTPUT
  #Print results to console
  #for (i in 1:ngroups){
  #  res <- obs.d[i,1]
  #  cat(sprintf("The length of vector %1.0f is: %f\n\n", i, res))
  #}
  #for (i in 1:length(obs.contrast)){
  #  res  <- obs.contrast[i,1]
  #  res2 <- p.contrast[i] 
  #  cat(sprintf("The absolute difference between two lengths is, |d1 - d2| = %f\n With a p-value = %f\n\n", res, res2))
  #}
  #for (i in 1:length(obs.angle)){
  #  res  <- obs.angle[i,1]
  #  res2 <- p.angle[i]
  #  cat(sprintf("The angle between vectors is %f\n With a p-value = %f\n\n", res, res2))
  #}
  
  # The next step should produce a table with every random distance, contrast, and angle.
  # Probabilities can be calculated manually from this table, or as above in R.
  
  ds    <- paste("d", 1:ncol(dist.d), sep = "")
  diffs <- paste("abs_diff", 1:ncol(dist.contrast), sep = "")
  angs  <- paste("angle", 1:ncol(dist.angle), sep = "")
  out.head <- c(ds, diffs, angs)
  iter <- array(1:permute)
  iter.lab <- "observed"
  iter.lab <- append(iter.lab,iter)
  out.mat  <- cbind(dist.d, dist.contrast, dist.angle)
  colnames(out.mat) <- out.head
  rownames(out.mat) <- iter.lab
  return(out.mat)
}

pairwise_pvals <- function(out.mat, angle = 1){
  #
  # Returns a table with the pairwise p-values from the run_vector_change functions
  #
  # Usage:
  #   out.mat: output matrix from run_vector_change containing the vector lengths, their pairwise difference, and the anlge difference
  #   angle: whether to return the angle pairwise table (default = 1). If angle = 0, then returns difference table
  #
  
  column_names <- attr(out.mat, "dimnames")[[2]] # Column names
  rows         <- nrow(out.mat) # Number of rows
  k            <- length(column_names[grep("^[d].*", column_names)]) # Number of groups
  n_comparisons <- k * (k-1) / 2  # Number of comparisons
  
  diff_table  <- matrix(0, nrow = k, ncol = k)
  angle_table <- matrix(0, nrow = k, ncol = k)
  pvals       <- matrix(0, nrow = 1, ncol = n_comparisons)
  col = k + 1
  
  for (i in 1:n_comparisons){
    pval <- sum(out.mat[1,col] < out.mat[2:rows,col])/rows
    col  <- col +1
    pvals[1,i] <- pval
  }
  
  diff_table[lower.tri(diff_table)] <- pvals
  diff_table <- as.matrix(Matrix::forceSymmetric(diff_table, uplo="L"))
  
  pvals <- matrix(0, nrow = 1, ncol = n_comparisons)
  for (i in 1:n_comparisons){
    pval <- sum(out.mat[1,col] < out.mat[2:rows,col])/rows
    col  <- col +1
    pvals[1,i] <- pval
  }
  
  angle_table[lower.tri(angle_table)] <- pvals
  angle_table <- as.matrix(Matrix::forceSymmetric(angle_table, uplo="L"))
  
  if (angle == 1){
    return(angle_table)
  } else if (angle == 0){
    return(diff_table)
  }
  
}

get_sex_vectors <- function(y.mat, covs, x.mat.full){
  #
  # Estimate the means of the effects of sex by population.
  #
  # Usage:
  #   y.mat = Y matrix to solve regression
  #   covs  = Covariates to control for during the parameter estimation, should contain the intercept as well.
  #           Make sure that all numeric covariates are scaled, so zero is the mean 
  #   x.mat.full = X matrix to solve regression, with interaction term of interest in the last columns.
  #                For both X matrices, make sure that the first column contains the two state condition, while groups are next
  #
  
  # 1. Setting data as matrices and setting general variables
  y.mat      <- as.matrix(y.mat)
  covs       <- as.matrix(covs)
  ncovs      <- ncol(covs) #number of covariates
  ngroups    <- ( ( ncol(x.mat.full) - 1 ) / 2 ) + 1 #number of groups
  x.mat.full <- cbind(covs, as.matrix(x.mat.full)) # Reduced X matrix
  
  # 2. Estimate parameters for full and reduced matrices
  b.mat.full <- solve((t(x.mat.full) %*% x.mat.full)) %*% (t(x.mat.full) %*% y.mat)
  
  # 3. Estimate means
  #Arbitrarily, a = groups (might be more than two), and b = two state condition (like sex)
  as_list  <- paste("a", 1:ngroups, sep="")
  bs_list  <- paste("b", 1:2, sep="")
  all_comb <- as.vector(outer(as_list, bs_list, function(x, y) paste(x, y, sep = ''))) # setting names of groups by state conditions
  
  base <- c(1, rep(0,ncovs-1)) # Base contains the intercept and zeros equal to the number of covariates 
  base <- t(replicate(length(all_comb), base))
  
  sex_groups <- expand.grid(rep(list(0:1), ngroups)) # Generating all possible combinations of state by group combinations
  keep       <- !(apply(sex_groups[,1:(ncol(sex_groups)-1)], 1, function(x) sum(x)) > 1)
  sex_groups <- sex_groups[keep,]
  sex_groups <- as.matrix(cbind(sex_groups[,ncol(sex_groups)], sex_groups[,-ncol(sex_groups)]))
  
  for(i in 2:ngroups){
    sex_groups <- cbind(sex_groups, sex_groups[,1] * sex_groups[,i]) # Generating the interaction term
  }
  
  x.ls.full  <- cbind(base, sex_groups)
  colnames(x.ls.full) <- NULL
  rownames(x.ls.full) <- all_comb
  obs.ls.full <- x.ls.full %*% b.mat.full  # Observed ls means (full model)
  sex   <- c(rep(1,ngroups), rep(2,ngroups))
  group <- rep(1:ngroups, 2)
  
  obs.ls.full <- as.data.frame(obs.ls.full)
  obs.ls.full$Sex <- sex
  obs.ls.full$Group <- group
  obs.ls.full$k <- ngroups
  
  return(obs.ls.full)
}