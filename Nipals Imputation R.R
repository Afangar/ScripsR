#Apply impute.nipals to impute NaNs in the database with different variables
library(mixOmics)

# Create a database where the imputed columns will be pasted
# Matriximputada is the new one and Matrixnipals is the database with the data to be imputed
Matriximputada <- as.data.frame(matrix(NA, nrow = nrow(Matrixnipals), ncol = ncol(Matrixnipals)))
colnames(Matriximputada) <- colnames(Matrixnipals)
row.names(Matriximputada) <- row.names(Matrixnipals)

# Total number of columns in the dataframe
total_columns <- ncol(Matrixnipals)

# Perform the imputation in blocks of 3 columns, in this case it is because there are 3 replicates of each experiment
for (i in seq(1, total_columns, by = 3)) {
  current_columns <- Matrixnipals[, i:min(i+2, total_columns)]  #select columns by 3
  filas_NaN <- apply(current_columns, 1, function(x) all(is.na(x))) #check which of the 3 rows are only NAN
  current_columns[filas_NaN, ] <- 0 #rows that are only NaN are changed to zeros
  imputed <- impute.nipals(current_columns, ncomp = 2) #carry out the imputation of the NaN
  Matriximputada[, i:min(i+2, total_columns)] <- imputed # put the 3 columns into the created database
}

#USE this code also in case you do not want it to be imputed if there were 2 NaN and only 1 data among the 3
#columns, what it does is return that data to 0 before imputing, to ignore false positives.
for (i in seq(1, total_columns, by = 3)) {
  current_columns <- Matrixnipals[, i:min(i+2, total_columns)]  
  
  # Return 0 if the other 2 replicas are NaN
  for (j in 1:nrow(current_columns)) {
    for (k in 1:ncol(current_columns)) {
      if (!is.na(current_columns[j, k]) && sum(is.na(current_columns[j, -k])) == (ncol(current_columns) - 1)) {
        current_columns[j, k] <- 0
      }
    }
  }
  
  # Identify rows where all replicas are NaN
  filas_NaN <- apply(current_columns, 1, function(x) all(is.na(x))) 
  
  # Return 0 to completely NaN rows to be able to impute
  current_columns[filas_NaN, ] <- 0 
  
  # Impute values using nipals
  imputed <- impute.nipals(current_columns, ncomp = 2)$completeObs
  
  # Save the imputed values in the resulting matrix
  Matriximputada[, i:min(i+2, total_columns)] <- imputed
}