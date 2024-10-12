##################################################################
##### COFUST Algorithm by Disegna et al. (2017)
##### Updated: October 8th, 2024
##### Contributors:  A. Benevento, F. Durante
##################################################################

##### General function to install packages if not present
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

##### Define bivariate copula functions
install_and_load("copula")

M.F=function(U){
return(pmin(U[,1],U[,2]))
}

W.F=function(U){
  return(pmax(U[,1]+U[,2]-1,0))
}

Pi.F=function(U){
  return(U[,1]*U[,2])
}

# Convex combination of M and W
MW.F=function(U,s){
  return((1-s)*M.F(U)+s*W.F(U))
}

# Convex combination of M and Pi
MPi.F=function(U,s){
  return((1-s)*M.F(U)+s*Pi.F(U))
}

##################################################################
# COFUST spatio-temporal dissimilarity function
# Input:  X.=time series obs (2 columns), 
#         s=normalized spatial distance
#         L=1,2 (L1 or L2 norm)
#         max.C.sp= "W" or ""Pi" # Spatial copula for maximal distance
#         f.type=1,2 (no exp, exp)
#         alpha=mixing parameter (weight of the spatial component)
##################################################################

COFUST.diss.ts.sp<-function(X.,s=0,max.C.sp="W",L=2,f.type=2,alpha=0){
  u=pobs(X.)
  
  # Precompute common terms
  M_val <- M.F(u)
  Cn_val <- C.n(u, X = u)
  MW_val <- MW.F(u, s)
  MPi_val <- MPi.F(u, s)
  
  # Calculate the values based on the `alpha` weight
  value11 <- max(M_val - ((1 - alpha) * Cn_val + alpha * MW_val))
  value12 <- mean(M_val - ((1 - alpha) * Cn_val + alpha * MW_val))
  value21 <- max(M_val - ((1 - alpha) * Cn_val + alpha * MPi_val))
  value22 <- mean(M_val - ((1 - alpha) * Cn_val + alpha * MPi_val))
  
  # Determine which value to return
  if (max.C.sp=="W") {
    result <- if (L == 1) value11 else value12
  } else {
    result <- if (L == 1) value21 else value22
  }
  
  # Apply transformation based on f.type
  if (f.type == 1) {
    return(result)
  } else {
    return(exp(result) - 1)
  }
}

##################################################################
# COFUST spatio-temporal dissimilarity matrix
# Input: X.=time series obs, S.=normalized spatial distance 
##################################################################
COFUST.Diss.ts.sp=function(X.,S.,max.C.sp= "W",L=2,f.type=2,alpha=0){
  # Use 'pobs' function to transform X.
  U. <- pobs(X.)
  
  # Preallocate the result matrix
  n <- ncol(U.)
  C.temp <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Compute dissimilarity for the pair (i, j) 
      C.temp[i, j] <- COFUST.diss.ts.sp(cbind(U.[, i], U.[, j]), s = S.[i, j], max.C.sp, L, f.type, alpha)
      C.temp[j, i] <- C.temp[i, j]
    }
  }
return(C.temp)
}