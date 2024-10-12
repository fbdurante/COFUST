##################################################################
### FUZZY k-MEDOIDS  
##################################################################
## Inputs
# Y: Matrix or data.frame (rows=objects to cluster)
# m: Parameter of fuzziness (usually, 1.5)
# k: number of clusters
# dist.matrix: dissimilarity matrix
# RS:	number of (random) starts (default: 100)
# max.iter: Maximum number of iterations (default: 100)
## Output
#hard_clustering: data frame with cluster assignments and related memberships  
#U: Membership degree matrix
#H: Prototype matrix, i.e. values (from input matrix Y) assumed by the medoids
#medoid: Vector containing the indexes of the medoid objects
#medoid.label: Vector containing the labels of the medoid objects
#dist2medoids: distance of each object to each medoids
#value: Vector containing the loss function values for the RS starts
#opt.val: optimal value
#membership: maximal membership degrees 
#clustering: crisp cluster assignment for each object
#distance: dissimilarity matrix
##################################################################
Fcmdd <- function(Y, m, k, dist.matrix, RS = 100, max.iter = 100){
  n <- nrow(Y)
  s <- ncol(Y)
  
  # Check row names and set if NULL
  if (is.null(row.names(Y))) row.names(Y) <- 1:n
  
  X <- as.matrix(Y)
  dist <- as.matrix(dist.matrix)
  rownames(dist) <- rownames(Y)
  colnames(dist) <- rownames(Y)
  
  # Initialize optimal objective function value and convergence criteria
  func.opt <- Inf
  conv <- 10^(-6)
  Val.opt <- numeric(RS)
  
  # Iterate RS times to reduce local optima risks
  for (rs in 1:RS) {
    U <- matrix(0, nrow = n, ncol = k)  # Initialize membership matrix
    medoid <- sort(sample(1:n, k, replace = FALSE))  # Random medoid selection
    medoid.old <- medoid + 1  # To ensure while loop runs at least once
    iter <- 1
    
    while (!identical(medoid, medoid.old) && iter <= max.iter) {
      iter <- iter + 1
      medoid.old <- medoid
      H <- X[medoid, ]  # Current medoid vectors
      
      # Compute distance matrix between data points and medoids
      dist.ic <- dist[, medoid]
      
      # Update membership matrix U
      for (i in 1:n) {
        dist_i <- dist.ic[i, ]  # Distance to each medoid
        identical_check <- apply(H, 1, function(h) identical(X[i, ], h))
        
        if (any(identical_check)) {
          cl <- which(identical_check)  # Find which medoid it is identical to
          U[i, ] <- 0
          U[i, cl] <- 1
        } else {
          dist_i[dist_i == 0 | is.na(dist_i)] <- 10^(-9)  # Avoid division by 0 or NA
          denom <- sum((1 / dist_i)^(1 / (m - 1)))  # Sum denominator once
          U[i, ] <- (1 / dist_i^(1 / (m - 1))) / denom  # Update membership degrees
        }
      }
      
      # Update medoids
      wdist <- dist %*% (U^m)
      for (c in 1:k) {
        medoid[c] <- which.min(wdist[, c])  # Choose the point with minimum weighted distance
      }
      medoid <- sort(medoid)  # Sort the medoids after the update
    } # End of while loop
    
    # Compute objective function
    func <- sum((U^m) * dist.ic)
    Val.opt[rs] <- func
    
    # Update optimal values
    if (func < func.opt) {
      U.opt <- U
      H.opt <- H
      func.opt <- func
      medoid.opt <- medoid
      dist.opt <- dist.ic
    }
  }
  
  # Assign row names and prepare output
  row.names(U.opt) <- row.names(Y)
  colnames(U.opt) <- paste("medoid", 1:k, sep = "")
  
  clustering <- apply(U.opt, 1, which.max)
  membership <- apply(U.opt, 1, max)
  
  hard_clustering <- data.frame(clustering, membership)
  rownames(hard_clustering) <- rownames(X)
  
  results <- list(
    U = U.opt,
    H = H.opt,
    medoid = medoid.opt,
    medoid.label = row.names(X)[medoid.opt],
    dist2medoids = dist.opt,
    value = Val.opt,
    opt.val = func.opt,
    membership = membership,
    hard_clustering = hard_clustering,
    clustering = clustering,
    distance = dist
  )
  
  return(results)
}


##################################################################
### FUZZY k-MEDOIDS: Graphical representation of spatial objects  
##################################################################
## Inputs
# S: Geographic coordinates
# Clust: Output of the Fcmdd algorithm
# map: file name "map.geo.json" of the underlying map
## Output
# Graphical representation in ggplot
##################################################################
library(ggplot2)
library(viridis)  

Graph.fcmdd <- function(S, Clust, map) {
  
  # Extract medoids and their clustering and membership information
  med = Clust$medoid
  Medoids = data.frame(
    Longitude = S[med, 1],  # Assuming longitude is in 1st column
    Latitude  = S[med, 2],  # Assuming latitude is in 2nd column
    group = Clust$hard_clustering$clustering[med],
    membership = Clust$hard_clustering$membership[med])
  
  # We relabel the clusters so that label=1 corresponds to
  # the cluster whose medoid is closest to north.
  old_labels = order(Medoids$Latitude)
  new_labels = 1:length(Medoids$Latitude)
  Clust$hard_clustering$clustering = new_labels[match(Clust$hard_clustering$clustering, old_labels)]
  
  # Create the plot using ggplot
  imp_map <- st_read(map)  ### Import geojson map
  
  gp = ggplot(data = imp_map) +
    geom_sf(fill = "whitesmoke", color = "black")+
    geom_point(data = data.frame(S, 
                                 group = Clust$hard_clustering$clustering,
                                 membership = Clust$hard_clustering$membership),
               aes(x = S[, 1], y = S[, 2], 
                   colour = factor(group),size = membership),
               alpha = 0.6) +  # Add transparency to points
    geom_point(data = Medoids[order(Medoids$Latitude),],
               aes(x = Longitude, y = Latitude),
               shape = 3, size = 8, colour = "red") +  # Medoids with a distinct shape and color
    scale_colour_viridis_d(option = "plasma", end = 0.9) +  # Use a colorblind-friendly viridis palette
    labs(colour = "Cluster", size = "Membership") +
    labs(x = "Longitude", y = "Latitude") +  
    theme_light() +
    coord_sf(xlim = c(min(Coord[, 1]) - 0.5, max(Coord[, 1]) + 0.5),  # Set limits based on the coordinates
             ylim = c(min(Coord[, 2]) - 0.5, max(Coord[, 2]) + 0.5))+
    guides(colour = guide_legend(order = 1), 
           size = guide_legend(order = 2)) 
  # Return the generated plot
  return(gp)
}

##################################################################
### FUZZY SILHOUETTE
### See Campello, R.; Hruschka, E.; 
### A fuzzy extension of the silhouette width criterion for cluster analysis. 
### Fuzzy Sets and Systems, v. 157, n. 21, p. 2858
##################################################################
# Input:
# dist: a dissimilarity matrix
# U: the membership matrix
# param: weighting parameter (default=1)
# Output:
# obj.sil: silhouette of each object 
# fuzzy.sil: fuzzy silhouette of the given fuzzy partition
# The index takes values between -1 and 1. The goal is to maximize it.
##################################################################

Fuzzy.Sil <- function(dist, U, param = 1) {
  U <- as.matrix(U)  # Membership matrix
  n <- nrow(U)       # Number of objects
  k <- ncol(U)       # Number of clusters
  
  D <- as.matrix(dist)  # Dissimilarity matrix
  
  # Determine the cluster where each object has the highest membership degree
  clusters.obj <- apply(U, 1, which.max)
  
  # Initialize vectors to store the results
  dif.degrees <- rep(0, n)
  a <- rep(0, n)  # Intra-cluster distance
  b <- rep(0, n)  # Inter-cluster distance
  obj.sil <- rep(0, n)  # Silhouette score for each object
  
  # Compute silhouette values for each object
  for (i in 1:n) {
    
    # Calculate the degree difference for object i
    dif.degrees[i] <- (max(U[i, ]) - max(U[i, ][-which.max(U[i, ])]))^param
    
    # Calculate inter-cluster distances
    B <- rep(0, k)
    
    # Identify indices of clusters different from the one of object i
    i2 <- which(clusters.obj != clusters.obj[i])  
    c2 <- unique(clusters.obj[i2])  # Unique clusters that object i does not belong to
    
    # For each of these clusters, compute the average distance
    for (c in c2) {
      i3 <- which(clusters.obj == c)  # Indices of objects in the c-th cluster
      B[c] <- mean(D[i, i3])  # Average distance to objects in cluster c
    }   
    
    # Minimum inter-cluster distance for object i (excluding its own cluster)
    b[i] <- min(B[-clusters.obj[i]])
    
    # Intra-cluster distance for object i (mean distance to other objects in the same cluster)
    m <- which(clusters.obj == clusters.obj[i])  # Indices of objects in the same cluster as i
    a[i] <- mean(D[i, m][-which(m == i)])  # Exclude itself from the average
    
    # Compute the silhouette score for object i if it has other members in its cluster
    if (length(m) > 1) {
      obj.sil[i] <- (b[i] - a[i]) / max(a[i], b[i])
    }
  }
  
  # Compute the final fuzzy silhouette score
  fuzzy.sil <- sum(dif.degrees * obj.sil) / sum(dif.degrees)
  
  # Prepare output
  out <- list()
  out$obj.sil <- obj.sil  # Individual silhouette scores
  out$fuzzy.sil <- fuzzy.sil  # Overall fuzzy silhouette score
  
  return(out)  # Return both the individual and the overall score
}

##################################################################
### FUZZY K-MEDOIDS with optimal K selected by Fuzzy Silhouette  
##################################################################
## Inputs
# Y: Matrix or data.frame (rows=objects to cluster)
# m: Parameter of fuzziness (usually, 1.5)
# k: number of clusters
# dist.matrix: dissimilarity matrix
# param: weighting parameter of the fuzzy silhoutte (default: 1)
# RS:	number of (random) starts (default: 100)
# max.iter: Maximum number of iterations (default: 100)
## Output
#FS.Values: values of the fuzzy silhoutte index for k=2,...,kmax
#Best.K: optimal number of cluster
#Best.C: optimal cluster output (as in Fcmdd)
##################################################################
Fcmdd.Sil <- function(Y, m, kmax, dist.matrix, param=1, RS = 100, max.iter = 100){
  n <- nrow(Y)
  s <- ncol(Y)
  
  # Check row names and set if NULL
  if (is.null(row.names(Y))) row.names(Y) <- 1:n
  
  #X <- as.matrix(Y)
  dist_matrix <- as.matrix(dist.matrix)
  rownames(dist_matrix) <- rownames(Y)
  colnames(dist_matrix) <- rownames(Y)
  
  # Preallocate vector for storing fuzzy silhouette values
  FS.vec <- numeric(kmax - 1)
  
  # Inizialization for best fuzzy silhoutte index
  Best.FS<- -Inf
  Best.Clust<-list()
  
  # Loop from 2 to kmax to compute fuzzy silhouette for each k
  for (k_val in 2:kmax) {
    # Call Fcmdd once for each k value, using precomputed distance matrix
    Aux <- Fcmdd(Y, m = m, k = k_val, dist.matrix = dist_matrix)
    
    # Compute fuzzy silhouette once and store it
    #Sil <- Fuzzy.Sil(Aux, param)$FS
    Sil <- Fuzzy.Sil(Aux$distance,Aux$U)$fuzzy.sil
    
    
    FS.vec[k_val - 1] <- Sil
    
    # Check if this is the best silhouette score
    if (Sil > Best.FS) {
      Best.FS <- Sil
      Best.Clust <- Aux
    }
  }
  
  # Combine k values with their corresponding fuzzy silhouette values
  Mat.aux <- cbind(2:kmax, FS.vec)
  colnames(Mat.aux) <- c("k", "fuzzy silhouette")
  
  return(list(FS.Values=Mat.aux,
              Best.K=which.max(FS.vec)+1,
              Best.C=Best.Clust))
}
