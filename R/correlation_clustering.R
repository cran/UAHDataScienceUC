#' @title Hierarchical Correlation Clustering
#' @description Performs hierarchical correlation clustering by applying weights, distance metrics, and other parameters
#' to analyze relationships between data points and a target.
#'
#' @param data A data frame containing the main data
#' @param target A data frame, numeric vector or matrix to use as correlation target. Default is NULL.
#' @param weight A numeric vector of weights. Default is empty vector.
#' @param distance_method A string specifying the distance metric to use. Options are:
#'   \itemize{
#'     \item "euclidean" - Euclidean distance
#'     \item "manhattan" - Manhattan distance
#'     \item "canberra" - Canberra distance
#'     \item "chebyshev" - Chebyshev distance
#'   }
#' @param normalize A boolean parameter indicating whether to normalize weights. Default is TRUE.
#' @param labels A string vector for graphical solution labeling. Default is NULL.
#' @param learn A boolean indicating whether to show detailed algorithm explanations. Default is FALSE.
#' @param waiting A boolean controlling pauses between explanations. Default is TRUE.
#'
#' @details This function executes the complete hierarchical correlation method in the following steps:
#' \enumerate{
#'   \item The function transforms data into useful objects
#'   \item Creates the clusters
#'   \item Calculates the distance from the target to every cluster using the specified distance metric
#'   \item Orders the distances in ascending order
#'   \item Orders the clusters according to their distance from the previous step
#'   \item Shows the sorted clusters and the distances used
#' }
#'
#' @return An R object containing:
#' \itemize{
#'   \item dendrogram - A hierarchical clustering dendrogram
#'   \item sortedValues - A data frame with the sorted cluster values
#'   \item distances - A data frame with the sorted distances
#' }
#'
#' @examples
#' data <- matrix(c(1,2,1,4,5,1,8,2,9,6,3,5,8,5,4), ncol=3)
#' dataFrame <- data.frame(data)
#' target1 <- c(1,2,3)
#' target2 <- dataFrame[1,]
#' weight1 <- c(1,6,3)
#' weight2 <- c(0.1,0.6,0.3)
#'
#' # Basic usage
#' correlation_clustering(dataFrame, target1)
#'
#' # With weights
#' correlation_clustering(dataFrame, target1, weight1)
#'
#' # Without weight normalization
#' correlation_clustering(dataFrame, target1, weight1, normalize = FALSE)
#'
#' # Using Canberra distance with weights
#' correlation_clustering(dataFrame, target1, weight2, distance = "canberra", normalize = FALSE)
#'
#' # With detailed explanations
#' correlation_clustering(dataFrame, target1, learn = TRUE)
#'
#' @author Original authors:
#' \itemize{
#'   \item Roberto Alcantara \email{roberto.alcantara@@edu.uah.es}
#'   \item Juan Jose Cuadrado \email{jjcg@@uah.es}
#'   \item Universidad de Alcala de Henares
#' }
#'
#' @export
correlation_clustering <- function(data, target = NULL, weight = c(),
                                   distance_method = "euclidean", normalize = TRUE,
                                   labels = NULL, learn = FALSE, waiting = FALSE) {

  res = list()

  # Normalize weights and provide explanation if learn=TRUE
  if(learn) {
    console.log("EXPLANATION:")
    console.log("")
    console.log("The Correlation Hierarchical Clustering algorithm is a classification technique that:")
    console.log("")
    console.log("    1. Initializes a cluster for each data point")
    console.log("    2. Calculates distances between clusters and a given target")
    console.log("    3. Applies weights to achieve weighted results")
    console.log("")

    if(normalize) {
      console.log("Due to normalize = TRUE, weights will be normalized to [0,1] range")
    } else {
      console.log("Due to normalize = FALSE, weights will remain unchanged")
    }
    console.log("")

    if(waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  weight <- if(learn) {
    normalizeWeight.details(normalize, weight, data)
  } else {
    normalizeWeight(normalize, weight, data)
  }

  if(learn) {
    hline()
    console.log("WEIGHTS:")
    console.log("")
    console.log("The following weights will be used:")
    print(weight)
    console.log("")

    # Add distance formula display
    display_distance_formula(distance_method)
    console.log("This distance metric will be used to:")
    console.log("    - Calculate distances between each cluster and the target")
    console.log("    - Sort clusters by their similarity to the target")
    console.log("")

    if(waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Input validation with explanations
  if(is.null(target) | is.null(data)) {
    if(learn) {
      hline()
      console.log("ERROR:")
      console.log("")
      console.log("Invalid input: target or data is NULL")
      console.log("")
    }
    return(NULL)
  }

  # Initialize data and target
  list <- if(learn) {
    initData.details(data)
  } else {
    initData(data)
  }
  target <- if(learn) {
    initTarget.details(target, data)
  } else {
    initTarget(target, data)
  }

  if(learn) {
    hline()
    console.log("INITIALIZATION:")
    console.log("")
    console.log("Initialized data:")
    print(list)
    console.log("")
    console.log("Initialized target:")
    print(target)
    console.log("")

    if(waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Calculate distances between clusters and target
  distances <- c()

  if(learn) {
    hline()
    console.log("DISTANCE CALCULATION:")
    console.log("")
    console.log(sprintf("Calculating distances between clusters and target using %s distance", distance_method))
    console.log("")

    if(waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Import proxy package
  if(!requireNamespace("proxy", quietly = TRUE)) {
    stop("Package 'proxy' is required for distance calculations")
  }

  distance_map <- c(
    "euclidean" = "euclidean",
    "manhattan" = "manhattan",
    "canberra" = "canberra",
    "chebyshev" = "maximum"
  )


  # Validate distance method
  valid_distances <- c("euclidean", "manhattan", "canberra", "chebyshev")
  distance <- match.arg(distance_method, valid_distances)
  proxy_distance <- distance_map[distance]

  # Convert target to matrix format for proxy::dist
  target_matrix <- as.matrix(target)

  for (index in seq_along(list)) {
    # Convert current cluster to matrix format
    cluster_matrix <- as.matrix(list[[index]])

    # Calculate distance using proxy package
    if(!is.null(weight)) {
      # Apply weights to the calculation
      dist_result <- proxy::dist(cluster_matrix, target_matrix,
                                 method = proxy_distance,
                                 weights = weight)
    } else {
      dist_result <- proxy::dist(cluster_matrix, target_matrix,
                                 method = proxy_distance)
    }

    # Extract the distance value
    dist <- as.numeric(dist_result)
    distances <- c(distances, dist)
  }

  if(learn) {
    hline()
    console.log("DISTANCES:")
    console.log("")
    console.log("Calculated distances:")
    print(distances)
    console.log("")
    console.log("Sorted distances:")
    print(sort(distances))
    console.log("")
  }

  # Sort distances and prepare clustering results
  sortedDistances <- sort(distances)
  values <- data.frame()
  names <- names(data)
  clusters <- list()
  dendrogramLabels <- c()
  clustersDendrogram <- data.frame()

  if(learn) {
    message("\n Then, using sorted distances, the function order the clusters. \n \n")
    if(waiting) {
      readline(prompt="Press [Enter] to continue...")
    }
  }

  # Initialize data structures
  sortedDistances <- sort(distances)
  values <- data.frame()
  names <- names(data)
  clusters <- list()
  dendrogramLabels <- c()
  clustersDendrogram <- data.frame()

  # Process each cluster based on sorted distances
  for (cluster in seq_along(list)) {
    distance <- sortedDistances[cluster]

    # Find cluster position based on distance
    clust <- getClusterPosition(distance, distances, list)
    list[[clust[[2]]]] <- data.frame(list[[clust[[2]]]])

    # Build dendrogram structure
    if(cluster == 1) {
      data <- c()
      data <- c(-clust[[2]])
    } else if(cluster == 2) {
      data <- c(data, -clust[[2]])
      clustersDendrogram <- rbind(clustersDendrogram,
                                  data.frame(matrix(data, ncol=2)))
    } else {
      clustersDendrogram <- rbind(clustersDendrogram,
                                  data.frame(matrix(c(-clust[[2]], cluster-2), ncol=2)))
    }

    # Update tracking structures
    dendrogramLabels <- c(dendrogramLabels, clust[[2]])
    values <- rbind(values, data.frame(clust[[1]]))
    clusters[[length(clusters) + 1]] <- data.frame(clust[[1]])
  }

  # Set column names for values
  names(values) <- names

  # Create distance results
  distances_df <- data.frame(cluster=dendrogramLabels, sortedDistances)

  # Create sorted values results
  sorted_values <- data.frame(cluster=dendrogramLabels, values)
  if(is.data.frame(data)) {
    names(sorted_values) <- names(data)
  }

  # Create dendrogram structure
  dendrogram <- list()
  dendrogram$merge <- as.matrix(clustersDendrogram)
  dendrogram$height <- seq_len(length(list)-1)
  dendrogram$order <- seq_len(length(list))
  dendrogram$labels <- if(!is.null(labels)) labels else seq_len(length(list))
  class(dendrogram) <- "hclust"

  # Prepare final result object
  result <- list(
    dendrogram = dendrogram,
    sortedValues = sorted_values,
    distances = distances_df
  )

  if(learn) {
    hline()
    console.log("RESULTS:")
    console.log("")
    console.log("Final sorted distances:")
    print(result$distances)
    console.log("")
    console.log("Final sorted clusters:")
    print(result$sortedValues)
    console.log("")
    console.log("Dendrogram visualization:")
    plot(result$dendrogram)
    console.log("")

    if(waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }

    hline()
  } else {
    plot(result$dendrogram)
  }

  return(result)
}

getClusterPosition <- function(distance,vector,list){
  found <- FALSE
  index <- 1
  while(!found & (index < length(vector))){

    if(!is.data.frame(list[[index]]) & vector[index] == distance){
      found <- TRUE
    } else {
      index <- index + 1
    }
  }
  cluster <- list[[index]]
  res <- list(cluster , index)
  res
}


initData <- function(data){
  solution <- list()
  for (row in c(1:nrow(data))) {
    rowData <- data[row,]
    values <- c()
    for (column in c(1:ncol(data))) {
      values <- c(values, rowData[,column])
    }
    listData <- matrix(values, ncol = ncol(data))
    solution[[length(solution) + 1]] <- listData
  }

  solution
}


initData.details <- function(data){
  console.log("DATA INITIALIZATION:")
  console.log("")
  console.log("Input data:")
  print(data)
  console.log("")

  solution <- list()
  for (row in c(1:nrow(data))) {
    rowData <- data[row,]
    values <- c()
    for (column in c(1:ncol(data))) {
      values <- c(values, rowData[,column])
    }
    listData <- matrix(values, ncol = ncol(data))
    solution[[length(solution) + 1]] <- listData
  }
  console.log("Each cluster is initialized as a matrix with one row and the same columns as the input data")
  console.log("")
  console.log("Initialized clusters:")
  print(solution)
  console.log("")
  solution
}


initTarget <- function(target,data){
  if(is.data.frame(target)){
    target <- initData(target)[[1]]
  } else if (is.matrix(target)){
    target <- target
  } else if (is.vector(target)){
    target <- matrix(target, ncol = length(target))
  }
  if(nrow(target) == 1 & ncol(target) == ncol(data)){
    target <- target
  } else {
    console.log('\nInvalid target dimensions!\n')
    target <- matrix(rep(0,ncol(data)),ncol = ncol(data))
  }
  target
}


initTarget.details <- function(target,data){
  hline()
  console.log("TARGET INITIALIZATION:")
  console.log("")
  console.log("Input target:")
  print(target)
  console.log("")

  print(target)
  if(is.data.frame(target)){
    console.log("Converting data frame target to matrix format")
    target <- initData(target)[[1]]
  } else if (is.matrix(target)){
    console.log("Target is already in matrix format")
    target <- target
  } else if (is.vector(target)){
    console.log("Converting vector target to matrix format")
    target <- matrix(target, ncol = length(target))
  }

  console.log("")
  console.log("Validating target dimensions:")

  if(nrow(target) == 1 & ncol(target) == ncol(data)){
    console.log("Target has one row")
    console.log("Target has same number of columns as input data")
    target <- target
  } else {
    console.log("! Invalid target dimensions")
    console.log("! Initializing default target with zeros")
    target <- matrix(rep(0,ncol(data)),ncol = ncol(data))
  }
  console.log("")
  console.log("Final target:")
  print(target)
  console.log("")

  target
}

normalizeWeight <- function(normalize,weight,data){
  if(is.null(weight)){
    weight <- c()
    value <- 1
    for (j in c(1:ncol(data))) {
      weight <- c(weight, value)
    }
  }
  res <- c()
  if(normalize){
    total <- sum(weight)
    for (index in c(1:length(weight))) {
      currentValue <- weight[index]
      normalizedValue <- currentValue/total
      res  <- c(res, normalizedValue)
    }

  } else {
    res <- weight
  }
  res
}


normalizeWeight.details <- function(normalize,weight,data){
  hline()
  console.log("WEIGHT NORMALIZATION:")
  console.log("")

  if(!is.null(weight)){
    console.log("Initial weights:")
    for (i in c(1:length(weight))) {
      console.log(sprintf("  Weight %d: %f", i, weight[i]))
    }
  }
  console.log("No weights provided")
  console.log(sprintf("Initializing %d weights with value 1", ncol(data)))
  if(is.null(weight)){
    weight <- c()
    value <- 1
    for (j in c(1:ncol(data))) {
      weight <- c(weight, value)
    }
  }

  console.log("")

  res <- c()
  if(normalize){
    total <- sum(weight)
    console.log("Normalizing weights:")
    console.log(sprintf("Total sum: %f", total))
    console.log("Formula: weight[i] = weight[i] / total")
    console.log("")
    for (index in c(1:length(weight))) {
      currentValue <- weight[index]
      normalizedValue <- currentValue/total
      res  <- c(res, normalizedValue)
    }

  } else {
    console.log("Weights will not be normalized")
    console.log("")
    console.log("Final weights:")
    res <- weight
  }
  message("\n These are the new weights: \n")
  for (i in c(1:length(res))) {
    console.log(sprintf("  Weight %d: %f", i, weight[i]))
  }

  console.log("")
  res
}
