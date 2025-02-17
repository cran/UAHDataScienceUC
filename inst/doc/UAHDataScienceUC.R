## ----eval = FALSE-------------------------------------------------------------
#  install.packages("UAHDataScienceUC")

## -----------------------------------------------------------------------------
# Load library
library(UAHDataScienceUC)

# Load data
data(db5)

# Create sample data
data <- db5[1:10, ]

## -----------------------------------------------------------------------------
# Perform k-means clustering
result <- kmeans_(data, centers = 3, max_iterations = 10)

# Plot results
plot(data, col = result$cluster, pch = 20)
points(result$centers, col = 1:3, pch = 8, cex = 2)

## -----------------------------------------------------------------------------
# Perform hierarchical clustering
result <- agglomerative_clustering(
  data,
  proximity = "single",
  distance_method = "euclidean",
  learn = TRUE
)

## -----------------------------------------------------------------------------
result <- dbscan(
  data,
  epsilon = 0.3,
  min_pts = 4,
  learn = TRUE
)

## -----------------------------------------------------------------------------
result <- gaussian_mixture(
  data,
  k = 3,
  max_iter = 100,
  learn = TRUE
)

# Plot results with contours
plot(data, col = result$cluster, pch = 20)

## -----------------------------------------------------------------------------
result <- genetic_kmeans(
  data,
  k = 3,
  population_size = 10,
  mut_probability = 0.5,
  max_generations = 10,
  learn = TRUE
)

## -----------------------------------------------------------------------------
# Create sample data
data <- matrix(c(1,2,1,4,5,1,8,2,9,6,3,5,8,5,4), ncol=3)
dataFrame <- data.frame(data)
target <- c(1,2,3)
weights <- c(0.1, 0.6, 0.3)

# Perform correlation clustering
result <- correlation_clustering(
    dataFrame,
    target = target,
    weight = weights,
    distance_method = "euclidean",
    normalize = TRUE,
    learn = TRUE
)

## -----------------------------------------------------------------------------
# Using different distance metrics
agglomerative_clustering(data, distance_method = "euclidean")
agglomerative_clustering(data, distance_method = "manhattan")
agglomerative_clustering(data, distance_method = "canberra")
agglomerative_clustering(data, distance_method = "chebyshev")

