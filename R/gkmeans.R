#' @title Genetic K-Means Clustering
#'
#' @description Performs Genetic K-Means clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param k the number of clusters.
#' @param population_size the number of individuals in the population.
#' @param mut_probability the probability of a mutation occurring.
#' @param max_generations the maximum number of iterations allowed.
#' @param learn a Boolean determining whether intermediate logs explaining how
#' the algorithm works should be printed or not.
#' @param waiting a Boolean determining whether the intermediate logs should be
#' printed in chunks waiting for user input before printing the next or not.
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @return A kmeans object as returned by the original kmeans function.
#'
#' @examples
#' ### Example 1: Simple usage with circles dataset
#' result1 <- genetic_kmeans(db1[1:20,], 2, learn = TRUE, waiting = FALSE)
#'
#' ### Example 2: Moons dataset with different population size
#' result2 <- genetic_kmeans(db2[1:20,], 2, population_size = 20,
#'                          learn = TRUE, waiting = FALSE)
#'
#' ### Example 3: Varying density clusters with different mutation probability
#' result3 <- genetic_kmeans(db3[1:20,], 3, mut_probability = 0.7,
#'                          learn = TRUE, waiting = FALSE)
#'
#' ### Example 4: Well-separated clusters with larger population
#' result5 <- genetic_kmeans(db5[1:20,], 3, population_size = 30,
#'                          mut_probability = 0.6, learn = TRUE, waiting = FALSE)
#'
#' ### Example 5: Using different parameters combinations
#' result6 <- genetic_kmeans(db1[1:20,], 2,
#'                          population_size = 15,
#'                          mut_probability = 0.8,
#'                          max_generations = 15,
#'                          learn = TRUE,
#'                          waiting = FALSE)
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
#' @importFrom stats runif sd prcomp predict
#' @importFrom graphics points legend
#' @export
genetic_kmeans <- function(
    data,
    k,
    population_size = 10,
    mut_probability = .5,
    max_generations = 10,
    learn = FALSE,
    waiting = TRUE,
    ...
) {
  # Helper function to round and format numbers
  format_number <- function(x) format(round(x, 2), nsmall = 2)

  if (learn) {
    hline()
    console.log("EXPLANATION:")
    console.log("")
    console.log("The Genetic K-Means algorithm combines the K-Means clustering method with genetic algorithm concepts.")
    console.log("It follows these main steps:")
    console.log("1. Initialize a population of random cluster assignments")
    console.log("2. Evaluate the fitness of each individual based on within-cluster variation")
    console.log("3. Select parents for the next generation based on fitness")
    console.log("4. Apply mutation and crossover to create new individuals")
    console.log("5. Repeat steps 2-4 for a specified number of generations")
    console.log("")
    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
    }
  }

  # Initialize the population
  population <- gka_initialization(nrow(data), population_size, k)

  if (learn) {
    hline()
    console.log("INITIALIZATION:")
    console.log("")
    console.log(paste("A population of", population_size, "individuals has been randomly initialized."))
    console.log("Each individual represents a possible clustering solution.")
    console.log("Here's a sample of the initial population (first 5 individuals, first 10 data points):")
    print(population[1:min(5, population_size), 1:min(10, ncol(population))])
    console.log("")
    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
    }
  }

  # Compute the fitness of the initial population
  centers <- gka_centers(data, k, population)
  twcv <- gka_twcv(data, k, population, centers)
  fitness <- gka_fitness(twcv)

  # Choose the best individual as the initial solution
  idx <- which.max(fitness)
  solution <- structure(
    list(
      cluster = population[idx, ],
      centers = centers[idx, , ],
      totss = NULL,
      withinss = NULL,
      tot.withinss = twcv[idx],
      betweenss = NULL,
      size = NULL,
      iter = NULL,
      ifault = 0
    ),
    class = "kmeans"
  )

  if (learn) {
    console.log("INITIAL FITNESS:")
    console.log(paste("Best fitness:", format_number(max(fitness))))
    console.log(paste("Average fitness:", format_number(mean(fitness))))
    console.log("")
    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
    }
  }

  # Run the algorithm for a maximum of max_generations generations
  for (generation in 1:max_generations) {
    if (learn) {
      hline()
      console.log(paste("GENERATION:", generation))
      console.log("")
    }

    # Selection
    idx <- gka_selection(population_size, fitness)
    chromosome <- population[idx, ]
    centers <- centers[idx, , ]

    if (learn) {
      console.log("SELECTION:")
      console.log("Parents for the next generation are selected based on their fitness.")
      console.log("Selected parent (first 10 data points):")
      print(chromosome[1:min(10, length(chromosome))])
      console.log("")
      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
      }
    }

    # Mutation
    prob <- gka_allele_mutation(data, k, centers, ...)
    for (i in seq_len(population_size)) {
      population[i, ] <- gka_mutation(chromosome, prob, k, mut_probability)
    }
    population <- gka_chromosome_fix(population, k)

    if (learn) {
      console.log("MUTATION:")
      console.log("Random mutations are applied to the chromosomes with probability", format_number(mut_probability))
      console.log("Sample of mutated population (first 5 individuals, first 10 data points):")
      print(population[1:min(5, population_size), 1:min(10, ncol(population))])
      console.log("")
      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
      }
    }

    # Crossover
    centers <- gka_centers(data, k, population)
    for (i in seq_len(population_size)) {
      population[i, ] <- gka_crossover(data, centers[i, , ])
    }
    population <- gka_chromosome_fix(population, k)

    if (learn) {
      console.log("CROSSOVER:")
      console.log("K-Means Operator (KMO) is applied as a form of crossover.")
      console.log("This reassigns each point to its nearest center.")
      console.log("Sample of population after crossover (first 5 individuals, first 10 data points):")
      print(population[1:min(5, population_size), 1:min(10, ncol(population))])
      console.log("")
      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
      }
    }

    # Compute the fitness of the new population
    centers <- gka_centers(data, k, population)
    twcv <- gka_twcv(data, k, population, centers)
    fitness <- gka_fitness(twcv)
    idx <- which.max(fitness)

    # Update the solution
    if (twcv[idx] < solution$tot.withinss) {
      idx <- which.max(fitness)
      solution$cluster <- population[idx, ]
      solution$centers <- centers[idx, , ]
      solution$tot.withinss <- twcv[idx]
    }

    if (learn) {
      console.log("FITNESS EVALUATION:")
      console.log("The fitness of each individual is calculated based on within-cluster variation.")
      console.log(paste("Best fitness in this generation:", format_number(max(fitness))))
      console.log(paste("Average fitness:", format_number(mean(fitness))))
      console.log(paste("Total Within-Cluster Variation of best solution:", format_number(solution$tot.withinss)))
      console.log("")
      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
      }
    }
  }

  # Finalize the solution
  row.names(solution$centers) <- seq_len(k)
  center <- apply(data, 2, mean)
  solution$totss <- sum(apply(data, 1, function(x) x - center) ^ 2)
  solution$withinss <- sapply(
    seq_len(k),
    function(cluster) {
      ccenter <- solution$centers[cluster, ]
      cdata <- data[solution$cluster == cluster, , drop = FALSE]
      sum(apply(cdata, 1, function(x) x - ccenter)^2)
    }
  )
  solution$betweenss <- solution$totss - solution$tot.withinss
  solution$size <- as.integer(table(factor(solution$cluster, levels = seq_len(k))))
  solution$iter <- max_generations

  if (learn) {
    hline()
    console.log("FINAL RESULTS:")
    console.log("")
    console.log(paste("Number of clusters:", k))
    console.log(paste("Total sum of squares:", format_number(solution$totss)))
    console.log(paste("Total within-cluster sum of squares:", format_number(solution$tot.withinss)))
    console.log(paste("Between-cluster sum of squares:", format_number(solution$betweenss)))
    console.log("Cluster sizes:")
    print(solution$size)
    console.log("Final cluster centers:")
    print(apply(solution$centers, c(1,2), format_number))
    console.log("")

    # Plotting the results
    console.log("Generating plot of clustering results...")

    # If data has more than 2 dimensions, use PCA for visualization
    if (ncol(data) > 2) {
      pca_result <- prcomp(data, scale. = TRUE)
      plot_data <- pca_result$x[, 1:2]
      x_label <- "PC1"
      y_label <- "PC2"
    } else {
      plot_data <- data
      x_label <- "X"
      y_label <- "Y"
    }

    # Create the plot
    plot(plot_data, col = solution$cluster, pch = 20,
         xlab = x_label, ylab = y_label,
         main = "Genetic K-Means Clustering Results")

    # Add cluster centers to the plot
    if (ncol(data) > 2) {
      centers_pca <- predict(pca_result, solution$centers)[, 1:2]
    } else {
      centers_pca <- solution$centers
    }
    points(centers_pca, col = 1:k, pch = 3, cex = 2, lwd = 2)

    # Add a legend
    legend("topright", legend = paste("Cluster", 1:k),
           col = 1:k, pch = 20, bty = "n")

    console.log("Plot generated. Check the graphics device.")
    console.log("")
  }

  solution
}

#' @title Initialization method
#'
#' @param n the number of observations in the data.
#' @param p the number of individuals in the population.
#' @param k the number of clusters.
#'
#' @return a matrix of size \code{p} by \code{n} with the cluster assignments
#' for each observation.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_initialization <- function(n, p, k) {
  population <- matrix(0, nrow = p, ncol = n)

  # For each member of the population
  for (i in seq_len(p)) {
    # Assign p observations to each cluster, where p is the floor of n / k
    quotient <- n %/% k
    pigeons <- sample(seq_len(n), quotient * k)
    for (j in seq_len(k)) {
      from <- (j - 1) * quotient + 1
      to   <- from    + quotient - 1
      population[i, pigeons[seq(from, to)]] <- j
    }

    # Assign the remaining observations to a random cluster
    remainder <- n %% k
    population[i, -pigeons] <- sample(seq_len(k), remainder)
  }

  population
}

#' @title Centroid computation
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation. The matrix is of size \code{n} by \code{m}.
#' @param k the number of clusters.
#' @param population a matrix of size \code{p} by \code{n} with the cluster
#' assignments for each observation.
#'
#' @return a 3D array of size \code{p} by \code{k} by \code{m}.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_centers <- function(data, k, population) {
  centers <- array(0, dim = c(nrow(population), k, ncol(data)))
  for (i in seq_len(nrow(population))) {
    for (j in seq_len(k)) {
      centers[i, j, ] <- colMeans(data[population[i, ] == j, ])
    }
  }
  centers
}

#' @title Total Within Cluster Variation (TWCV) computation
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation. The matrix is of size \code{n} by \code{m}.
#' @param k the number of clusters.
#' @param population a matrix of size \code{p} by \code{n} with the cluster
#' assignments for each observation.
#' @param centers a 3D array of size \code{p} by \code{k} by \code{m} with the
#' cluster centers for each individual in the population.
#'
#' @return a vector of size \code{p} with the total within cluster variation of
#' each individual in the population.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_twcv <- function(data, k, population, centers) {
  twcv <- numeric(nrow(population))
  for (i in seq_len(nrow(population))) {
    twcv[i] <- sum(
      sapply(
        seq_len(k),
        function(j) sum((t(data[population[i, ] == j, ]) - centers[i, j, ]) ^ 2)
      )
    )
  }
  twcv
}

#' @title Fitness function
#'
#' @param twcv a vector of size \code{p} with the total within cluster
#' variation of each individual in the population.
#'
#' @return a vector of size \code{p} with the fitness of each individual in the
#' population.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_fitness <- function(twcv) {
  f <- -twcv
  g <- f - (mean(f) - 2 * sd(f))
  g[g <= 0] <- 1e-6
  g
}

#' @title Selection method
#'
#' @param p the number of individuals in the population.
#' @param fitness a vector of size \code{p} with the fitness of each individual
#' in the population.
#'
#' @return the index of the individual selected for reproduction.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_selection <- function(p, fitness) {
  sample(seq_len(p), 1, prob = fitness / sum(fitness))
}

#' @title Allele mutation probability computation
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation. The matrix is of size \code{n} by \code{m}.
#' @param k the number of clusters.
#' @param centers a matrix of size \code{k} by \code{m} with the cluster centers
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @return a matrix of size \code{n} by \code{k} with the probability of each
#' allele mutating to a specific cluster.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
gka_allele_mutation <- function(data, k, centers, ...) {
  d <- proxy::dist(data, centers, ...)
  dmax <- apply(d, 1, max)
  tmp1 <- matrix(0, nrow = nrow(data), ncol = k)
  for (j in seq_len(k)) {
    tmp1[, j] <- 1.5 * dmax - d[, j]
  }
  tmp1[tmp1 <= 0] <- 1e-6
  tmp2 <- rowSums(tmp1)
  prob <- tmp1 / tmp2
  prob
}

#' @title Mutation method
#'
#' @param chromosome a vector of size \code{n} with the cluster assignments for
#' each observation.
#' @param prob a matrix of size \code{n} by \code{k} with the probability of
#' each allele mutating to a specific cluster.
#' @param k the number of clusters.
#' @param mut_probability the probability of a mutation occurring.
#'
#' @return a vector of size \code{n} with the cluster assignments for each
#' observation i.e. a new chromosome.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom stats runif
gka_mutation <- function(chromosome, prob, k, mut_probability) {
  new_chromosome <- chromosome
  for (i in seq_len(length(new_chromosome))) {
    if (runif(1) < mut_probability) {
      new_chromosome[i] <- sample(seq_len(k), 1, prob = prob[i, ])
    }
  }
  new_chromosome
}

#' @title Crossover method i.e. K-Means Operator
#'
#' @description K-Means Operator (KMO) which replaces the crossover operator in
#' the Genetic K-Means algorithm (GKA).
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation. The matrix is of size \code{n} by \code{m}.
#' @param centers a matrix of size \code{k} by \code{m} with the cluster centers
#' for a specific individual in the population.
#'
#' @return a vector of size \code{n} with the cluster assignments for each
#' observation i.e. a new chromosome.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
gka_crossover <- function(data, centers) {
  d <- proxy::dist(data, centers)
  apply(d, 1, which.min)
}

#' @title Chromosome fixing method
#'
#' @description This method fixes chromosomes which do not have at least one
#' observation assigned to each cluster.
#'
#' @param population a matrix of size \code{p} by \code{n} with the cluster
#' assignments for each observation.
#' @param k the number of clusters.
#'
#' @return a matrix of size \code{p} by \code{n} with the cluster assignments
#' for each observation.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_chromosome_fix <- function(population, k) {
  for (i in seq_len(nrow(population))) {
    sdiff <- setdiff(seq_len(k), unique(population[i, ]))
    if (length(sdiff) <= 0)
      next

    for (cluster in sdiff) {
      n <- ncol(population)
      population[i, sample(seq_len(n), n %/% k)] <- cluster
    }
  }
  population
}
