#' Generates simulated data with clusters of co-varying data.
#'
#' @param nsample The number of loci to simulate (e.g. 20 genes).
#' @param nfeature The number of features to generate (e.g. the
#'  number of RNA-seq samples).
#' @param clusters A matrix/data frame containing the cluster
#'  boundaries and, optionally, the cluster co-variance.
#'  Must have the columns c('start', 'end', 'sigma')
#'  sigma is only necessary if specifying the cluster
#'  covariance. If this column is missing or the value
#'  is NA, sigma_intraclust is used.
#' @param sigma The diagonal variance.
#' @param sigma_intraclust The global within cluster co-variance.
#'  Can be over-ridden using the clusters 'sigma' column.
#' @param sigma_interclust The covariance between all clusters.
#' @param seed Seed the random number generator before running.
#'
#' @returns List with covariance matrix, nsample by nfeature matrix
#'  of features with covariance following covariance matrix, and
#'  the true cluster boundaries (unchanged).
#'
#' For more complex clusters and covariance structures, it will be easier
#' to generate data yourself with a designed covariance matrix and the
#' MASS::mvrnorm function.
#'
#' @export
generate_simple_clusters <- function(
    nsample,
    nfeature,
    clusters,
    sigma = 1,
    sigma_intraclust = 0.99,
    sigma_interclust = 0.0,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  covariance <- matrix(sigma_interclust, nrow = nsample, ncol = nsample)

  has_indiv_sigma <- 'sigma' %in% colnames(clusters)

  for (i in seq_len(nrow(clusters))) {

    start <- clusters[i, 'start']
    end <- clusters[i, 'end']

    if (has_indiv_sigma) {
      cluster_sigma <- clusters[i, 'sigma']
      if (is.na(cluster_sigma)) {
        cluster_sigma <- sigma_intraclust
      }
    } else {
      cluster_sigma <- sigma_intraclust
    }

    covariance[start:end, start:end] <- cluster_sigma
  }

  diag(covariance) <- sigma

  mu <- rnorm(nsample, mean = 0, sd = 1)
  features <- MASS::mvrnorm(n = nfeature, mu = mu, Sigma = covariance)

  return(list(
    clusters = clusters,
    covariance = covariance,
    features = features
  ))
}


#' Here's the case of no-clusters (all points in separate clusters)
#' @param seed Initialise the random number generator
no_clusters <- function(seed = NULL) {
  clusters <- matrix(
    sort(rep(1:20, 2)),
    nrow = 20,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:20, c("start", "end"))
  )

  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}


#' here's the case of one single cluster spanning the full range of the window.
#' @param seed Initialise the random number generator
one_cluster <- function(seed = NULL) {
  clusters <- matrix(
    c(1, 20),
    nrow = 1,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1, c("start", "end"))
  )

  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}


#' here's the case of two completely separate clusters, directly next to each
#' other
#' @param seed Initialise the random number generator
two_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      1, 10,
      11, 20
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:2, c("start", "end"))
  )


  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}


#' here's the case of five completely separate clusters, directly next to each
#' other
#' @param seed Initialise the random number generator
five_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      1, 10,
      11, 20,
      21, 30,
      31, 40,
      41, 50
    ),
    nrow = 5,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:5, c("start", "end"))
  )


  return(generate_simple_clusters(
      50,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}

#' Case of two clusters directly next to each other,
#' with some correlation between them.
#' @param seed Initialise the random number generator
two_covarying_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      1, 10,
      11, 20
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:2, c("start", "end"))
  )

  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.6,
      seed = seed
  ))
}


#' Case of five clusters directly next to each other,
#' with some correlation between them.
#' @param seed Initialise the random number generator
five_covarying_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      1, 10,
      11, 20,
      21, 30,
      31, 40,
      41, 50
    ),
    nrow = 5,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:5, c("start", "end"))
  )

  return(generate_simple_clusters(
      50,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.6,
      seed = seed
  ))
}


#' here's the case of one single cluster in the center of the window, with
#' some space on either side.
#' @param seed Initialise the random number generator
one_gappy_cluster <- function(seed = NULL) {
  clusters <- matrix(
    c(6, 15),
    nrow = 1,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1, c("start", "end"))
  )

  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}


#' Case with two 6-clusters with gap at each end and between each other.
#' @param seed Initialise the random number generator
two_gappy_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      3, 8,
      13, 18
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:2, c("start", "end"))
  )

  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}


#' Case with five 6-clusters each with gaps between them.
#' @param seed Initialise the random number generator
five_gappy_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      3, 8,
      13, 18,
      23, 28,
      33, 38,
      43, 48
    ),
    nrow = 5,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:5, c("start", "end"))
  )

  return(generate_simple_clusters(
      50,
      30,
      clusters,
      sigma = 1,
      sigma_intraclust = 0.99,
      sigma_interclust = 0.0,
      seed = seed
  ))
}


#' Case of Single incompletely linked cluster in center of window.
#' @param seed Initialise the random number generator
one_noisy_cluster <- function(seed = NULL) {
  clusters <- matrix(
    c(6, 15),
    nrow = 1,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1, c("start", "end"))
  )

  covariance <- matrix(rbeta(20 * 20, shape1 = 0.5, shape2 = 4), nrow = 20, ncol = 20)
  covariance[6:15, 6:15] <- rbeta(10 * 10, shape1 = 4, shape2 = 0.5)
  # Force it to be symmetric
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  diag(covariance) <- 2

  mu <- rnorm(20, mean = 0, sd = 1)
  features <- MASS::mvrnorm(n = 30, mu = mu, Sigma = covariance)

  return(list(
    clusters = clusters,
    covariance = covariance,
    features = features
  ))
}


#' Case of two incompletely linked clusters with random noise.
#' @param seed Initialise the random number generator
two_noisy_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(3, 8,
      13, 18),
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:2, c("start", "end"))
  )

  covariance <- matrix(rbeta(20 * 20, shape1 = 0.5, shape2 = 4), nrow = 20, ncol = 20)
  covariance[3:8, 3:8] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[13:18, 13:18] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)

  # Force it to be symmetric
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  diag(covariance) <- 2


  mu <- rnorm(20, mean = 0, sd = 1)
  features <- MASS::mvrnorm(n = 30, mu = mu, Sigma = covariance)

  return(list(
    clusters = clusters,
    covariance = covariance,
    features = features
  ))
}


#' Case of five incompletely linked clusters with random noise.
#' @param seed Initialise the random number generator
five_noisy_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(3, 8,
      13, 18,
      23, 28,
      33, 38,
      43, 48),
    nrow = 5,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:5, c("start", "end"))
  )

  covariance <- matrix(rbeta(50 * 50, shape1 = 0.1, shape2 = 6), nrow = 50, ncol = 50)
  covariance[3:8, 3:8] <- rbeta(6 * 6, shape1 = 6, shape2 = 0.1)
  covariance[13:18, 13:18] <- rbeta(6 * 6, shape1 = 6, shape2 = 0.1)
  covariance[23:28, 23:28] <- rbeta(6 * 6, shape1 = 6, shape2 = 0.1)
  covariance[33:38, 33:38] <- rbeta(6 * 6, shape1 = 6, shape2 = 0.1)
  covariance[43:48, 43:48] <- rbeta(6 * 6, shape1 = 6, shape2 = 0.1)

  # Force it to be symmetric
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  diag(covariance) <- 2

  mu <- rnorm(50, mean = 0, sd = 1)
  features <- MASS::mvrnorm(n = 50, mu = mu, Sigma = covariance)

  return(list(
    clusters = clusters,
    covariance = covariance,
    features = features
  ))
}


#' Case of five incompletely linked clusters with random noise and some
#' correlation between clusters.
#' @param seed Initialise the random number generator
five_noisy_covarying_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(3, 8,
      13, 18,
      23, 28,
      33, 38,
      43, 48),
    nrow = 5,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:5, c("start", "end"))
  )

  covariance <- matrix(rbeta(50 * 50, shape1 = 0.5, shape2 = 4), nrow = 50, ncol = 50)
  covariance[3:8, 3:8] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[13:18, 13:18] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[23:28, 23:28] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[33:38, 33:38] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[43:48, 43:48] <- rbeta(6 * 6, shape1 = 4, shape2 = 0.5)

  covariance[3:8, 13:18] <-  rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[3:8, 33:38] <-  rbeta(6 * 6, shape1 = 4, shape2 = 0.5)
  covariance[13:18, 33:38] <-  rbeta(6 * 6, shape1 = 4, shape2 = 0.5)

  # Force it to be symmetric
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  diag(covariance) <- 2.5

  heatmap(covariance, Rowv = NA, Colv = NA, symm = TRUE)

  mu <- rnorm(50, mean = 0, sd = 1)
  features <- MASS::mvrnorm(n = 30, mu = mu, Sigma = covariance)

  return(list(
    clusters = clusters,
    covariance = covariance,
    features = features
  ))
}


#' Case with two 6-clusters with gap at each end and between each other.
#' @param seed Initialise the random number generator
two_overlapping_gappy_clusters <- function(seed = NULL) {
  clusters <- matrix(
    c(
      3, 11,
      10, 18
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(1:2, c("start", "end"))
  )

  return(generate_simple_clusters(
      20,
      30,
      clusters,
      sigma = 3,
      sigma_intraclust = 0.9,
      sigma_interclust = 0.0,
      seed = seed
  ))
}
