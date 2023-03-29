#' Generate correlation matrix given a set of gene features.
#'
#' @param data A matrix or dataframe of features for genes.
#'  Features should be in the columns, genes on the rows.
#' @param gene_names [Optional] A vector of gene names to
#'  associate with matrices to assist interpretation later.
#'  if not provided, rownames of data is used.
#'
#' @returns An ngene x ngene matrix of correlation values,
#'  with the lower triangle as NAs.
#' @export
gene_corr <- function(data, gene_names = NULL) {

  if (is.null(gene_names)) {
    gene_names <- rownames(data)
  }

  mat <- data.matrix(data)
  mat <- Hmisc::rcorr(as.matrix(t(mat)))

  # store raw correlation coefficient matrix
  mat <- mat$r
  # Keeping only upper triangle
  # NB, this still keeps the memory allocated
  # Need to look at array based sparse matrix.
  # Scipy dist object-like
  mat[lower.tri(mat)] <- NA

  # Load Gene name Data
  rownames(mat) <- gene_names
  colnames(mat) <- gene_names

  return(mat)
}


#' Finds the number of off diagonal cells in an upper triangle
#' given a grid position.
triangle_size <- function(i, j) {
  offset <- abs(i - j) + 1
  off_diagonal <- offset * (offset - 1)
  upper <- off_diagonal / 2
  return(upper)
}


#' Calculates the sum of correlations for a cell given
#' it's predecessors in the upper triangle of a square matrix.
sum_child_nodes <- function(sums, i, j) {

  # Making this access indices by table would be very slightly faster.
  # [c1, target]
  # [c3 , c2]
  c1 <- sums[i, j - 1]
  c2 <- sums[i + 1, j]

  # We need to subtract c3 to avoid double counting it.
  c3 <- sums[i + 1, j - 1]
  if (is.na(c3) || is.null(c3)) {
    c3 <-  0
  }
  return(c1 + c2 - c3)
}


#' Calculates the positions of the upper triangle in diagonal order.
#' Intended for use with iteration.
#'
#' @param index The index to find a grid coordinate for.
#' @param size The size of the matrix you're iterating over.
#'
#' @returns A list with the grid coordinates (i, j), and the
#'  distance from the center diagonal (k).
index_to_diag_ij <- function(index, size, include_diag = FALSE) {
  if (include_diag && (index <= size)) {
    return(list(i = index, j = index, k = 0))
  } else if (include_diag) {
    index <- index - size
  }
  i <- size - 1 - floor(sqrt(- 8 * index + 4 * size * (size - 1) + 1) / 2 - 0.5)
  j <- index + i + ((size - i + 1) * (size - i) - size * (size - 1)) / 2
  i <- j - i
  k <- abs(i - j)
  return(list(i = i, j = j, k = k))
}


#' Generate average correlation matrix.
#'
#' For each off diagonal cell in the upper triangle of a square matrix
#' calculates the average of all cells in the triangle below
#' (i.e. sum(mat[i:j, i:j]) if the diagonal and lower triangle are all 0).
#'
#' @param data An n x n matrix of correlation values.
#'  The lower triangle should be set to NA (it will raise an error if not).
#' @param bandwith [Optional] The maximum off-diagonal distance to calculate the
#'  averages for. This effectively sets the maximum cluster size but can
#'  significantly speed up calculation.
#'
#' @returns An n x n matrix with the triangle averages.
#' @export
average_corr <- function(data, bandwidth = NULL) {

  averages <- sums <- matrix(
    nrow = nrow(data),
    ncol = nrow(data)
  )

  diag(data) <- NA
  size <- nrow(data)
  n_cells <- (size * (size - 1)) / 2

  if (is.null(bandwidth)) {
    bandwidth <- size
  }

  for (h in seq_len(n_cells)) {
    indices <- index_to_diag_ij(h, nrow(data))
    i <- indices$i
    j <- indices$j
    k <- indices$k

    if (k > bandwidth) {
      break
    }

    this_corr <- data[i, j]

    # On the first diagonal, we just want the raw values.
    if (k == 1) {
      this_sum <- this_corr
    } else {
      this_sum <- sum_child_nodes(sums, i, j) + this_corr
    }

    sums[i, j] <- this_sum
    averages[i, j] <- this_sum / triangle_size(i, j)
  }
  return(averages)
}


#' Creates matrix with 1 where correlation exceeds thresholds and 0 otherwise
threshold_matrix <- function(data, threshold) {
  thresholded <- apply(
    data >= threshold,
    MARGIN = 2,
    FUN = as.numeric
  )
  thresholded[is.na(thresholded)] <- 0
  return(thresholded)
}


#' Rotates a matrix -90 degrees.
#' Contrasts with t() which rotates 90 degrees.
t2 <- function(mat) {
  return(apply(t(mat), 2, rev))
}


#' Check if a matrix is square.
is_square <- function(mat) {
  return(nrow(mat) != ncol(mat))
}


#' Find a square matrix's size
#' Asserts that the matrix is square.
find_matrix_size <- function(mat) {
  if (!is_square(mat)) {
    stop("Matrix must be square.")
  }

  return(nrow(mat))
}


#' Find the coordinates of diagonals in a square matrix.
#' Can only get upper diagonal coordinates.
#'
#' @param mat the matrix that you want to get coords from.
#'  This is only used to get the size of the matrix.
#' @param offset the diagonal to get. 0 is the main diagonal,
#'  1 will be the adjacent n-1 diagonal etc.
#'
#' @returns An n-offset by 2 matrix of coordinates for the
#'  diagonal positions.
diag_indices <- function(mat, offset = 0) {
  nrows <- nrow(mat)

  if (nrows <= 1) {
    stop("Matrix is too small.")
  }

  if (offset < 0) {
    stop("We only support positive offsets.")
  }

  len <- nrows - offset

  indices <- matrix(0, nrow = len, ncol = 2)
  indices[, 1] <- 1:len
  indices[, 2] <- indices[, 1] + offset
  return(indices)
}


#' Finds gaps in the first offset diagonal of a matrix.
#' @param mat The matrix to find gaps in. This must be a binary indicator
#'  matrix where 1 says the cell has passed a threshold. We don't check that
#'  this is true so make sure it is in the function passing the matrix.
#' @returns An (n-1) vector with gaps marked as 1 and non-gaps as 0.
find_gaps <- function(mat) {
  out <- (1 - mat[diag_indices(mat, offset=1)])
  return(out)
}


#' Finds the scanning cumulative sum of a vector, resetting the
#' counter whenever it encounters a 0.
#' @returns An (n-1) vector with the cumulative sum of gaps since the last
#'  non-gap. E.g. 0, 1, 2, 0 means that the second and third positions are gaps
#'  and the first and 4th are non-gaps.
resetting_cumsum <- function(gaps) {
  Reduce(
    function(i, j) ifelse(j == 0, 0, i + j),
    x = gaps,
    accumulate = TRUE
  )
}


#' Given a vector from resetting_cumsum returns the largest
#' gap for a slice (i..j) of the vector.
#'
#' @returns An integer
find_max_gap <- function(gaps, i, j) {
  vec <- gaps[i:(j - 1)]
  return(max(vec))
}


#' Find the edge score for a given triangle in a square matrix.
#' The edge score is the minimum sum of cells on the outer row or column
#' of the triangle where the cell passed the threshold.
#'
#' @param mat A binary indicator matrix, where 1 indicates that the cell passed
#'  the threshold. Note that only the upper triangle is used.
#' @param bandwidth The maximum triangle size that will be considered.
#'  Used to reduce run-time and memory requirements.
#'
#' @returns A matrix with the upper triangle populated with edge scores for
#'  each potential cluster.
find_edges <- function(mat, bandwidth = NULL) {
  sums_i <- matrix(0, nrow(mat), ncol(mat))
  sums_j <- matrix(0, nrow(mat), ncol(mat))
  mins <- matrix(0, nrow(mat), ncol(mat))

  size <- nrow(mat)
  n_cells <- (size * (size - 1)) / 2

  if (is.null(bandwidth) || (bandwidth > size)) {
    b_cells <- 0
    bandwidth <- size
  } else {
    b_cells <- ((size - bandwidth) * (size - bandwidth - 1)) / 2
  }
  for (h in seq_len(n_cells - b_cells)) {
    indices <- index_to_diag_ij(h, nrow(mat))

    this <- mat[indices$i, indices$j]
    if (indices$k <= 1) {
      child_i <- child_j <- 0
    } else if (indices$k >= bandwidth) {
      break
    } else {
      child_j <- sums_j[indices$i, indices$j - 1]
      child_i <- sums_i[indices$i + 1, indices$j]
    }

    sums_j[indices$i, indices$j] <- child_j + this
    sums_i[indices$i, indices$j] <- child_i + this
    mins[indices$i, indices$j] <- min(c(child_j, child_i)) + this
  }
  return(mins)
}



#' Find clusters satisfying threshold criteria.
#'
#' @param edge_mins
#' @param surface_means The frequency within triangle of averaged
#'  elements passing triangle (from average_corr on the filtered matrix).
#' @param edge_mins For each position, the minimum edge score (from find_edges).
#' @param edge_fill_rate The threshold that the edge_mins[i, j] for a cluster
#'  must pass.
#' @param gap_penalty The penalty for the maximum gap size in the cluster.
#' @param surface_fill_rate The minimum frequency of cells in the cluster
#'  triangle passing the threshold.
#' @param min_cluster The minimum cluster size (an integer).
#' @param max_cluster [Optional] The maximum cluster size (an integer).
#'  Helpful to speed up. If not specified, there is no limit on cluster size.
#'
#' @returns A dataframe containing the positions of clusters and their scores.
#' @export
find_threshold_clusters <- function(
  adjacency,
  edge_mins,
  surface_means,
  gaps,
  edge_fill_rate = 0.5,
  gap_penalty = 0.5,
  surface_fill_rate = 0.5,
  min_cluster = 4,
  max_cluster = NULL
) {

  # TODO: split this out into more sub-functions.
  #       It's not too hard to understand but may be hard to test.

  size <- nrow(edge_mins)
  if (is.null(max_cluster)) {
    max_cluster <- size - 1
  }

  offset_max <- min(c(size - 1, max_cluster))
  offset_min <- max(c(min_cluster, 2))

  if (min_cluster > max_cluster) {
    stop("max_cluster must be greater than or equal to min_cluster.")
  }

  cluster_index <- 1
  clusters <- list()

  skip <- logical(size)

  for (k in offset_max:offset_min) {
    indices <- diag_indices(edge_mins, k)
    for (h in seq_len(nrow(indices))) {
      i <- indices[h, 1]
      j <- indices[h, 2]

      # This means that a larger cluster has already included this range.
      if (all(skip[j:i])) {
        next
      }

      surface_mean <- surface_means[i, j]

      if (surface_means[i, j] < surface_fill_rate) {
        next
      }


      edge_score <- edge_mins[i, j] / (k + 1)
      max_gap <- find_max_gap(gaps, i, j)

      if (edge_score < edge_fill_rate) {
        next
      } else if (max_gap >= (gap_penalty * (1 - edge_fill_rate) * (k + 1))) {
        next
      }

      # The average on edges less than threshold
      edge_i_nok <- mean(adjacency[i:(j - 1), j], na.rm = TRUE)
      edge_i_nok <- edge_i_nok < surface_fill_rate * edge_fill_rate
      edge_j_nok <- mean(adjacency[i, (i + 1):j], na.rm = TRUE)
      edge_j_nok <- edge_j_nok < surface_fill_rate * edge_fill_rate

      if (edge_i_nok || edge_j_nok) {
        next
      }

      clusters[[cluster_index]] <- as.data.frame(t(c(
        start = i,
        end = j,
        edge_score = edge_score,
        max_gap = max_gap,
        surface_mean = surface_mean,
        size = k + 1,
        merged = FALSE
      )))
      cluster_index <- cluster_index + 1

      skip[i:j] <- TRUE
    }
  }

  clusters <- do.call(rbind, clusters)

  if (is.null(clusters)) {
    # If there are no clusters, do.call will return NULL
    # It should return an empty table with the right columns

    cols <- c(
      "start", "end", "edge_score",
      "max_gap", "surface_mean", "size", "merged"
    )
    clusters <- setNames(
      data.frame(matrix(ncol = length(cols), nrow = 0)),
      cols
    )
  } else {
    clusters <- clusters[order(clusters$start, clusters$end), ]
    rownames(clusters) <- NULL
  }

  return(clusters)
}


#' Calculates the frequency that an intersecting triangle of two
#' clusters passes the threshold.
#' Used to decide if cluster should be merged.
#'
#' @param l The left cluster row (from find_clusters).
#' @param r The right ...
#' @thresholded The raw correlation matrix (not averaged) with
#'  threshold applied.
#'
#' @returns A float.
intersect_score <- function(l, r, thresholded) {
  m <- thresholded[r$start:l$end, r$start:l$end]
  m[lower.tri(m)] <- NA
  return(mean(m, na.rm = TRUE))
}


#' Calculates the frequency that the square of non-overlapping parts
#' between two clusters passes the threshold.
#' Used to decide if cluster should be merged.
#'
#' @param l The left cluster row (from find_clusters).
#' @param r The right ...
#' @thresholded The raw correlation matrix (not averaged) with
#'  threshold applied.
#'
#' @returns A float.
union_diff_score <- function(l, r, thresholded) {
  m <- thresholded[l$start:(r$start - 1), (l$end + 1):r$end]
  return(mean(1 - m, na.rm = TRUE))
}


#' Merge overlapping clusters if they pass thresholds.
#'
#' @param clusters The dataframe returned from find_clusters.
#' @param thresholded The raw correlation matrix (not averaged) with
#'  threshold applied.
#' @param surface_means The frequency within triangle of averaged
#'  elements passing triangle (from average_corr on the filtered matrix).
#' @param edge_mins For each position, the minimum edge score (from find_edges).
#' @param gaps A vector from find_gaps.
#' @param merge_threshold The threshold for cluster merging.
#'  Cluster merges are scored from "thresholed" as
#'  (frequency intersecting) / (frequency 1-not-intersecting).
#'  Basically, the numerator is the score to merge, the denominator is the
#'  score to keep separate. A threshold of 1 means that if this ratio is
#'  1 or more, the clusters will be merged.
#'
#' @returns An updated dataframe containing positions and scores of clusters.
#' @export
merge_threshold_clusters <- function(
  clusters,
  adjacency,
  surface_means,
  edge_mins,
  gaps,
  merge_threshold = 0.9
) {
  if (nrow(clusters) <= 1) {
    # If only one cluster, cannot possibly merge
    return(clusters)
  }

  clusters <- clusters[order(clusters$start, clusters$end), ]

  new_clusters <- list()
  h <- 1
  i <- 2

  lag <- clusters[1, ]

  while (i <= nrow(clusters)) {
    cl <- clusters[i, ]
    if (lag$end < cl$start) {
      new_clusters[[h]] <- lag
      h <- h + 1
      lag <- cl
      i <- i + 1
      next
    }

    iscore <- intersect_score(lag, cl, adjacency)
    uscore <- union_diff_score(lag, cl, adjacency)

    if ((iscore / uscore) >= merge_threshold) {
      lag$end <- cl$end
      lag$max_gap <- find_max_gap(gaps, lag$start, lag$end)
      k <- abs(lag$start - lag$end) + 1
      lag$edge_score <- edge_mins[lag$start, lag$end] / k
      lag$surface_mean <- surface_means[lag$start, lag$end]
      lag$size <- k
      lag$merged <- TRUE
    } else {
      new_clusters[[h]] <- lag
      h <- h + 1
      lag <- cl
    }

    i <- i + 1
  }

  new_clusters[[h]] <- lag
  new_clusters <- do.call(rbind, new_clusters)

  if (is.null(new_clusters)) {
    stop(paste(
      "Somehow we've got zero clusters after merging.",
      "This shouldn't be possible."
    ))
  }

  new_clusters <- new_clusters[order(new_clusters$start, new_clusters$end), ]
  rownames(new_clusters) <- NULL
  return(new_clusters)
}


#' Cluster an adjacency matrix using the thresholding method.
#'
#' @param adjacency A matrix of 0 or 1 where 1 indicates an
#'  between two genes.
#' @param threshold = 0.8
#' @param edge_fill_rate = 0.8
#' @param surface_fill_rate = 0.6
#' @param gap_penalty = 0.5
#' @param min_cluster = 4
#' @param max_cluster = NULL
#' @param smooth = "average"
#' @param merge = TRUE
#' @param merge_threshold = 1.1
#' @param verbose = FALSE
#'
#' @returns
#'
#' @export
threshold_clustering <- function(
  adjacency,
  threshold = 0.8,
  edge_fill_rate = 0.8,
  surface_fill_rate = 0.6,
  gap_penalty = 0.5,
  min_cluster = 4,
  max_cluster = NULL,
  smooth = "average",
  merge = TRUE,
  merge_threshold = 1.1,
  verbose = FALSE
) {
  # We want to avoid mutating anything.
  adjacency <- adjacency

  if (verbose && (smooth %in% c("average", "TOM"))) {
    print(sprintf("- Running %s smoothing", smooth))
  }
  if (smooth == "average") {
    smoothed <- average_corr(adjacency, bandwidth = max_cluster)
  } else if (smooth == "TOM") {
    smoothed <- WGCNA::TOMsimilarity(adjacency)
  } else {
    smoothed <- adjacency
  }

  adjacency[lower.tri(adjacency)] <- NA
  smoothed[lower.tri(smoothed)] <- NA

  smoothed <- threshold_matrix(smoothed, threshold = threshold)

  if (verbose) {
    print("- Finding edges")
  }
  edge_mins <- find_edges(smoothed, bandwidth = max_cluster)

  if (verbose) {
    print("- Finding surface means")
  }
  surface_means <- average_corr(smoothed, bandwidth = max_cluster)
  gaps <- find_gaps(smoothed)

  if (verbose) {
    print("- Finding clusters")
  }
  clusters1 <- find_threshold_clusters(
    adjacency = adjacency, 
    edge_mins = edge_mins,
    surface_means = surface_means,
    gaps = gaps,
    edge_fill_rate = edge_fill_rate,
    gap_penalty = gap_penalty,
    surface_fill_rate = surface_fill_rate,
    min_cluster = min_cluster,
    max_cluster = max_cluster
  )

  if (merge) {
    if (verbose) {
      print("- Merging clusters")
    }
    clusters2 <- merge_threshold_clusters(
      clusters = clusters1,
      adjacency = adjacency,
      surface_means = surface_means,
      edge_mins = edge_mins,
      gaps = gaps,
      merge_threshold = merge_threshold
    )
  }

  out <- list(
    clusters = clusters2,
    smoothed = smoothed,
    surface_means = surface_means,
    edge_mins = edge_mins
  )
  return(out)
}
