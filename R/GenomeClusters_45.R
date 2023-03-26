#' Generate correlation matrix
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
#' @param index The index to find a grid coordinate for.
#' @param size The size of the matrix you're iterating over.
#'
#' Returns a list with the grid coordinates (i, j), and the
#' distance from the center diagonal (k).
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
#' @param mat the matrix that you want to get coords from.
#'  This is only used to get the size of the matrix.
#' @param offset the diagonal to get. 0 is the main diagonal,
#'  1 will be the adjacent n-1 diagonal etc.
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
  indices[,1] <- 1:len
  indices[,2] <- indices[,1] + offset
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
    function(i, j) ifelse(j==0, 0, i + j),
    x = gaps,
    accumulate = TRUE
  )
}


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

  if (is.null(bandwidth)) {
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


#' Cluster detection (2) - cluster surface check function:
#' @export
calc_surface_sum <- function(data, i, j) {

  #checking number of cells in triangle above threshold
  surface_sum <- 0
  size <- sqrt(length(data))

  for (e in j:(size - i)) {
    for (e2 in i:(size - j)) {
      surface_sum <- surface_sum + data[e2, e]
    }
  }
  return(surface_sum)
}


#' Cluster detection (3) - cell skipping function:
#' @export
skip_matrix_cell <- function(
  threshold,
  data,
  i,
  j,
  min_clust_size
) {
  size <- sqrt(length(data))
  if (data[i, j] < threshold) {
    # skip matrix cells with correlation below threshold
    skip <- TRUE
  } else if ((size - i - j) < (min_clust_size - 1)) {
    # skip matrix cells too close to give a cluster larger than min_clust_size
    skip <- TRUE
  } else if ((i + min_clust_size) > (size - 1)) {
    # skip matrix cells too close to end of contig
    skip <- TRUE
  } else {
    skip <- FALSE
  }
  return(skip)
}


#' Cluster detection (4) - initial cluster rendering function
#' WARNING! Mutates matrix
#' @export
gen_clust_matrix <- function(data, i, j) {

  size <- sqrt(length(data))

  # Generates cluster matrix (data); scan columns
  for (e in i:(size - j)) {

    # Generates cluster matrix (data); scan lines
    for (e2 in i:(size - j)) {
      # TODO: find non-mutating version if possible
      data[e2, e] <- 1
    }
  }
  data[i, (size - j)] <- 0.5
  return(data)
}


#' Cluster detection (5) - initial cluster detection function
#' @export
first_pass_clusters <- function(
  filtered,
  clustered,
  gap_penalty,
  edge_fill_rate,
  surface_fill_rate,
  min_clust_size,
  threshold
) {
  # Initialize variables
  # number of genes in contig = size of the correlation matrix
  size <- sqrt(length(filtered))
  ikeep <- 0 # store line position when cluster found
  jkeep <- 0 # store column position when cluster found
  # store number of cells above corelation threshold in a putative cluster
  surfacesum <- 0
  i <- 0
  j <- 0
  surfacemax <- 0

  # counter scanning lines
  for (i in 2:(size)) {
    # counter scanning columns
    for (j in 1:(size - (i - 1))) {
      skip_cell <- skip_matrix_cell(
        threshold,
        filtered,
        i,
        j,
        min_clust_size
      )

      if (skip_cell) {
        next
      }

      edges <- calc_edges_and_gaps(filtered, i, j)

      # edges and gaps checking
      if (
        (edges[1] > (edge_fill_rate * (size - i - j))) &&
        (edges[2] < (gap_penalty * (1 - edge_fill_rate) * (size - i - j)))
      ) {
        next
      }

      surfacesum <- calc_surface_sum(filtered, i, j)
      surfacemax <- 0

      # calculates number of cells in cluster triangle
      for (e in 1:(1 + size - i - j)) {
        surfacemax <- surfacemax + e
      }

      if (surfacesum > (surface_fill_rate * surfacemax)) {
        next
      }
      ikeep <- i
      jkeep <- j
      clustered <- gen_clust_matrix(
        clustered,
        ikeep,
        jkeep
      )
    }
  }
  return(clustered)
}


#' Cluster detection (6) - join clusters & clean up function
#' @export
merge_clusters <- function(clustered, ave_corr, threshold) {
  size <- sqrt(length(ave_corr))
  clustered <- t2(clustered)
  i <- 0
  j <- 0

  # Cells of the diagonal within a cluster turned to
  # 1's = join overlaping and contiguous clusters
  for (i in seq_len(size - 1)) {
    if (
      (clustered[i + 1, i] == 1) &&
      (ave_corr[i, i + 1] > threshold)
    ) {
      clustered[i, i] <- 1
    }
  }

  # Cells of the diagonal at the edge of cluster with correlation
  # threshold turned to 0's
  stop_cleaning <- FALSE
  count_rounds <- 0
  repeat {
    count_rounds <- count_rounds + 1
    cluster_cleaned <- FALSE
    for (i in 1:(size - 1)) {
      if (
        (clustered[i, i] == 1) &&
        (clustered[i + 1, i + 1] == 0) &&
        (ave_corr[i, i + 1] < threshold)
      ) {
        clustered[i, i] <- 0
        cluster_cleaned <- TRUE
      }
    }
    for (j in 2:size) {
      if (
        (clustered[j, j] == 1) &&
        (clustered[j - 1, j - 1] == 0) &&
        (ave_corr[j - 1, j] < threshold)
      ) {
        clustered[j, j] <- 0
        cluster_cleaned <- TRUE
      }
    }
    if (cluster_cleaned == FALSE) {
      stop_cleaning <- TRUE
    } else {
      next
    }
    if (stop_cleaning == TRUE) {
      break
    }
  }
  return(clustered)
}


#' Cluster detection (7) - finalize cluster rendering function
#' @export
finalize_cluster_matrix <- function(clustered, start, end) {
  clustered[start:end, start:end] <- 1
  return(clustered)
}


#' Cluster detection (8) - Clean up and scan final clusters function
#' @export
second_pass_clusters <- function(
  clustered,
  min_clust_size
) {
  cluster_found <- FALSE # boolean storing position in/out of cluster
  cluster_counter <- 0 # counter for the total number of clusters in contig
  # counter for the total number of genes in cluster in contig
  clustered_genes <- 0
  clustered_all <- 0
  start_i <- 0
  end_i <- 0
  i <- 0
  j <- 0
  size <- sqrt(length(clustered))

  for (i in seq_len(size)) {
    if (
      (clustered[i, i] == 1) &&
      (!cluster_found)
    ) { # first gene of a cluster
      cluster_found <- TRUE
      clustered_genes <- 1
      start_i <- i - 1
    } else if (
      (clustered[i, i] == 0) &&
      cluster_found &&
      ((i + 1 - start_i) > min_clust_size)
    ) { # reached end of cluster
      cluster_found <- FALSE
      cluster_counter <- cluster_counter + 1
      clustered_genes <- i - start_i
      clustered_all <- clustered_all + clustered_genes
      end_i <- i - 1
      clustered <- finalize_cluster_matrix(clustered, start_i, end_i)

      j <- i
      repeat {
        clustered[j, i] <- 0
        j <- j + 1
        if (j == size + 1) {
          break
        }
      }
      k <- i
      repeat {
        clustered[i, k] <- 0
        k <- k - 1
        if (k < 2) {
          break
        }
      }
    } else if (
      (i == size) &&
      (cluster_found == TRUE) &&
      ((i + 1 - start_i) > min_clust_size)
    ) { # reched end of contig within a cluster
      cluster_found <- FALSE
      cluster_counter <- cluster_counter + 1
      clustered_genes <- (i - start_i)
      clustered_all <- clustered_all + clustered_genes
      end_i <- i
      clustered <- finalize_cluster_matrix(clustered, start_i, end_i)

      j <- i

      repeat {
        clustered[j, i] <- 0
        j <- j + 1

        if (j == size + 1) {
          break
        }
      }
      k <- i
      repeat {
        clustered[i, k] <- 0
        k <- k - 1
        if (k < 2) {
          break
        }
      }
    } else if (
      (clustered[i, i] == 0) &&
      (cluster_found == TRUE) &&
      ((i + 1 - start_i) < min_clust_size)
    ) {
      # cluster too small = clean up matrix with 0's
      cluster_found <- FALSE

      for (e in start_i:i) {
        j <- e
        repeat {
          clustered[j, e] <- 0
          j <- j + 1
          if (j == size + 1) {
            break
          }
        }
        k <- e
        repeat {
          clustered[e, k] <- 0
          k <- k - 1
          if (k < 2) {
            break
          }
        }
      }
    } else {
      # out of a cluster or if cluster too small, clean up cluster matrix with 0's
      j <- i
      repeat {
        clustered[j, i] <- 0
        j <- j + 1
        if (j == size + 1) {
          break
        }
      }
      k <- i
      repeat {
        clustered[i, k] <- 0
        k <- k - 1
        if (k < 2) {
          break
        }
      }
    }
  }
  return(clustered)
}


#' Run the master clustering pipeline
#'
#' @param MyData The dataset loaded from input file.
#' @param AnalysisName The name of the analysis/dataset.
#' @param MyCut correlation threshold [0 to 1] from the averaged matrix
#'  (see [*1])
#' @param UseQuantile, boolean, whether to use quantile (TRUE) or corr coeff
#'  (FALSE) of averaged matrix as threshold
#' @param Edge.Fill.rate, % of cells above correlation threshold in the edge
#'  of a submatrix to proceed further in cluster detection [0 to 1]
#' @param Edge.Fill.rate Defines the % of matrix cells in the horizontal and
#'  vertical edges of cluster triangles that need to be above the correlation
#'  threshold. Numeric between 0 (any number of cells) and 1 (all cells). Self
#'  correlations are taken into account.
#' @param Surface.Fill.rate, % of cells above correlation threshold in a
#'  submatrix call a cluster [0 to 0.5]
#' @param Surface.Fill.rate Defines the % of the cells of the cluster triangle
#'  that must be above the correlation threshold. Numeric between 0
#'  (any number of cells) and 1 (all cells). Self correlations are ignored.
#' @param Min.Clus.size, Minimum number of genes in a cluster during first
#'  pass detection
#' @param Min.Clus.size defines the minimal size above which clusters are
#'  detected. This parameter is used both on first and second pass detection,
#'  meaning that the final list of clusters does not necessarily reflect how
#'  they were initially detected. As a result you might find clusters of
#'  size eg 20 with a small Min.Clus.size value, while no cluster will be
#'  detected if Min.Clus.size is set to 20.
#' @param Gap.penalty, Reduces the relative max gap length [positive] (see [*2])
#' @param Gap.penalty Sets the % of matrix cells in the diagonal edge of
#'  cluster triangles (not counting self-correlations) that need to be above
#'  the correlation threshold. 1 is neutral, numeric between 0 and 1 will relax
#'  threshold, numeric above 1 will increase stringency.
#' @param Report.type none, mini, genelist, normal or verbose
#' @returns
#'
#' @examples
#' # todo
#'
#' @export
detect_clusters <- function(
  data,
  threshold,
  use_quantile,
  edge_fill_rate,
  surface_fill_rate,
  min_clust_size,
  gap_penalty,
  gene_names = NULL
) {
  gap_penalty <- (1 / gap_penalty)

  # Running cluster detection
  correlated <- gene_corr(data, gene_names)
  size <- sqrt(length(correlated))

  averaged <- average_corr(correlated)

  if (use_quantile) {
    threshold <- quantile(averaged, threshold, na.rm = TRUE)
  }

  filtered <- threshold_matrix(averaged, threshold)

  # Unsure why this is necessary?
  # rotate matrix 90? counter clockwise to facilitate cluster detection
  filtered <- t2(filtered)

  # Matrix in which the position of detected clusters will be stored in
  # the form of 1s
  clustered <- matrix(0, nrow(filtered), ncol(filtered))

  clustered <- first_pass_clusters(
    filtered,
    clustered,
    gap_penalty,
    edge_fill_rate,
    surface_fill_rate,
    min_clust_size,
    threshold
  )

  clustered <- merge_clusters(clustered, averaged, threshold)

  clustered <- second_pass_clusters(
    clustered,
    min_clust_size
  )

  cluster_counter <- 0
  for (i in 1:(size - 1)) {
    if (
      (clustered[i, i] == 0) &&
      (clustered[i + 1, i + 1] == 1)
    ) {
      cluster_counter <- cluster_counter + 1
    }
    if (clustered[1, 1] == 1) {
      cluster_counter <- cluster_counter + 1
    }
  }

  # Combine averaged matrix with detected cluster matrix and save as files
  final <- t(averaged) + 0.3 * clustered
  final[clustered == 0] <- 0
  diag(final) <- diag(clustered)

  final[upper.tri(final)] <- correlated[upper.tri(correlated)]
}


#' Find clusters satisfying threshold criteria.
#'
#' @param edge_mins
#' @param surface_means
#' @param gaps
#' @param edge_fill_rate
#' @param gap_penalty
#' @param surface_fill_rate
#' @param min_cluster
#' @param max_cluster
#'
#' @returns
find_clusters <- function(
  edge_mins,
  surface_means,
  gaps,
  edge_fill_rate = 0.5,
  gap_penalty = 0.5,
  surface_fill_rate = 0.5,
  min_cluster = 4,
  max_cluster = NULL
) {
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
  #clusters <- matrix(FALSE, nrow = nrow(edge_mins), ncol = ncol(edge_mins))

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

      edge_score <- edge_mins[i, j] / (k + 1)
      max_gap <- find_max_gap(gaps, i, j)

      if (edge_score < edge_fill_rate) {
        next
      } else if (max_gap >= (gap_penalty * (1 - edge_fill_rate) * (k + 1)) ) {
        next
      }

      surface_mean <- surface_means[i, j]

      if (surface_means[i, j] < surface_fill_rate) {
        next
      }

      clusters[[cluster_index]] <- as.data.frame(t(c(
        start=i,
        end=j,
        edge_score=edge_score,
        max_gap=max_gap,
        surface_mean=surface_mean,
        size=k + 1
      )))
      cluster_index = cluster_index + 1

      skip[i:j] <- TRUE
    }
  }
  clusters <- do.call(rbind, clusters)
  return(clusters)
}
