###########################################################################
# This file contains a series of R functions for analyzing co-varying
# gene clusters in genomes.
#
# Copyright (c) 2021 Sylvain Raffaele
#
# This is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# Sylvain Raffaele
# sylvain.raffaele@inrae.fr
# Laboratoire des Interactions Plantes-Microbes-Environnement
# 31326 Castanet Tolosan, France
#
# version 0.2.29 (Jul 2021)
###########################################################################

###########################################################################
# USAGE
#
# Load required libraries:
# library('rlang')
# library('ggplot2')
# library("Hmisc")
# library('gplots')
#
# Define your favorite color palette:
# mypalette<-colorRampPalette(c("darkblue", "blue", "cornflowerblue", "lightblue", "white", "white", "goldenrod1", "orangered", "darkred"), space="rgb")
#
# If needed, set working directory containing your input files
# setwd("C:/MyWorkingDir/R")
#
# Load this set of functions by typing in:
# >source("GenomeClusters.R")
#
# Load your input file containing expression data under the name "MyData"
# MyData<-read.table("MyDataFile.tab", sep="\t", header=F)
#
# Load your input file containing ordered list of genes on contig
# MyGeneNames<-scan("MyData_Genelist.txt", what="", sep="\n")
#
# Name your analysis by typing in:
# AnalysisName<-"MyAnalysis"
#
# Run the master function by typing in:
# >RunClusterDetection(MyData, AnalysisName, MyCut=xx, Edge.Fill.rate=xx, Surface.Fill.rate=xx, Min.Clus.size=x, Gap.penalty=x, verbose=x)
#
# DESCRIPTION OF PARAMETERS
#
# MyData passes on the dataset loaded from input file
#
# AnalysisName passes on the name defined by previous instruction
#
# Edge.Fill.rate defines the % of matrix cells in the horizontal and
# vertical edges of cluster triangles that need to be above the correlation
# threshold. Numeric between 0 (any number of cells) and 1 (all cells). Self
# correlations are taken into account.
#
# Surface.Fill.rate defines the % of the cells of the cluster triangle
# that must be above the correlation threshold. Numeric between 0
# (any number of cells) and 1 (all cells). Self correlations are ignored.
#
# Min.Clus.size defines the minimal size above which clusters are detected.
# This parameter is used both on first and second pass detection, meaning
# that the final list of clusters does not necessarily reflect how they
# were initially detected. As a result you might find clusters of size eg 20
# with a small Min.Clus.size value, while no cluster will be detected if
# Min.Clus.size is set to 20.
#
# Gap.penalty is used to set the % of matrix cells in the diagonal edge of
# cluster triangles (not counting self-correlations) that need to be above
# the correlation threshold. 1 is neutral, numeric between 0 and 1 will relax
# threshold, numeric above 1 will increase stringency.
#
# verbose: whether of not to print additional information in the txt report
# and save intermediate matrices, usually for debugging. Either TRUE or FALSE
#
###########################################################################

###########################################################################
# Generate pairwise correlation matrix function :
###########################################################################

RawCorMat <- function (MyData, MyGeneNames) {

  mat.MyData <- data.matrix(MyData)
  MycorMyData <- rcorr(as.matrix(t(mat.MyData)))
  rmatMyData <- MycorMyData$r # store raw correlation coefficient matrix
  rmatMyData[lower.tri(rmatMyData)] <- NA # Keeping only upper triangle

  rownames(rmatMyData) <- MyGeneNames # Load Gene name Data
  colnames(rmatMyData) <- MyGeneNames

  return(rmatMyData)
}

###########################################################################
# Generate local average matrix (suffix 'A') function:
###########################################################################

AveCorMat <- function (rmatMyData) {

  rmatMyDataA <- SumMatrix <- CountMatrix <- matrix(
    nrow=nrow(rmatMyData),
    ncol=nrow(rmatMyData)
  )

  for (i in 1:nrow(rmatMyData)) {
    rmatMyData[i, i] <- NA
  }

  for (i in 1:nrow(rmatMyDataA)) {
    for (j in 1:nrow(rmatMyDataA)) {
      if (i>=j) {
        next
      }
      TempMatrix <- rmatMyData[min(i, j):max(i, j), min(i, j):max(i, j)]
      CountMatrix[i, j] <- sum(!is.na(TempMatrix))
      SumMatrix[i, j] <- sum(TempMatrix, na.rm = TRUE)
      rmatMyDataA[i, j] <- SumMatrix[i, j] / CountMatrix[i, j]
    }
  }
  rmatMyDataA[lower.tri(rmatMyDataA)] <- NA # Keeping only upper triangle

  # erase self correlations
  for (i in 1:nrow(rmatMyDataA)) {
    rmatMyDataA[i, i] <- NA
  }

  return(rmatMyDataA)
}

###########################################################################
# Cluster detection (1) - cluster edges and gaps function:
###########################################################################

# takes rmatMyDatafilt matrix and i and j coordinates as input
CalculateEdgesAndGaps <- function (rmatMyDatafilt, i, j) {

  e <- 1
  MySize <- sqrt(length(rmatMyDatafilt))
  # checking % of cells at the edge of cluster triangle above correlation threshold
  horizedge <- 0
  vertedge <- 0
  diaggap <- 0
  d <- 0
  for (e in 0:(MySize - (i + j))) {
    horizedge <- horizedge + rmatMyDatafilt[i, j + e]
    vertedge <- vertedge + rmatMyDatafilt[i + e, j]
    if (rmatMyDatafilt[i + e, (MySize - i - e)] == 1) {
      d <- 0
    } else {
      d <- d + 1
    }
    diaggap <- c(diaggap, d)
  }

  # only diag considered for gap checking
  EdgesAndGaps <- c(min(horizedge, vertedge), max(diaggap))
  return(EdgesAndGaps)
}

###########################################################################
# Cluster detection (2) - cluster surface check function:
###########################################################################

CalculateSurfaceSum <- function(rmatMyDatafilt, i, j) {

  #checking number of cells in triangle above threshold
  SurfaceSum <- 0
  MySize <- sqrt(length(rmatMyDatafilt))

  for (e in j:(MySize - i)) {
    for (e2 in i:(MySize - j)) {
      SurfaceSum <- SurfaceSum + rmatMyDatafilt[e2, e]
    }
  }
  return(SurfaceSum)
}


###########################################################################
# Cluster detection (3) - cell skiping function:
###########################################################################

SkipMatrixCell <- function(
  corr.threshold,
  rmatMyDatafilt,
  i,
  j,
  ikeep,
  jkeep,
  Min.Clus.size,
  report_name,
  Report.type
) {
  MySize <- sqrt(length(rmatMyDatafilt))
  # skip matrix cells with correlation below threshold
  if (rmatMyDatafilt[i, j] < corr.threshold) {
    SkipCell <- TRUE
    if (Report.type == "verbose") {
      cat(
        paste("skip: low cor. ", i, ",", j),
        file = report_name,
        sep = "\n",
        append = TRUE
      )
    }
  } else {
    # skip matrix cells too close to give a cluster larger than Min.Clus.size
    if ((MySize - i - j) < (Min.Clus.size - 1)) {
      SkipCell <- TRUE
      if (Report.type == "verbose") {
        cat(
          paste("skip: too few genes", i, ",", j),
          file = report_name,
          sep = "\n",
          append = TRUE
        )
      }
    } else {
      # skip matrix cells too close to end of contig
      if ((i + Min.Clus.size) > (MySize - 1)) {
        SkipCell <- TRUE
        if (Report.type == "verbose") {
          cat(
            paste("skip: too close to chr. end ", i, ",", j),
            file = report_name,
            sep = "\n",
            append = TRUE
          )
        }
      } else {
        SkipCell <- FALSE
      }
    }
  }
  return(SkipCell)
}


###########################################################################
# Cluster detection (4) - initial cluster rendering function:
###########################################################################

GenerateClusterMatrix <- function(rmatMyDataClus, ikeep, jkeep) {

  MySize <- sqrt(length(rmatMyDataClus))

  # Generates cluster matrix (rmatMyDataClus); scan columns
  for (e in ikeep:(MySize - jkeep)) {

    # Generates cluster matrix (rmatMyDataClus); scan lines
    for (e2 in ikeep:(MySize - jkeep)) {
      rmatMyDataClus[e2, e] <- 1
    }
  }
  rmatMyDataClus[ikeep, (MySize - jkeep)] <- 0.5
  return(rmatMyDataClus)
}


###########################################################################
# Cluster detection (5) - initial cluster detection function:
###########################################################################

FirstPassClusters <- function(
  rmatMyDatafilt,
  rmatMyDataClus,
  Gap.penalty,
  Edge.Fill.rate,
  Surface.Fill.rate,
  Min.Clus.size,
  report_name,
  Report.type,
  corr.threshold
) {
  # Initialize variables
  # number of genes in contig = size of the correlation matrix
  MySize <- sqrt(length(rmatMyDatafilt))
  ikeep <- 0 # store line position when cluster found
  jkeep <- 0 # store column position when cluster found
  surfacesum <- 0 # store number of cells above corelation threshold in a putative cluster
  Cluster.Found <- FALSE # boolean storing position in/out of cluster
  i <- 0
  j <- 0
  surfacemax <- 0

  # counter scanning lines
  for (i in 2:(MySize)) {
    # counter scanning columns
    for (j in 1:(MySize - (i - 1))) {
      SkipCell <- SkipMatrixCell(
        corr.threshold,
        rmatMyDatafilt,
        i,
        j,
        ikeep,
        jkeep,
        Min.Clus.size,
        report_name,
        Report.type
      )

      if (SkipCell == FALSE) {
        EdgesAndGaps <- CalculateEdgesAndGaps(rmatMyDatafilt, i, j)
        if (Report.type == "verbose") {
          cat(
            paste(
              i , j, "edge length", (1 + MySize - i - j),
              " minimal edge support : ", EdgesAndGaps[1],
              ">?", (Edge.Fill.rate*(1+MySize-i-j)),
              " - maximal diag gap : ", EdgesAndGaps[2],
              "<?", (Gap.penalty * Edge.Fill.rate * (1 + MySize - i - j))),
            file = report_name,
            sep = "\n",
            append = TRUE
          )
        }

        # [*2] Edges and gaps checking
        if (
        (EdgesAndGaps[1] > (Edge.Fill.rate * (MySize - i - j))) &
        (EdgesAndGaps[2] < (Gap.penalty * (1 - Edge.Fill.rate) * (MySize - i - j)))
        ) {
          surfacesum <- CalculateSurfaceSum(rmatMyDatafilt, i, j)
          surfacemax <- 0

          # calculates number of cells in cluster triangle # rev 2021_07_30
          for (e in 1:(1 + MySize - i - j)) {
            surfacemax <- surfacemax + e
          }

          if (Report.type == "verbose") {
            cat(
              paste(
                "checking surface ", i, ", ", j,
                ": ", surfacesum,
                "/", (Surface.Fill.rate * surfacemax)
              ),
              file = report_name,
              sep = "\n",
              append = TRUE
            )
          }
          if (surfacesum > (Surface.Fill.rate * surfacemax)) {
            ikeep <- i
            jkeep <- j
            rmatMyDataClus <- GenerateClusterMatrix(
              rmatMyDataClus,
              ikeep,
              jkeep
            )
            if (Report.type == "verbose") {
              cat(
                paste(
                  "Edges & surface ok, cluster summit: ", ikeep,
                  " , ", jkeep,
                  " - ", (2 + MySize - jkeep - ikeep), " genes\n"
                ),
                file = report_name,
                sep = "\n",
                append = TRUE
              )
            }
          } else {
            next
          }
        } else {
          next
        }
      } else {
        next
      }
    }
  }
  return(rmatMyDataClus)
}


###########################################################################
# Cluster detection (6) - join clusters & clean up function:
###########################################################################

JoinClusters <- function(rmatMyDataClus, rmatMyDataA, corr.threshold) {
  MySize <- sqrt(length(rmatMyDataA))
  rmatMyDataClus <- apply(t(rmatMyDataClus), 2, rev) # rotate 90? counter-clockwise [note: switch the order of t() and apply()]
  rmatMyDataClus <- apply(t(rmatMyDataClus), 2, rev)
  i <- 0
  j <- 0

  #Cells of the diagonal within a cluster turned to 1's = join overlaping and contiguous clusters
  for (i in 1:(MySize - 1)) {
    if (
      (rmatMyDataClus[i + 1, i] == 1) &&
      (rmatMyDataA[i, i + 1] > corr.threshold)
    ) {
      rmatMyDataClus[i, i] <- 1
    } else {
      next
    }
  }

  # Cells of the diagonal at the edge of cluster with correlation below threshold turned to 0's
  stop_cleaning <- FALSE
  count_rounds <- 0
  repeat {
    count_rounds <- count_rounds + 1
    cluster_cleaned <- FALSE
    for (i in 1:(MySize - 1)) {
      if (
        (rmatMyDataClus[i, i] == 1) &&
        (rmatMyDataClus[i + 1, i + 1] == 0) &&
        (rmatMyDataA[i, i + 1] < corr.threshold)
      ) {
        rmatMyDataClus[i, i] <- 0
        cluster_cleaned <- TRUE
      } else {
        next
      }
    }
    for (j in 2:MySize) {
      if (
        (rmatMyDataClus[j, j] == 1) &&
        (rmatMyDataClus[j - 1, j - 1] == 0) &&
        (rmatMyDataA[j - 1, j] < corr.threshold)
      ) {
        rmatMyDataClus[j, j] <- 0
        cluster_cleaned <- TRUE
      } else {
        next
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
  cat(count_rounds, " rounds ")
  return(rmatMyDataClus)
}

###########################################################################
# Cluster detection (7) - finalize cluster rendering function:
###########################################################################

FinalizeClusterMatrix <- function(rmatMyDataClus, start_i, end_i) {
  # MySize <- nrow(rmatMyDataClus)
  for (e in start_i:end_i) { #Generates cluster matrix (rmatMyDataClus); scan columns
    for (f in start_i:end_i) { #Generates cluster matrix (rmatMyDataClus); scan lines
      rmatMyDataClus[f, e] <- 1
    }
  }
  return(rmatMyDataClus)
}

###########################################################################
# Cluster detection (8) - Clean up and scan final clusters function:
###########################################################################

SecondPassClusters<-function(rmatMyDataClus, rmatMyData, report_name, Min.Clus.size, Report.type) {
  Cluster.Found <- FALSE # boolean storing position in/out of cluster
  cluster_counter <- 0 # counter for the total number of clusters in contig
  clustered_genes <- 0 # counter for the total number of genes in cluster in contig
  clustered_all <- 0
  start_i <- 0
  end_i <- 0
  cluster_start <- NULL
  cluster_end <- NULL
  i <- 0
  j <- 0
  MySize <- sqrt(length(rmatMyDataClus))

  for (i in 1:MySize) {
    if (
      (rmatMyDataClus[i,i] == 1) &&
      (Cluster.Found == FALSE)
     ) { # first gene of a cluster
      Cluster.Found <- TRUE
      #clustered_genes<-1
      cluster_start <- rownames(rmatMyData)[i - 1]
      start_i <- i - 1
    } else if (
      (rmatMyDataClus[i, i] == 0) &&
      (Cluster.Found == TRUE) &&
      ((i + 1 - start_i) > Min.Clus.size)
    ) { # reached end of cluster
      Cluster.Found <- FALSE
      cluster_counter <- cluster_counter + 1
      clustered_genes <- (i - start_i)
      cluster_end <- rownames(rmatMyData)[i - 1]
      clustered_all <- clustered_all + clustered_genes
      end_i <- i - 1
      rmatMyDataClus <- FinalizeClusterMatrix(rmatMyDataClus, start_i, end_i)
      if (
        (Report.type == "normal") ||
        (Report.type == "verbose") ||
        (Report.type == "mini")
      ) {
        cat(
          paste(
            "Cluster ", cluster_counter, " :\t", cluster_start, "\t",
            cluster_end, "\t", clustered_genes, "\tgenes in cluster"
          ),
          file = report_name,
          sep = "\n",
          append = TRUE
        )
      }
      cat(
        "\nCluster ", cluster_counter,
        " detected: ", "start ", cluster_start, " - end ", cluster_end,
        " - ", clustered_genes, " genes in cluster"
      )

      j <- i
      repeat {
        rmatMyDataClus[j, i] <- 0
        j <- j + 1
        if (j == MySize + 1) {
          break
        }
      }
      k <- i
      repeat {
        rmatMyDataClus[i, k] <- 0
        k <- k - 1
        if (k < 2) {
          break
        }
      }
    } else if (
      (i == MySize) &&
      (Cluster.Found == TRUE) &&
      ((i + 1 - start_i) > Min.Clus.size)
    ) { # reched end of contig within a cluster
      Cluster.Found <- FALSE
      cluster_counter <- cluster_counter + 1
      clustered_genes <- (i - start_i)
      cluster_end <- rownames(rmatMyData)[i]
      clustered_all <- clustered_all + clustered_genes
      end_i <- i
      rmatMyDataClus <- FinalizeClusterMatrix(rmatMyDataClus, start_i, end_i)

      cat(
        "\nCluster ", cluster_counter, " detected: ",
        "start ", cluster_start,
        " - end ", cluster_end,
        " - ", clustered_genes, " genes in cluster"
      )

      if (
        (Report.type == "normal") ||
        (Report.type == "verbose") ||
        (Report.type == "mini")
      ) {
        cat(
          paste(
            "Cluster ", cluster_counter, " detected: ",
            "start ", cluster_start,
            " - end ", cluster_end,
            " - ", clustered_genes, " genes in cluster"
          ),
          file = report_name,
          sep = "\n",
          append = TRUE
        )
      }

      j <- i

      repeat {
        rmatMyDataClus[j, i] <- 0
        j <- j + 1

        if (j == MySize + 1) {
          break
        }
      }
      k <- i
      repeat {
        rmatMyDataClus[i, k] <- 0
        k <- k - 1
        if (k < 2) {
          break
        }
      }
    } else if (
      (rmatMyDataClus[i, i] == 0) &&
      (Cluster.Found == TRUE) &&
      ((i + 1 - start_i) < Min.Clus.size)
    ) {		# cluster too small = clean up matrix with 0's
      Cluster.Found <- FALSE

      for (e in start_i:i) {
        j <- e
        repeat {
          rmatMyDataClus[j, e] <- 0
          j <- j + 1
          if (j == MySize + 1) {
            break
          }
        }
        k <- e
        repeat {
          rmatMyDataClus[e, k] <- 0
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
        rmatMyDataClus[j, i] <- 0
        j <- j + 1
        if (j == MySize + 1) {
          break
        }
      }
      k <- i
      repeat {
        rmatMyDataClus[i, k] <- 0
        k <- k - 1
        if (k < 2) {
          break
        }
      }
    }
  }
  return(rmatMyDataClus)
}


###########################################################################
# Master function:
###########################################################################

RunClusterDetection <- function(
  MyData,
  MyGeneNames,
  AnalysisName,
  MyCut, # correlation threshold [0 to 1] from the averaged matrix (see [*1])
  UseQuantile, # boolean, whether to use quantile (TRUE) or corr coeff (FALSE) of averaged matrix as threshold
  Edge.Fill.rate, # % of cells above correlation threshold in the edge of a submatrix to proceed further in cluster detection [0 to 1]
  Surface.Fill.rate, # % of cells above correlation threshold in a submatrix call a cluster [0 to 0.5]
  Min.Clus.size, # Minimum number of genes in a cluster during first pass detection
  Gap.penalty, # Reduces the relative max gap length [positive] (see [*2])
  Report.type
) { # none, mini, genelist, normal or verbose
  StartAll.time <- Sys.time()
  Gap.penalty <- (1 / Gap.penalty)

  if ((Report.type == "list") || (Report.type == "mini")) {
    report_name <- paste(AnalysisName, "_GenomeClustR report.txt", sep = "")
  } else if ((Report.type == "normal") || (Report.type == "verbose")) {
    report_name <- paste0(
      AnalysisName,
      as.character(format(Sys.time(), "%Y%m%d%H%M%S")),
      "_GenomeClustR report.txt",
    )
    cat(
      paste(
        "GenomeClustR analysis ", AnalysisName,
        as.character(format(Sys.time(), "%Y_%m_%d %H:%M:%S"))
      ),
      file = report_name,
      sep = "\n",
      append = TRUE
    )
    cat(
      paste(
        "\nEdge.Fill.rate: ",
        Edge.Fill.rate,
        "\nSurface.Fill.rate: ",
        Surface.Fill.rate,
        "\nMin.Clus.size:",
        Min.Clus.size
      ),
      file = report_name,
      sep = "\n",
      append = TRUE
    )
    cat(
      paste(
        "\nGap.penalty: ", Gap.penalty,
        " - longest gap allowed ", (100 * (Gap.penalty * (1 - Edge.Fill.rate))),
        "% of cluster size"
      ),
      file = report_name,
      sep = "\n",
      append = TRUE
    )
  } else if (Report.type == "none" ) {
    cat("\nNo report file mode")
  }

  # Running cluster detection

  cat("\nRunning pairwise correlation calculation... ")
  StartStep.time <- Sys.time()
  rmatMyData <- RawCorMat(MyData, MyGeneNames)
  MySize <- sqrt(length(rmatMyData))
  StepDuration <- difftime(Sys.time(), StartStep.time, units="secs")
  if ((Report.type == "normal") || (Report.type=="verbose")) {
    cat(
      paste("1. Pairwise correlation calculation runtime: ", StepDuration),
      file = report_name,
      sep = "\n",
      append = TRUE
    )
  }
  cat("done, runtime ", StepDuration, " secs")

  if (Report.type == "verbose") {
    write.table(
      rmatMyData,
      file = paste(report_name, "_1_Raw_correlation_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )

    pdf(file = paste(report_name, "_1_Raw_correlation_matrix.pdf", sep = ""))
    heatmap.2(
      rmatMyData,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = TRUE
    )
    dev.off()
  }

  cat("\nRunning correlation matrix local averaging... ")
  StartStep.time <- Sys.time()
  rmatMyDataA <- AveCorMat(rmatMyData)
  StepDuration <- difftime(Sys.time(), StartStep.time, units = "secs")
  cat("done, runtime ", StepDuration, " secs")

  # preparing analysis report part1

  if (UseQuantile == TRUE) {
    cat(
      paste(
        "Gene clusters analysis ", AnalysisName,
        "cor. threshold:", quantile(rmatMyDataA, MyCut, na.rm = TRUE)
      ),
      file = report_name,
      sep = "\n",
      append = TRUE
    )
  } else {
    cat(paste(
      "Gene clusters analysis ",
      AnalysisName,
      "cor. threshold:",
      MyCut,
      file = report_name,
      sep = "\n",
      append = TRUE
    ))
  }

  if (Report.type == "verbose") {
    write.table(
      rmatMyDataA,
      file = paste(report_name, "_2_Averaged_correlation_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )
    pdf(file = paste(report_name, "_2_Averaged_correlation_matrix.pdf", sep = ""))
    heatmap.2(
      rmatMyDataA,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = TRUE
    )
    dev.off()
  }

  cat("\nPreparing temporary matrices... ")
  StartStep.time <- Sys.time()
  if (UseQuantile == TRUE) {
    corr.threshold <- quantile(rmatMyDataA, MyCut, na.rm = TRUE) 	# [*1] 2021_04_26 changed raw matrix quantile to averaged matrix quantile
  } else {
    corr.threshold <- MyCut
  }
  rmatMyDatafilt <- rmatMyDataA # filtered matrix receives 1 in cells above correlation threshold, 0 otherwise
  rmatMyDatafilt[rmatMyDatafilt >= corr.threshold] <- 1
  rmatMyDatafilt[rmatMyDatafilt < corr.threshold] <- 0
  rmatMyDatafilt[is.na(rmatMyDatafilt) == TRUE] <- 0

  if (Report.type == "verbose") {
    write.table(
      rmatMyDatafilt,
      file = paste(report_name, "_3_Raw_thresholded_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )

    pdf(file = paste(report_name, "_3_Raw_thresholded_matrix.pdf", sep = ""))
    heatmap.2(
      rmatMyDatafilt,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = FALSE
    )
    dev.off()
  }

  rmatMyDatafilt <- apply(t(rmatMyDatafilt),2, rev) # rotate matrix 90? counter clockwise to facilitate cluster detection
  rmatMyDataClus <- rmatMyDatafilt # Matrix in which the position of detected clusters will be stored (in the form of '1's)
  rmatMyDataClus[rmatMyDatafilt == 1] <- 0 # initialize all positions to 0's
  StepDuration <- difftime(Sys.time(), StartStep.time, units="secs")
  cat("done, runtime ", StepDuration, " secs")

  # preparing analysis report part2

  if ((Report.type == "normal") || (Report.type == "verbose")) {
    if (UseQuantile == TRUE) {
      cat(
        paste(
          "correlation threshold: ", corr.threshold,
          "(", MyCut, "th quantile)"
        ),
        file = report_name,
        sep = "\n",
        append = TRUE
      )
    } else {
      cat(
        paste(
          "correlation threshold: ", corr.threshold,
          " (quantile thresholding disabled)"
        ),
        file = report_name,
        sep = "\n",
        append = TRUE
      )
    }
  }
  cat("###############", file = report_name, sep = "\n", append = TRUE)

  cat("\nRunning first pass cluster detection... ")
  StartStep.time <- Sys.time()
  rmatMyDataClus <- FirstPassClusters(
    rmatMyDatafilt,
    rmatMyDataClus,
    Gap.penalty,
    Edge.Fill.rate,
    Surface.Fill.rate,
    Min.Clus.size,
    report_name,
    Report.type,
    corr.threshold
  )

  if (Report.type == "verbose") {
    write.table(
      rmatMyDataClus,
      file = paste(report_name, "_4_Rotated_1stPass_thresholded_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )
    pdf(file = paste0(report_name, "_4_Rotated_1stPass_thresholded_matrix.pdf"))
    heatmap.2(
      rmatMyDataClus,
      Rowv = NA, Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = FALSE
    )
    dev.off()

    rmatTemp <- t(apply(rmatMyDataClus, 2, rev))
    rmatTemp <- t(apply(rmatTemp, 2, rev))
    rmatTemp <- t(rmatTemp)

    write.table(
      rmatTemp,
      file = paste(report_name, "_5_Native_1stPass_thresholded_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )
    pdf(file = paste0(report_name, "_5_Native_1stPass_thresholded_matrix.pdf"))
    heatmap.2(
      rmatTemp,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = FALSE
    )
    dev.off()
  }
  StepDuration <- difftime(Sys.time(), StartStep.time, units = "secs")
  cat("done, runtime ", StepDuration, " secs")

  cat("\nRunning cluster refinement... ")
  StartStep.time <- Sys.time()
  rmatMyDataClus <- JoinClusters(rmatMyDataClus, rmatMyDataA, corr.threshold)

  if (Report.type == "verbose") {
    write.table(
      rmatMyDataClus,
      file = paste(report_name, "_6_Refined_thresholded_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )

    pdf(file = paste0(report_name, "_6_Refined_thresholded_matrix.pdf"))
    heatmap.2(
      rmatMyDataClus,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = FALSE
    )
    dev.off()
  }

  StepDuration <- difftime(Sys.time(), StartStep.time, units = "secs")
  cat("done, runtime ", StepDuration, " secs")

  cat("\nRunning second pass cluster detection...\n")
  StartStep.time <- Sys.time()
  rmatMyDataClus <- SecondPassClusters(
    rmatMyDataClus, rmatMyData,
    report_name, Min.Clus.size, Report.type
  )
  clustered_all <- 0

  if (Report.type == "verbose") {
    cat("All clustered genes: ", file = report_name, sep = "\n", append = TRUE)
  }
  for (i in 1:MySize) {
    if (rmatMyDataClus[i, i] == 1) {
      clustered_all <- clustered_all + 1
      if (
        (Report.type == "verbose") ||
        (Report.type == "list") ||
        (Report.type == "normal")
      ) {
        cat(
          paste(rownames(rmatMyData)[i]),
          file = report_name,
          sep = "\n",
          append = TRUE
        )
      }
    } else {
      next
    }
  }
  cluster_counter <- 0
  for (i in 1:(MySize - 1)) {
    if (
      (rmatMyDataClus[i, i] == 0) &&
      (rmatMyDataClus[i + 1, i + 1] == 1)
    ) {
      cluster_counter <- cluster_counter + 1
    }
    if (rmatMyDataClus[1, 1] == 1) {
      cluster_counter <- cluster_counter + 1
    }
  }

  if (Report.type == "verbose") {
    write.table(
      rmatMyDataClus,
      file = paste(report_name, "_7_Final_thresholded_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )

    pdf(file = paste0(report_name, "_7_Final_thresholded_matrix.pdf"))
    heatmap.2(
      rmatMyDataClus,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE, cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none",
      key = FALSE
    )
    dev.off()
  }

  StepDuration <- difftime(Sys.time(), StartStep.time, units = "secs")
  cat("\nsecond pass complete, runtime ", StepDuration, " secs")
  cat(
    "\nTotal ", cluster_counter,
    "clusters found encompassing", clustered_all,
    "/", MySize, "genes\n"
  )

  # Combine averaged matrix with detected cluster matrix and save as files
  rmatMyDatafinal <- t(rmatMyDataA) + 0.3 * rmatMyDataClus
  for (i in 1:MySize) { # Mask averaged correlations outside clusters
    for (j in 1:MySize) {
      if (rmatMyDataClus[i, j] == 0) {
        rmatMyDatafinal[i, j] <- 0
      }
    }
  }
  for (i in 1:MySize) { # Display clusters along the first diagonal
    rmatMyDatafinal[i, i] <- rmatMyDataClus[i, i]
  }

  rmatMyDatafinal[upper.tri(rmatMyDatafinal)] <- rmatMyData[upper.tri(rmatMyDataClus)] # Keeping only upper triangle

  if (Report.type == "verbose") {
    write.table(
      rmatMyDatafinal,
      file = paste(report_name, "_8_Final_summary_matrix.tab"),
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t"
    )
  }

  if ((Report.type == "normal") || (Report.type == "verbose")) {
    pdf(file = paste0(report_name, ".pdf"))
    heatmap.2(
      rmatMyDatafinal,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none"
    )
    dev.off()

    jpeg(
      file = paste0(report_name, "_clusters.jpg"),
      width = 4800,
      height = 4800
    )

    heatmap.2(
      rmatMyDatafinal,
      Rowv = NA,
      Colv = NA,
      labRow = rownames(rmatMyData),
      labCol = colnames(rmatMyData),
      col = mypalette(100),
      na.rm = TRUE,
      cexRow = 1,
      cexCol = 1,
      dendrogram = "none",
      trace = "none"
    )
    dev.off()

    cat("pdf and jpg image files saved, txt report saved\n")
  }

  StepDuration <- difftime(Sys.time(), StartAll.time, units = "secs")

  #Finalizing report
  cat("###############", file = report_name, sep = "\n", append = TRUE)
  cat("Total runtime: ", StepDuration, " secs\n")
}
