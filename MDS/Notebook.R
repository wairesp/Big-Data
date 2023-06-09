get_procrustes_parameters <- function(x, target, translation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (translation) {
    c.target <- scale(target, center = TRUE, scale = FALSE)
    m.target <- attr(c.target, "scaled:center")

    c.x <- scale(x, center = TRUE, scale = FALSE)
    m.x <- attr(c.x, "scaled:center")

    matrix_prod <- t(c.target) %*% x
    svd_results <- svd(matrix_prod)
    rotation_matrix <- svd_results$v %*% t(svd_results$u)

    translation_vector <- m.target - t(rotation_matrix) %*% m.x

  } else {
    matrix_prod <- t(target) %*% x
    svd_results <- svd(matrix_prod)
    rotation_matrix <- svd_results$v %*% t(svd_results$u)
    translation_vector <- matrix(data = 0, nrow = n_col, ncol = 1)
  }

  return(list(rotation_matrix = rotation_matrix, translation_vector = translation_vector))
}

perform_procrustes <- function(x, target, matrix_to_transform, translation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (n_row != nrow(target)) {
    stop("\"x\" and \"target\" do not have the same number of rows")
  }

  if (n_col != ncol(target)) {
    stop("\"x\" and \"target\" do not have the same number of columns")
  }

  if (n_col != ncol(matrix_to_transform)) {
    stop("\"x\" and \"matrix_to_transform\" do not have the same number of columns")
  }

  procrustes_parameters <- get_procrustes_parameters(x = x, target = target, translation = translation)

  ones_vector <- matrix(data = 1, nrow = nrow(matrix_to_transform), ncol = 1)
  translation_matrix <- matrix(data = 0, nrow = nrow(matrix_to_transform), ncol = ncol(matrix_to_transform))
  translation_matrix <- translation_matrix + ones_vector %*% t(procrustes_parameters$translation_vector)

  return(matrix_to_transform %*% procrustes_parameters$rotation_matrix + translation_matrix)
}

classical_mds <- function(x, k, dist_fn, return_distance_matrix = FALSE, ...) {

  if (!is.function(dist_fn)) {
    stop("\"dist_fn\" must be a function")
  } else if (any(is.na(x))) {
    stop("There are some \"NA\" values in the data. Please remove them")
  }

  mds <- list()
  dist_matrix <- dist_fn(x, ...)
  mds_result <- stats::cmdscale(d = dist_matrix, k = k, eig = TRUE)

  mds$points <- mds_result$points
  mds$eigen <- mds_result$eig[1:k]
  mds$GOF <- mds_result$GOF

  if (return_distance_matrix) {
    mds$distance <- as.matrix(dist_matrix)
  }

  return(mds)
}

get_partitions_for_divide_conquer <- function(n, l, c_points, r) {

  if (l-c_points <= 0) {
    stop("\"l\" must be greater than \"c_points\"")
  } else if (l-c_points <= c_points) {
    stop("\"l-c_points\" must be greater than \"c_points\"")
  } else if (l-c_points <= r) {
    stop ("\"l-c_points\" must be greater than \"r\"")
  }

  permutation <- sample(x = n, size = n, replace = FALSE)

  if (n<=l) {
    list_indexes <- list(permutation)
  } else {
    min_sample_size <- max(r+2, c_points)
    p <- 1 + ceiling((n-l)/(l-c_points))
    last_partition_sample_size <- n - (l + (p-2) * (l-c_points))

    if (last_partition_sample_size < min_sample_size & last_partition_sample_size > 0) {
      p <- p - 1
      last_partition_sample_size <- n - (l + (p-2) * (l-c_points))
    }

    first_parition <- permutation[1:l]
    last_partition <- permutation[(n-last_partition_sample_size+1):n]
    list_indexes <- split(x = permutation[(l+1):(n-last_partition_sample_size)], f = 1:(p-2))
    names(list_indexes) <- NULL
    list_indexes[[p-1]] <- list_indexes[[1]]
    list_indexes[[p]] <- last_partition
    list_indexes[[1]] <- first_parition
  }

  return(list_indexes)
}

split_matrix <- function(matrix, num_points) {
  x_1 <- matrix[1:num_points, , drop = FALSE]
  x_2 <- matrix[(num_points + 1):nrow(matrix), , drop = FALSE]
  return(list(x_1, x_2))
}

main_divide_conquer_mds <- function(idx, x, x_sample_1, r, original_mds_sample_1, dist_fn, ...) {
  # Filter the matrix
  x_filtered <- x[idx, , drop = FALSE]

  # Join with the sample from the first partition
  x_join_sample_1 <- rbind(x_sample_1, x_filtered)

  # Perform MDS
  mds_all <- classical_mds(x = x_join_sample_1, k = r, dist_fn = dist_fn, ...)
  mds_points <- mds_all$points
  mds_eigen <- mds_all$eigen
  mds_GOF <- mds_all$GOF

  # Split MDS into two parts: the first part is the MDS for the sampling points of the
  # first partition and the second part is the MDS for the points of the partition
  mds_split <- split_matrix(matrix = mds_points, num_points = nrow(x_sample_1))
  mds_sample_1 <- mds_split[[1]]
  mds_partition <- mds_split[[2]]
  mds_procrustes <- perform_procrustes(x = mds_sample_1,
                                       target = original_mds_sample_1,
                                       matrix_to_transform = mds_partition,
                                       translation = FALSE)

  mds_eigen <- mds_eigen/length(idx)
  mds_GOF <- mds_GOF*length(idx)

  return(list(points = mds_procrustes, eigen = mds_eigen, GOF = mds_GOF))
}

divide_conquer_mds <- function(x, l, c_points, r, n_cores = 1, dist_fn = stats::dist, ...) {

  n_row_x <- nrow(x)

  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, k = r, dist_fn = dist_fn, ...)
    mds_to_return$eigen <- mds_to_return$eigen/n_row_x
    mds_to_return$GOF <- mds_to_return$GOF

  } else {

    mds_matrix <- matrix(data = NA, nrow = n_row_x, ncol = r)

    # Generate indexes list. Each element corresponds to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, c_points = c_points, r = r)
    num_partitions <- length(idx)
    length_1 <- length(idx[[1]])

    # Perform MDS for the first partition
    x_1 <- x[idx[[1]], , drop = FALSE]
    mds_1 <- classical_mds(x = x_1, k = r, dist_fn = dist_fn, ...)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/length_1
    mds_1_GOF <- mds_1$GOF*length_1

    # Take a sample from the first partition
    sample_1 <- sample(x = length_1, size = c_points, replace = FALSE)
    x_sample_1 <- x_1[sample_1, ,drop = FALSE]
    mds_sample_1 <- mds_1_points[sample_1, , drop = FALSE]

    mds_others_results <- parallel::mclapply(idx[2:num_partitions],
                                             main_divide_conquer_mds,
                                             x = x,
                                             x_sample_1 = x_sample_1,
                                             r = r,
                                             original_mds_sample_1 = mds_sample_1,
                                             dist_fn = dist_fn,
                                             mc.cores = n_cores,
                                             ...)

    # Obtain points
    mds_others_points <- do.call(rbind, parallel::mclapply(mds_others_results, function(x) x$points, mc.cores = n_cores))
    mds_matrix[1:length_1, ] <- mds_1_points
    mds_matrix[(length_1 + 1):n_row_x, ] <- mds_others_points
    order_idx <- do.call(c, idx)
    order_idx <- order(order_idx)
    mds_matrix <- mds_matrix[order_idx, , drop = FALSE]
    mds_matrix <- apply(mds_matrix, MARGIN = 2, FUN = function(x) x - mean(x))
    mds_matrix <- mds_matrix %*% eigen(stats::cov(mds_matrix))$vectors

    # Obtain eigenvalues
    eigen <- parallel::mclapply(mds_others_results, function(x) x$eigen, mc.cores = n_cores)
    eigen[[num_partitions]] <- mds_1_eigen
    eigen <- Reduce(`+`, eigen)
    eigen <- eigen/num_partitions

    # Obtain GOF
    GOF <- parallel::mclapply(mds_others_results, function(x) x$GOF, mc.cores = n_cores)
    GOF[[num_partitions]] <- mds_1_GOF
    GOF <- Reduce(`+`, GOF)
    GOF <- GOF/n_row_x
    mds_to_return <- list(points = mds_matrix, eigen = eigen, GOF = GOF)
  }

  return(mds_to_return)
}

get_partitions_for_interpolation <- function(n, n_obs, l, r) {

  if (l<=r) {
    stop("\"l\" must be greater than \"r\"")
  }

  if (n<=l) {
    p <- 1
  } else{
    p <- 1 + ceiling((n - l)/n_obs)
    n_last <- n - (l + (p-2)*n_obs)
  }

  permutation <- sample(x = n, size = n, replace = FALSE)

  if (p>1) {
    first_part <- permutation[1:l]
    middle_part <- permutation[(l+1):(n-n_last)]
    last_part <- permutation[(n-n_last+1):n]

    list_index <- split(middle_part, 1:(p-2))
    names(list_index) <- NULL
    list_index[[p-1]] <- list_index[[1]]
    list_index[[1]] <- first_part
    list_index[[p]] <- last_part

  } else {
    list_index <- list(permutation)
  }

  return(list_index)
}

get_P_matrix <- function(n_row) {
  identity_matrix <- diag(x = 1, nrow = n_row, ncol = n_row)
  one_vector <- matrix(data = 1, nrow = n_row, ncol = 1)
  P <- identity_matrix - 1/n_row * one_vector %*% t(one_vector)
  return(P)
}

interpolation_mds_main <- function(idx, x, data_1, x_1, n_row_1, q_vector, x_1__s_1__inv, dist_fn, ...) {

  # Filter the matrix
  x_other <- x[idx, , drop = FALSE]
  n_row_other <- nrow(x_other)
  n_row_1 <- nrow(data_1)

  # Get A matrix
  full_A <- dist_fn(x = rbind(x_other, data_1), ...)
  full_A <- as.matrix(full_A)
  n_col_A <- ncol(full_A)
  A <- full_A[1:n_row_other, (n_col_A-n_row_1+1):n_col_A, drop = FALSE]

  # Get delta matrix
  A_sq <- A^2

  # One vecotr
  one_vector_other <- matrix(data = 1, nrow = n_row_other, ncol = 1)

  # Get coordinates
  x_2 <- 1/(2*n_row_1) * (one_vector_other %*% t(q_vector) - A_sq) %*% x_1__s_1__inv

  return(x_2)
}

interpolation_mds <- function(x, l, r, n_cores = 1, dist_fn = stats::dist,...) {

  n <- nrow(x)
  n_row_partition <- l
  indexes_partition <- get_partitions_for_interpolation(n = n, n_obs = n_row_partition, l = l, r = r)
  num_partitions <- length(indexes_partition)

  if (num_partitions <= 1) {

    # It is possible to run MDS directly
    mds <- classical_mds(x = x, k = r, dist_fn = dist_fn, ...)
    points <- mds$points
    eigen_v <- mds$eigen/n
    GOF <- mds$GOF
    list_to_return <- list(points = points, eigen = eigen_v, GOF = GOF)

  } else {

    # Get the first group
    n_row_1 <- length(indexes_partition[[1]])

    # Obtain MDS for the first group
    data_1 <- x[indexes_partition[[1]], ,drop = FALSE]
    mds_eig <- classical_mds(x = data_1, k = r, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    distance_matrix <- as.matrix(mds_eig$distance)
    X_1 <- mds_eig$points
    eigen_v <- mds_eig$eigen/nrow(X_1)
    GOF <- mds_eig$GOF

    # Get P matrix
    P <- get_P_matrix(n_row = n_row_1)
    Q <- -1/2 * P %*% distance_matrix^2 %*% t(P)
    q_vector <- diag(Q)
    S <- 1 / (n_row_1-1) * t(X_1) %*% X_1
    x_1__s_1__inv <- X_1 %*% solve(S)

    # Calculations needed to do Gower interpolation
    mds_others <- parallel::mclapply(indexes_partition[2:num_partitions],
                                     interpolation_mds_main,
                                     x = x,
                                     data_1 = data_1,
                                     x_1 = X_1,
                                     n_row_1 = n_row_1,
                                     q_vector = q_vector,
                                     x_1__s_1__inv = x_1__s_1__inv,
                                     dist_fn = dist_fn,
                                     mc.cores = n_cores,
                                     ...)

    mds_points <- matrix(data = NA, nrow = n, ncol = r)
    mds_points[1:n_row_1, ] <- X_1
    mds_points[(n_row_1+1):n, ] <- do.call(rbind, mds_others)
    idx_all <- do.call(c, indexes_partition)
    idx_all <- order(idx_all)
    mds_points <- mds_points[idx_all, , drop = FALSE]
    mds_points <- apply(mds_points, MARGIN = 2, FUN = function(y) y - mean(y))
    mds_points <- mds_points %*% eigen(stats::cov(mds_points))$vectors
    list_to_return <- list(points = mds_points, eigen = eigen_v, GOF = GOF)
  }

  return(list_to_return)
}

get_partitions_for_fast <- function(n, l, s_points, r) {

  if (n < s_points) {
    stop("Number of rows of \"x\" must be greater than \"s_points\"")
  } else if (n*s_points/l < r) {
    stop("Number of rows of \"x\" must be greater than \"r\" x \"l\"/\"s_points\"")
  }

  p <- floor(l/s_points)
  min_sample_size <- max(r+2, s_points)
  size_partition <- floor(n/p)
  last_sample_size <- n - (p-1) * size_partition

  if (last_sample_size < min_sample_size & last_sample_size > 0) {
    p <- p - 1
  }

  permutation <- sample(x = n, size = n, replace = FALSE)
  permutation_all <- permutation[1:((p-1)*size_partition)]
  permutation_last <- permutation[((p-1)*size_partition+1):n]
  list_indexes <- split(x = permutation_all, f = 1:(p-1))
  names(list_indexes) <- NULL
  list_indexes[[p]] <- permutation_last

  return(list_indexes)
}

main_fast_mds <- function(idx, matrix, l, s_points, r, n_cores, dist_fn, ...) {

  # Partition the matrix
  x_partition <- matrix[idx, , drop = FALSE]

  # Apply the method
  mds <- fast_mds(x = x_partition, l = l, s_points = s_points, r = r, n_cores = n_cores, dist_fn = dist_fn, ...)

  return(mds)
}


fast_mds <- function(x, l, s_points, r, n_cores = 1, dist_fn = stats::dist, ...) {
  n <- nrow(x)

  if (n <= l) {

    mds <- classical_mds(x = x, k = r, dist_fn = dist_fn, ...)
    mds$eigen <- mds$eigen / nrow(x)

    return(mds)

  } else {

    # Split x
    index_partition <- get_partitions_for_fast(n = n, l = l, s_points = s_points, r = r)
    num_partition <- length(index_partition)

    # Apply MDS to all the partitions
    mds_partition <- parallel::mclapply(index_partition,
                                        main_fast_mds,
                                        matrix = x,
                                        l = l,
                                        s_points = s_points,
                                        r = r,
                                        n_cores = n_cores,
                                        dist_fn = dist_fn,
                                        mc.cores = n_cores,
                                        ...)

    mds_partition_points <- parallel::mclapply(mds_partition, function(x) x$points, mc.cores = n_cores)
    mds_partition_eigen <- parallel::mclapply(mds_partition, function(x) x$eigen, mc.cores = n_cores)
    mds_GOF <- parallel::mclapply(mds_partition, function(x) x$GOF, mc.cores = n_cores)

    # take a sample for each partition
    length_partition <- parallel::mclapply(index_partition, length, mc.cores = n_cores)
    sample_partition <- parallel::mclapply(length_partition, sample, size = s_points, replace = FALSE, mc.cores = n_cores)
    indexes_filtered <- parallel::mcmapply(function(idx, sample) idx[sample],
                                           idx = index_partition,
                                           sample = sample_partition,
                                           SIMPLIFY = FALSE,
                                           mc.cores = n_cores)

    length_sample <- parallel::mclapply(sample_partition, length, mc.cores = n_cores)

    indexes_scaled <- parallel::mcmapply(function(i, long) ((i-1)*long + 1):(i*long),
                                         i = 1:num_partition,
                                         long = length_sample,
                                         SIMPLIFY = FALSE,
                                         mc.cores = n_cores)

    # Join all the points
    x_partition_sample <- parallel::mclapply(indexes_filtered,
                                             function(index_partitions, matrix) { matrix[index_partitions, , drop = FALSE] },
                                             matrix = x,
                                             mc.cores = n_cores)

    x_M <- do.call(rbind, x_partition_sample)

    # Apply MDS to the subsampling points
    mds_M <- classical_mds(x = x_M, k = r, dist_fn = dist_fn, ...)
    mds_M_points <- mds_M$points

    # Extract the MDS configuration for the sampling points from mds_M_points
    mds_M_sampling_points <- parallel::mclapply(indexes_scaled,
                                                function(indexes_scaled, matrix) { matrix[indexes_scaled, , drop = FALSE] },
                                                matrix = mds_M_points,
                                                mc.cores = n_cores)

    # Extract the MDS configuration for the sampling points from mds_partition_points
    mds_partition_sampling_points <- parallel::mcmapply(function(matrix, index_partitions, idx) { matrix[idx, , drop = FALSE] },
                                                        matrix = mds_partition_points,
                                                        idx = sample_partition,
                                                        SIMPLIFY = FALSE,
                                                        mc.cores = n_cores)

    # Apply Procrustes
    procrustes <- parallel::mcmapply(perform_procrustes,
                                     x = mds_partition_sampling_points,
                                     target = mds_M_sampling_points,
                                     matrix_to_transform = mds_partition_points,
                                     translation = FALSE,
                                     SIMPLIFY = FALSE,
                                     mc.cores = n_cores)

    # Build the list to be returned
    idx_order <- Reduce(c, index_partition)
    idx_order <- order(idx_order)
    mds <-do.call(rbind, procrustes)
    mds <- mds[idx_order, ,drop = FALSE]
    mds <- apply(mds, MARGIN = 2, FUN = function(y) y - mean(y))
    mds <- mds %*% eigen(stats::cov(mds))$vectors
    eigen <- Reduce(`+`, mds_partition_eigen)/num_partition

    # Build GOF metric
    GOF <- parallel::mcmapply(function(x, y) x*length(y),
                              x = mds_GOF,
                              y = index_partition,
                              SIMPLIFY = FALSE,
                              mc.cores = n_cores)
    GOF <- Reduce(`+`, GOF)/n

    return(list(points = mds, eigen = eigen, GOF = GOF))
  }
}

#
# function fast_mds(x, l, s_points, r, n_cores, dist_fn, ...)
#   n = número de filas en x
#   if n <= l
#     realizar MDS clásico en x
#     return resultado de MDS
#   else
#     particionar x en subconjuntos más pequeños usando get_partitions_for_fast
#     for cada partición
#       aplicar fast_mds recursivamente
#       store resultado de MDS para la partición
#     end for
#     tomar una muestra de puntos de cada partición
#     realizar MDS clásico en los puntos de muestra para obtener configuración de alineación
#     for cada partición
#       alinear configuración de MDS con configuración de alineación usando perform_procrustes
#       store resultado de MDS alineado para la partición
#     end for
#       combinar resultados de MDS alineados de todas las particiones
#       return resultado de MDS alineado combinado
#   end if
# end function
