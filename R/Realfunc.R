
#all the commenting/documentation will be added later
#real funs of the package 
generate_partitions <- function(n) {
  partitions <- list(list())
  for (i in 1:n) {
    new_partitions <- list()
    for (p in partitions) {
      for (j in 1:length(p)) {
        if (i >= p[[j]][length(p[[j]])]) {
          new_partitions <- c(new_partitions, list(c(p[[j]], i)))
        }
      }
      new_partitions <- c(new_partitions, list(list(i)))
    }
    partitions <- new_partitions
  }
  return(partitions)
}

is_standard_tableau <- function(tableau) {
  rows <- length(tableau)
  cols <- max(sapply(tableau, length))
  for (i in 1:cols) {
    col_values <- sapply(tableau, function(row) sum(row == i))
    if (!identical(col_values, sort(col_values))) {
      return(FALSE)
    }
  }
  for (i in 1:rows) {
    row_values <- sapply(1:cols, function(j) sum(tableau[[i]] == j))
    if (!identical(row_values, sort(row_values))) {
      return(FALSE)
    }
  }
  return(TRUE)
}

generate_standard_tableaux <- function(partition) {
  tableaux <- list()
  stack <- list(list(partition, list()))
  while (length(stack) > 0) {
    top <- stack[[length(stack)]]
    remaining <- top[[1]]
    tableau <- top[[2]]
    stack <- stack[-length(stack)]
    if (length(remaining) == 0) {
      if (is_standard_tableau(tableau)) {
        tableaux <- c(tableaux, list(tableau))
      }
    } else {
      row <- ifelse(length(tableau) == 0, 0, tail(tableau, 1))
      for (i in 1:remaining[[1]]) {
        stack <- c(stack, list(list(remaining[-1], c(tableau, row + i))))
      }
    }
  }
  return(tableaux)
}

is_lattice_permutation <- function(word) {
  counts <- integer(max(word))
  for (i in 1:length(word)) {
    counts[word[i]] <- counts[word[i]] + 1
    for (j in 1:(word[i] - 1)) {
      if (counts[j] < counts[j + 1]) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

generate_lattice_permutations <- function(digits) {
  lattice_perms <- list()

  generate_permutations <- function(digits, current_permutation) {
    if (length(current_permutation) == length(digits)) {
      if (is_lattice_permutation(current_permutation)) {
        lattice_perms <<- c(lattice_perms, list(current_permutation))
      }
      return
    }
    for (digit in unique(digits)) {
      if (sum(current_permutation == digit) < sum(digits == digit)) {
        generate_permutations(digits, c(current_permutation, digit))
      }
    }
  }

  generate_permutations(digits, integer(0))

  return(lattice_perms)
}

order_vector <- function(input_vector) {
  ordered_vector <- sort(input_vector, decreasing = TRUE)
  unique_elements <- unique(ordered_vector)
  output_vector <- numeric(sum(unique_elements))
  start_index <- 1
  for (i in 1:length(unique_elements)) {
    count <- sum(ordered_vector == unique_elements[i])
    output_vector[start_index:(start_index + count - 1)] <- rep(unique_elements[i], count)
    start_index <- start_index + count
  }
  return(output_vector)
}

generate_partitions <- function(n) {
  partitions <- list()
  generate_partitions_recursive(n, integer(0), partitions)
  return(partitions)
}

generate_partitions_recursive <- function(n, partition, partitions) {
  if (n == 0) {
    partitions[[length(partitions) + 1]] <- partition
    return
  }
  for (k in 1:n) {
    new_partition <- c(partition, k)
    generate_partitions_recursive(n - k, new_partition, partitions)
  }
}

decompose <- function(n) {
  if (length(q) >= n && !is.null(q[[n]])) {
    return(q[[n]])
  }

  result <- list(list(n))

  for (i in 1:(n - 1)) {
    a <- n - i
    R <- decompose(i)
    for (r in R) {
      if (r[[1]] <= a) {
        result <- c(result, list(c(a, r)))
      }
    }
  }

  q[[n]] <- result
  return(unlist(result))
}

generate_permutations <- function(n) {
    permutations <- list()

    generate_permutations_recursive <- function(perm, used, n) {
        if (length(perm) == n) {
            permutations <<- c(permutations, list(perm))
            return
        }
        for (i in 1:n) {
            if (!i %in% used) {
                generate_permutations_recursive(c(perm, i), c(used, i), n)
            }
        }
    }

    generate_permutations_recursive(integer(0), integer(0), n)
    return(permutations)
}

generate_transpositions <- function(perm) {
  transpositions <- list()
  n <- length(perm)

  for (i in 1:(n - 1)) {
    if (perm[i] != i) {
      j <- which(perm == i)[1]
      transpositions <- c(transpositions, list(c(i, j)))
      perm[j] <- perm[i]
    }
  }

  return(transpositions)
}

generate_transpositions2 <- function(permutation) {
  transpositions <- list()
  n <- length(permutation)

  for (i in 1:(n - 1)) {
    transposition <- c(permutation[i], permutation[i + 1])
    transpositions <- c(transpositions, list(transposition))
  }

  return(transpositions)
}

permutation_to_cycle <- function(permutation) {
  visited <- rep(FALSE, length(permutation))
  cycle <- c()

  for (i in 1:length(permutation)) {
    if (!visited[i]) {
      j <- i
      while (!visited[j]) {
        visited[j] <- TRUE
        cycle <- c(cycle, j)
        j <- permutation[j]
      }
    }
  }

  return(cycle)
}

pairs_from_vector <- function(vector) {
  pairs <- list()
  n <- length(vector)

  for (i in 1:(n - 1)) {
    pair <- c(vector[i], vector[i + 1])
    pairs <- c(pairs, list(pair))
  }

  return(pairs)
}

convert_to_vector <- function(pair) {
  if (length(pair) != 2) {
    stop("Input vector must have length 2.")
  }

  start <- min(pair)
  end <- max(pair)
  vector <- c(start, (start + 1):(end - 1), end, (end - 1):start)
  return(vector)
}

calculate_special_distance <- function(matrix, i) {
  # Find position of i in the matrix
  position_i <- which(matrix == i, arr.ind = TRUE)

  if (nrow(position_i) == 0) {
    stop("Integer i not found in the matrix.")
  } else if (nrow(position_i) > 1) {
    stop("Integer i occurs more than once in the matrix.")
  }

  # Find position of (i+1) in the matrix
  position_i_plus_one <- which(matrix == (i + 1), arr.ind = TRUE)

  if (nrow(position_i_plus_one) == 0) {
    stop("Integer (i+1) not found in the matrix.")
  } else if (nrow(position_i_plus_one) > 1) {
    stop("Integer (i+1) occurs more than once in the matrix.")
  }

calculate_rho <- function(sigma, rho) {
  # Calculate rho(sigma)
  result <- numeric(length(sigma))
  for (i in 1:(length(sigma) - 1)) {
    temp_sigma <- sigma
    temp_sigma[i:(i+1)] <- sigma[(i+1):i]
    result[i] <- rho(temp_sigma)
  }
  # Apply f(rho) = X
  X <- mean(result)
  return(X)
}

f_rho <- function(rho) {
  f_rho_values <- numeric(factorial(length(rho)))
  permutations <- permutations(length(rho))

  for (i in 1:nrow(permutations)) {
    sigma <- permutations[i, ]
    f_rho_values[i] <- calculate_rho(sigma, rho)
  }

  return(f_rho_values)
}

associate_permutations <- function(n, numbers) {
  if (length(numbers) != factorial(n)) {
    stop("Length of input vector must be equal to n!")
  }

  permutations <- unique(permutations(n))
  associated_list <- list()

  for (i in 1:length(numbers)) {
    associated_list[[i]] <- list(permutation = permutations[[i]], number = numbers[i])
  }

  return(associated_list)
}


#i need to edit this function a little bit more
#not final

calculate_fourier_transform <- function(f) {
  n <- length(f)
  permutations <- permutations(n)
  fourier_transform <- array(0, dim = c(n, n))

  for (i in 1:nrow(permutations)) {
    sigma <- permutations[i, ]
    for (j in 1:n) {
      fourier_transform[j, ] <- fourier_transform[j, ] + f(sigma) * sigma[j]
    }
  }

  return(fourier_transform)
}
