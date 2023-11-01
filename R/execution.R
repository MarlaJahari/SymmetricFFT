#provided with the problem size, a partition, and a permutation
#computes the young's orthogonal representation of the specified permutation for the given partition
#prints the degree of the yamanouchi Symbols and the representation matrix for the permutation within the given partition

example2<-function(N,partition,p1,p2){
  #call the partitions function and get P and WI
  result <- partitions(N)
  P <- result[[1]]
  WI <- result[[2]]

  p<-1
  while (!setequal(P[[p]], partition)) {
    p<-p+1
  }

  #call the sn_multiply function
  pm <- sn_multiply(p1, p2)

  #call the yor function and get YOR and PT
  result <- yor(N)
  YOR <- result$YOR
  PT <- result$PT

  YORnp <- YOR[[N]][[p]]

  # Call the yor_permutation function for p1 and p2
  RM1 <- yor_permutation(p1, YORnp)
  RM2 <- yor_permutation(p2, YORnp)

  # Call the yor_permutation function for pm
  RMm <- yor_permutation(pm, YORnp)

  error <- RM1 %*% RM2 - RMm

  cat("If a and b are permutations and c = a * b, demonstrates that the representation of c is the representation of a multiplied by the representation of b\n\n")

  ST1 <- permutation_string(p1)
  cat("Representation of ", ST1, "\n")
  print(round(RM1, 5))
  cat("\n")

  ST2 <- permutation_string(p2)
  cat("Representation of ", ST2, "\n")
  print(round(RM2, 5))
  cat("\n")

  STm <- permutation_string(pm)
  cat(ST1, " * ", ST2, " is ", STm, "\n\n")

  cat("Representation of ", STm, "\n")
  print(round(RMm, 5))
  cat("\n\n")

  cat("Maximum error: \n")
  print(max(abs(error)))
}

example2(5, c(3, 2), c(1, 2, 4, 3, 5), c(2, 3, 1, 4, 5))


