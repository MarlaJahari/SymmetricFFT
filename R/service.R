# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

sn_convolution <- function(FFT1, FFT2) {
  l1<- length(FFT1)
  l2<- length(FFT2)
  L<- min(l1, l2)
  Convolution<- vector("list", L)  # Create a list to hold the 2D arrays

  while (L>0) {
    Convolution[[L]]<- FFT1[[l1]]*FFT2[[l2]]
    l1<- l1-1
    l2<- l2-1
    L<- L-1
  }

  return(Convolution)
}

sn_correlation<- function(FFT1, FFT2) {
  l1<- length(FFT1)
  l2<- length(FFT2)
  L<- min(l1, l2)
  Correlation<- vector("list", L)  # Create a list to hold the 2D arrays

  while(L> 0) {
    Correlation[[L]] <- FFT1[[l1]] * t(FFT2[[l2]])  # Use t() to transpose the array
    l1<- l1-1
    l2<- l2-1
    L<- L-1
  }

  return(Correlation)
}

Counter<-function(N) {
  structure(list(N = N), class = "Counter")
}

sn_multiply <- function(P1, P2) {
  N <- length(P1)
  Prod <- integer(N)  # Create an integer vector to hold the product

  for (i in 1:N) {
    Prod[i] <- P1[P2[i]]
  }

  return(Prod)
}

sn_inverse<- function(P) {
  N<-length(P)
  Inv<-integer(N)  # Create an integer vector to hold the inverse

  for (i in 1:N) {
    Inv[P[i]]<-i
  }

  return(Inv)
}

#calculates the permutation corresponding to a given index for a permutation of size N
index_permutation<-function(N, Index) {
  ccf<- index_ccf(N,Index)
  permutation<-ccf_permutation(ccf)
  return(permutation)
}

sn_p<- function(N) {
  index <- sample(1:factorial(N), 1)  # Generate a random index
  P<- index_permutation(N, index)
  return(P)
}

sn_cc <-function(N,LB,UB) {
  CC <-integer(N)  # an integer vector to hold the cyclically shifted values

  for (i in 1:(LB-1)) {
    CC[i] <-i
  }
  for (i in LB:UB) {
    CC[i] <-i+1
  }
  CC[UB] <-LB

  for (i in (UB+1):N) {
    CC[i] <-i
  }
  return(CC)
}

sn_at<- function(N, K) {
  AT<- integer(N)  #an integer vector to hold the patterned values

  for (i in 1:(K-1)) {
    AT[i] <-i
  }

  AT[K]<-K+1
  AT[K+1]<-K

  for (i in (K+2):N) {
    AT[i]<- i
  }

  return(AT)
}

sn_at<-function(N) {
  k<- sample(1:(N-1),1)  #generate a random value for k
  AT<- integer(N)  #create an integer vector to hold the transposition

  AT[1:k]<- 1:k
  AT[k]<- k + 1
  AT[k+1]<- k

  if (k< N - 1) {
    AT[(k + 2):N] <- (k + 2):N
  }

  return(AT)
}

sn_t<-function(N, I, J) {
  Tr <-integer(N)  #an integer vector to hold the transposition

  for (i in 1:N) {
    Tr[i]<-i
  }

  Tr[I]<-J
  Tr[J]<-I

  return(Tr)
}


sn_t <-function(N) {
  i<-sample(1:N,1)  #generate a random value for i
  j<-sample(1:N,1)  #generate a random value for j
  Tr<-integer(N)  #create an integer vector to hold the transposition

  for (k in 1:N){
    Tr[k]<-k
  }

  Tr[i]<-j
  Tr[j]<-i

  return(Tr)
}

#the contiguous cycle factorization is a way to represent a permutation as a product of adjacent transpositions
#In a permutation, an adjacent transposition<->swapping two adjacent elements.
#CCF represents a permutation as a product of such transpositions that, when applied successively,
#reconstruct the original permutation

#computes the CCF of a permutation P

permutation_ccf <-function(P) {
  N<-length(P)
  CCF<- integer(N - 1)  #create an integer vector to hold the CCF

  i<- 1
  for (j in N:2) {
    CCF[i]<- P[j]
    cc<- sn_cc(N, P[j], j)
    cc_inv<- sn_inverse(cc)
    P<- sn_multiply(cc_inv, P)
    i<- i+1
  }
  return(CCF)
}

#calculates the unique index of a permutation corresponding to a given contiguous cycle factorization
ccf_index <-function(CCF) {
  N<-length(CCF)+1
  Index<-1

  for (i in 1:(N-1)){
    N<- N-1
    if (CCF[i]!=1) {
      Index <-Index +(CCF[i]-1)*factorial(N)
    }
  }
  return(Index)
}

#calculates the unique index that a given permutation P maps to
permutation_index<- function(P) {
  ccf<- permutation_ccf(P)
  index<- ccf_index(ccf)
  return(index)
}

#calculates the contiguous cycle factorization corresponding to a permutation with a given index for a permutation of size N
index_ccf<-function(N, Index) {
  CCF<-integer(N-1)  # Create an integer vector to hold the CCF
  Index<-Index-1

  for(i in 1:(N-1)) {
    q <- floor(Index/factorial(N - i))
    Index <- Index- q*factorial(N - i)
    CCF[i] <- q+1
  }

  return(CCF)
}

#calculates the permutation corresponding to a contiguous cycle factorization
ccf_permutation<-function(CCF){
  N<-length(CCF)+1
  P<-integer(N)  # create an integer vector to hold the permutation

  for(i in 1:N) {
    P[i]<-i
  }
  for (i in 1:(N-1)) {
    cc<-sn_cc(N, CCF[i],N+1-i)
    P<-sn_multiply(P, cc)
  }
  return(P)
}


#big function alert
#calculates the adjacent transposition factorization of a permutation P
permutation_atf<-function(P) {
  N<-length(P)
  CCF<- permutation_ccf(P)
  dATF<- vector("list", length(CCF))
  L<- 0

  for(i in 1:length(CCF)) {
    P<- CCF[[i]]
    D<- N-P
    P<- P-1
    L<- L+D
    N<- N-1
    dATFs<-integer(D)

    for(j in 1:D) {
      dATFs[j]<-P+j
    }
    dATF[[i]]<-dATFs
  }
  ATF<-integer(L)
  index<- 1

  for(i in 1:length(CCF)) {
    for(j in 1:length(dATF[[i]])) {
      ATF[index]<- dATF[[i]][j]
      index<- index+1
    }
  }
  return(ATF)
}

#we'll need the Matrix library for the next functions sparse matrix manipulation
if (!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}

#calculates Young's Orthogonal Representation of a permutation P corresponding to a partition using provided matrices YORnp
yor_permutation<-function(P,YORnp){
  ATF<-permutation_atf(P)

  if(length(ATF)==0){
    Dim<-nrow(YORnp[[1]])
    RM<-Diagonal(Dim)
  } else{
    RM<-Matrix(1, nrow=nrow(YORnp[[ATF[1]]]), ncol=ncol(YORnp[[ATF[1]]]), sparse=TRUE)

    for (i in 1:length(ATF)){
      RM<-RM %*% YORnp[[ATF[i]]]
    }
  }
  return(RM)
}


