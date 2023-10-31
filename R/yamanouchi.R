#checks if YSymbol2 is the result of applying the adjacent transposition (K, K+1) to YSymbol1
#returns 1 if true, 0 else

istranspose<-function(YSymbol1, YSymbol2, K) {
  for(i in 1:(K - 1)) {
    if(YSymbol1[i]!=YSymbol2[i]) {
      return(0)
    }
  }

  if(YSymbol1[K]!=YSymbol2[K+1] || YSymbol1[K+1]!=YSymbol2[K]) {
    return(0)
  }

  for(i in (K+2):length(YSymbol1)) {
    if (YSymbol1[i]!=YSymbol2[i]) {
      return(0)
    }
  }

  return(1)
}

#computes the index array IA, the result of applying adjacent transpositions to the Yamanouchi symbols stored in the YSymbols array

ys_indice<-function(N,YSymbols,DA) {
  L<-length(YSymbols)
  IA<-matrix(0, nrow=L, ncol=N-1)
  PF<-matrix(0, nrow=L, ncol=N-1)

  for (i in 1:L) {
    YSymbol<-YSymbols[[i]]
    k<-1
    ti<- i+1

    while(k<N){
      D<- DA[i,k]

      if(D!=1 && D!=-1 && PF[i,k]==0){
        while(ti<= L) {
          if(istranspose(YSymbol, YSymbols[[ti]], k)==1){
            IA[i,k]<- ti
            IA[ti,k]<- i
            PF[i,k]<- 1
            PF[ti,k]<- 1
            break
          }
          ti<-ti+1
        }
      }
      k<-k+1
    }
  }

  return(IA)
}

#computes the distance array DA based on the content array CA
#the distance array represents the distances of the i_th standard tableau for each adjacent transposition (k,k+1)
#k ranges from 1 to N-1

ys_distance <- function(CA) {
  L<- nrow(CA)
  N<- ncol(CA)
  DA<- matrix(0, nrow = L, ncol = N - 1)

  for(i in 1:L) {
    for(k in 1:(N - 1)) {
      DA[i,k]<-CA[i, k+1]-CA[i,k]
    }
  }
  return(DA)
}

#computes the content of each element for a set of yamanouchi symbols associated with a given partition
#content of an element represents the number of elements in the same row with smaller row numbers in the yamanouchi symbol

ys_content<-function(N, R, YSymbols) {
  num_symbols<-length(YSymbols)
  CA<-matrix(0, nrow = num_symbols, ncol = N)

  for(i in 1:num_symbols) {
    YSymbol<-YSymbols[[i]]
    PF<-rep(0, R)

    for(n in 1:N) {
      row<-YSymbol[[n]]
      PF[row]<- PF[row]+1
      CA[i, n]<- PF[row]-row
    }
  }

  return(CA)
}

#to obtain both the distance array (DA) and the index array (IA) for a given set of Yamanouchi symbols
#corresponding to a partition

ys_information<-function(N,R,YSymbols){
  CA<- ys_content(N,R,YSymbols)
  DA<- ys_distance(CA)
  IA<- ys_indices(N, YSymbols, DA)
  return(list(DA = DA, IA = IA))
}

#big function alert
#actual yamagouchi func incoming
#calculates yamanouchi symbols for partitions of a given size N
#calculations are based on the partition tree and decomposition information

ys_symbols<-function(N,P,PT){
  YS<-vector("list", N)

  # initialize yamanouchi symbols for N=1
  YSn<-list(list(list(1)))
  YS[[1]]<-YSn

  for(n in 2:N) {
    Pn<-P[[n]]
    PTn<-PT[[n]]
    Pn_L<-length(Pn)
    YSn<- vector("list", Pn_L)

    for(p in 1:Pn_L) {
      Pnp<-Pn[[p]]
      Pnp_L<-length(Pnp)
      PTnp<-PTn[[p]]
      PTnp_L<-length(PTnp)
      c<-1
      row<-1
      YSnp<-vector("list", degree(n, Pnp))

      for (d in 1:PTnp_L) {
        while(row < Pnp_L && Pnp[[row]] <= Pnp[[row + 1]]) {
          row<- row+1
        }
        YSd<- YS[[n-1]][[PTnp[[d]]]]
        YSd_L<- length(YSd)

        for (i in 1:YSd_L) {
          YSnpc<-integer(n)
          YSdi<-YSd[[i]]

          for (j in 1:(n - 1)) {
            YSnpc[j]<-YSdi[j]
          }
          YSnpc[n]<-row
          YSnp[[c]]<-YSnpc
          c<- c+1
        }
        row<-row+1
      }
      YSn[[p]]<-YSnp
    }
    YS[[n]]<-YSn
  }
  return(YS)
}

#compute youngs orthogonal representation (YOR) for the adjacent transposition (K, K+1) for a given partition P

yor_pk<-function(DA,IA,L,k){
  n<-L
  for(i in 1:L) {
    if(IA[i, k]!=0){
      n<- n+1
    }
  }

  colptr<- integer(L+1)
  colptr[1]<-1
  rowval<-integer(n)
  nzval<-double(n)
  cp<-1

  for (i in 1:L) {
    D<-1/DA[i, k]
    I<-IA[i, k]

    if(I==0){
      colptr[i+1]<-colptr[i]+1
      rowval[cp]<-i
      nzval[cp]<-D
      cp<-cp+1
    } else{
      colptr[i+1]<-colptr[i]+2

      if (I>i){
        rowval[cp]<-i
        nzval[cp]<-D
        cp<-cp+1
        rowval[cp]<-I
        nzval[cp]<- sqrt(1- D*D)
        cp<- cp+1
      } else {
        rowval[cp]<-I
        nzval[cp]<- sqrt(1 - D*D)
        cp <-cp+1
        rowval[cp] <- i
        nzval[cp] <- D
        cp <- cp + 1
      }
    }
  }

  YORpk <- Matrix::sparseMatrix(i = rowval, p = colptr, x = nzval, dims = c(L, L))

  return(YORpk)
}


#calculates youngs orthogonal representation for a given partition P using yamanouchi symbols

yor_p<- function(N, R, YSymbols) {
  YORp<- vector("list", N - 1)
  L<-length(YSymbols)
  info<-ys_information(N, R, YSymbols)
  DA<-info$DA
  IA<-info$IA

  for (k in 1:(N - 1)) {
    YORp[[k]]<-yor_pk(DA, IA, L, k)
  }

  return(YORp)
}

#computes young's orthogonal representations (YOR) for partitions of a given problem size N

yor<- function(N) {

  #generate partitions and width information
  partitions <- partitions(N)
  P <- partitions[[1]]
  WI <- partitions[[2]]

  # Compute the Partition Tree
  PT <- partition_tree(N, P, WI)

  # Generate Yamanouchi Symbols
  YS <- ys_symbols(N, P, PT)

  # Initialize the YOR array
  YOR <- list()

  # Iterate over each value of n
  for (n in 1:N) {
    YSn <- YS[[n]]
    Pn <- P[[n]]
    Pn_L <- length(Pn)
    YORn <- list()

    # Compute YOR representations for each subpartition
    for (p in 1:Pn_L) {
      YSnp <- YSn[[p]]
      Pnp <- Pn[[p]]
      R <- length(Pnp)
      YORn[[p]] <- yor_p(n, R, YSnp)
    }

    # Store YOR representations for the current n
    YOR[[n]] <- YORn
  }

  # Return the YOR array and the Partition Tree
  return(list(YOR = YOR, PT = PT))
}

