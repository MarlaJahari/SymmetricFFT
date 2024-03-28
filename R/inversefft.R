#update the necessary element in an array of FFTs for left cosets.
#fill in the FFT components for the left coset determined by n
#coms added

update_sfft <- function(N, n, sFFT, FFTp, YORnp, YORd, PTnp) {
  Dim<-dim(YORnp[[1]])[1]
  M<-diag(Dim)
  for (ccn in n:(N - 1)) {
    M<-M %*% YORnp[[ccn]]
  }
  M<-t(M)
  M<-M %*% FFTp
  lb<-1
  for (d in 1:length(PTnp)) {
    index<-PTnp[[d]]
    sDim<-dim(YORd[[index]])[1]
    ub<-lb+sDim -1
    sFFT[[index]]<-sFFT[[index]] + (Dim/(sDim*N))*M[lb:ub, lb:ub]
    lb<-ub+1
  }
}

#computing the inverse fast fourier transform (IFFT) for functions defined over the symmetric group Sn
#populates an array SNF with values representing the IFFT of the input (as in the past tense) FFT components
#recursively computes the IFFT for left cosets and partitions using young's orthogonal representations.

compute_ifft<-function(N,SNF,FFT,YOR,PT,C){
  if(N==2) {
    YORn<-YOR[[2]]

    #do not let %*% intimidate you, its the matrix multiplication operator

    sFFT1<- 0.5*(YORn[[1]][[1]] %*% FFT[[1]]+YORn[[2]][[1]] %*% FFT[[2]])
    SNF[C$N]<-sFFT1[1, 1]
    C$N<- C$N+1
    sFFT2<- 0.5*(FFT[[1]]+FFT[[2]])
    SNF[C$N]<-sFFT2[1,1]
    C$N<- C$N+1
  } else {
    YORn<-YOR[[N]]
    NPn<-length(YORn)
    YORd<-YOR[[N-1]]
    NPd<-length(YORd)
    PTn<-PT[[N]]
    sFFT<-vector("list", length=NPd)

    for(n in 1:N) {
      for(p in 1:NPd) {
        Dim<-nrow(YORd[[p]][[1]])
        sFFT[[p]]<-matrix(0,Dim,Dim)
      }
      for(p in 1:NPn){
        sFFT[[p]]<-update_sfft(N,n,sFFT[[p]],FFT[[p]],YORn[[p]],YORd,PTn[[p]])
      }

      compute_ifft(N-1,SNF,sFFT,YOR,PT,C)
    }
  }
}

compute_sifft<-function(N, n, FFT, YOR, PT) {
  pSNF<-numeric(factorial(N - 1))
  C<-1

  YORn<-YOR[[N]]
  NPn<-length(YORn)
  YORd<-YOR[[N - 1]]
  NPd<-length(YORd)
  PTn<-PT[[N]]

  sFFT<-vector("list", NPd)

  for(p in 1:NPd) {
    Dim<-dim(YORd[[p]][[1]])[1]
    sFFT[[p]]<-matrix(0, nrow=Dim, ncol=Dim)
  }

  for(p in 1:NPn) {
    update_sfft(N,n,sFFT,FFT[[p]],YORn[[p]],YORd,PTn[[p]])
  }

  compute_ifft(N-1,pSNF,sFFT,YOR,PT,C)
  return(pSNF)
}

#calculate the inverse fast fourier transform of a given input fast fourier transform (FFT)
#using tools in algebraic combinatorics and representation theory

sn_ifft<- function(N, FFT, YOR, PT) {
  np<- 1
  if (np==1 || N<10) {
    SNF<-numeric(factorial(N))
    compute_ifft(N, SNF, FFT, YOR, PT, list(N = 1))
    return(SNF)
  } else {
    pSNFA<- vector("list", N)
    RR_FFT<- vector("list", np)
    RR_YOR<- vector("list", np)
    RR_PT<- vector("list", np)
    for(p in 1:np) {
      if(p != myid()) {
        RR_FFT[[p]]<- FFT
        RR_YOR[[p]]<- YOR
        RR_PT[[p]]<- PT
      }
    }
    i<- 1
    nextidx<-function(){
      idx<- i
      i <<- i + 1
      return(idx)
    }
    for (p in 1:np) {
      if (p!= myid()) {
        while (TRUE) {
          idx<-nextidx()
          if (idx > N) {
            break
          }
          pSNFA[[idx]]<-compute_sifft(N, idx, RR_FFT[[p]], RR_YOR[[p]], RR_PT[[p]])
        }
      }
    }
    SNF<-numeric(factorial(N))
    BS<-factorial(N - 1)
    i<-1
    for(n in 1:N) {
      pSNF<-pSNFA[[n]]
      for(si in 1:BS) {
        SNF[i]<- pSNF[si]
        i<- i+1
      }
    }
    return(SNF)
  }
}

