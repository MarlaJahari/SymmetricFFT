#calculates the direct sum matrix (DSM) based on the input parameters Dim, sFFTn, and PTnp

dsm <-function(Dim,sFFTn,PTnp) {
  DSM<-matrix(0, nrow=Dim, ncol=Dim)
  offset<-0

  for(i in 1:length(PTnp)) {
    FC<-sFFTn[[PTnp[i]]]
    sDim<-nrow(FC)

    for(r in 1:sDim){
      for(c in 1:sDim){
        v<- FC[r, c]
        if (v!= 0.0) {
          DSM[offset+r, offset+c]<-v
        }
       }
    }

    offset<-offset+sDim
  }
  return(DSM)
}

#calculates the fourier coefficient (FC) corresponding to a partition P of N
fc<-function(N, YORnp, PTnp, sFFT){

  Dim<-nrow(YORnp[[1]])
  FC<-dsm(Dim, sFFT[[N]], PTnp)
  CCM<-Diagonal(Dim)

  for(n in (N-1):1) {
    CCM<-YORnp[[n]] %*% CCM
    DSM<-dsm(Dim, sFFT[[n]], PTnp)
    FC_n<-CCM %*% DSM
    FC<-FC + FC_n
  }
  return(FC)
}

# fetches data from remote references (RR_YOR, RR_PT, RR_sFFT) and uses the fetched data to calculate the fourier coefficient
fc_remote <- function(N, p, YOR, PT, sFFT) {
  FC<- fc(N, YOR[[N]][[p]], PT[[N]][[p]], sFFT)
  return(FC)
}

#calculates the fast fourier transform (FFT) for a group based on the input parameters N, YORn, PTn, and sFFT
combine_sfft<-function(N, YORn, PTn, sFFT) {
  NP<-length(YORn)
  FFT<-list()

  for(p in 1:NP) {
    YORnp<-YORn[[p]]
    PTnp<-PTn[[p]]
    FC<-fc(N, YORnp, PTnp, sFFT)
    FFT[[p]]<-FC
  }

  return(FFT)
}

# calculates the fast fourier transform (FFT) based on the input parameters N, SNF, YOR, PT, and C
compute_fft<- function(N, SNF, YOR, PT, C) {
  if(N==1) {
    sFFT<-list(matrix(SNF[C$N], nrow = 1))
    C$N<-C$N+1
    return(sFFT)
  }
  sFFT<-vector("list", N)
  for(n in 1:N) {
    sFFT[[n]]<- compute_fft(N - 1, SNF, YOR, PT, C)
  }
  FFT<- combine_sfft(N, YOR[[N]], PT[[N]], sFFT)
  return(FFT)
}

#we wont really need this, will delete later
#function that fetches data from remote references (RR_YOR and RR_PT) and uses the fetched data to calculate the fast fourier transform (FFT)
compute_fft_remote<-function(N,SNF,YOR,PT,C){
  FFT<-compute_fft(N, SNF, YOR, PT, C)
  return(FFT)
}

# Parameters:
# N: Int
# SNF: Array of Float64
# YOR: Array of YOR data
# PT: Array of PT data
# Returns: FFT result

#implement the Fast Fourier Transform (FFT) algorithm without parallel computing, following the Clausen's algorithmic approach
sn_fft <- function(N, SNF, YOR, PT) {
  if (N==1) {
    #base case
    FFT<-matrix(SNF, nrow = 1)
  } else{
    sFFT<-vector("list", N)
    for(n in 1:N){
      sFFT[[n]]<-sn_fft(N-1, SNF, YOR, PT)
    }
    FFT<-combine_sfft(N, YOR[[N]], PT[[N]], sFFT)
  }
  return(FFT)
}

# Parameters:
#	N - the problem size
#	PA - a list of permutations, each represented as a vector of integers
#	VA - a vector of values associated with each permutation
# Return Values:
#	SNF - a vector of values associated with permutations indexed by permutation_index()
#	- This is the format for the SNF parameter of sn_fft()

#	- any permutation of N not represented in PA will be assigned a value of zero


#returns a vector SNF that represents the signature normal form of permutations
snf<-function(N, PA, VA) {
  SNF<-numeric(factorial(N))
  for(i in 1:length(PA)) {
    Permutation<-PA[[i]]
    Index<-permutation_index(Permutation)
    SNF[Index]<-VA[i]
  }
  return(SNF)
}


#calculates the number of Standard Tableaux of the Young Diagram of a given partition P with a given size N
degree<-function(N, P) {
  if (N>20){
    NST<-2432902008176640000
    for(i in 1:length(P)) {
      for(j in 1:P[i]) {
        NST<- NST/hook_length(P,i,j)
      }
    }
    for(n in 21:N) {
      NST<- NST*n
    }
    NST<-round(NST)
    return(NST)
  } else {
    NST<-factorial(N)
    NST<-NST/hook_product(P)
    NST<-round(NST)
    return(NST)
  }
}

#calculates the hook product of a given partition P
hook_product<-function(P) {
  HP<- 1
  for(i in 1:length(P)) {
    for(j in 1:P[i]) {
      HP<- HP*hook_length(P,i,j)
    }
  }
  return(HP)
}

#calculates the hook length of a box in the i_th row and j_th column of the young diagram of a given partition P
hook_length<-function(P,i,j) {
  R<- P[i]-j
  B<-0
  i<-i+1
  while(i<=length(P) && j<=P[i]) {
    i<- i+1
    B<- B+1
  }
  HL<- 1+R+B
  return(HL)
}



