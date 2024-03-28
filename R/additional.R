#find the indices of partitions in the Pd array corresponding to the partitions in the PDA array
#searching for each partition in PDA within the Pd array

dia1<- function(PDA, Pd, lb) {
  PDA_L<- length(PDA)
  DIA<- numeric(PDA_L)
  c<- 1
  p<- lb
  while(TRUE) {
    if(identical(Pd[[p]], PDA[[c]])) {
      DIA[c]<- p
      c<- c+1
      if(c>PDA_L){
        return(DIA)
      }
    }
    p<- p+1
  }
}

dia2<-function(PDA, Pd, lb1, lb2) {
  PDA_L<-length(PDA)
  DIA<-integer(PDA_L)

  p<-lb1
  while (TRUE) {
    if(identical(Pd[[p]],PDA[[1]])) {
      DIA[1]<-p
      break
    }
    p<- p+1
   }
  p<-lb2
  c<-2
  while(TRUE) {
    if (identical(Pd[[p]], PDA[[c]])) {
      DIA[c]<-p
      c<-c+1
      if(c>PDA_L) {
        return(DIA)
      }
    }
    p<-p+1
  }
}

#generating a partition decomposition array PDA from a given partition P
#detects where the partition P splits into smaller partitions based on the values in the partition
pda<- function(P) {
  P_L<-length(P)
  num<-1
  for(i in 1:(P_L - 1)) {
    if(P[i] > P[i + 1]) {
      num<- num + 1
    }
  }
  PDA<- vector("list", num)
  c<- 1
  for(i in 1:(P_L - 1)) {
    if(P[i]>P[i + 1]) {
      D<-as.integer(P)
      D[i]<- D[i]-1
      PDA[[c]]<-D
      c <- c+1
    }
  }
  if(P[P_L]==1) {
    D<- as.integer(P[1:(P_L - 1)])
    PDA[[c]] <- D
  } else{
    D<- as.integer(P)
    D[P_L]<- D[P_L] - 1
    PDA[[c]]<- D
  }
  return(PDA)
}

#generates a partition tree PT based on the input data P and WI
#the PT structure represents how each partition in level n decomposes into partitions at level  n-1

partition_tree<-function(N,P,WI){
  PT<-vector("list",N)
  PT[[1]]<-vector("list",0)

  for(n in 2:N){
    Pd<-P[[n-1]]
    Pn<-P[[n]]
    PTn<-vector("list", length(Pn))

    for(p in 2:length(Pn)){
      PDA<-pda(Pn[[p]])
      if(length(PDA)==1||PDA[[1]][1]==PDA[[2]][1]) {
        lb<-1

        if(PDA[[1]][1]!=1) {
          lb<-WI[n-1, PDA[[1]][1]-1]+1
        }

        PTn[[p]]<- dia1(PDA,Pd,lb)
       }else{
        lb1<-1

        if (PDA[[1]][1]!=1) {
          lb1<-WI[n-1, PDA[[1]][1]-1]+1
        }

        lb2<- WI[n-1, PDA[[2]][1]-1]
        PTn[[p]]<-dia2(PDA, Pd, lb1, lb2)
      }
    }

    PT[[n]]<-PTn
  }

  return(PT)
}

#generates partitions of N and gives width in the form of WI.
#resulting partitions are organized in a nested list structure


partitions<- function(N) {
  part<- list()
  WI <- matrix(0, N, N)
  for (n in 1:N) {
    WI[n, 1] <- 1
  }
  for (w in 2:N-1) {
    WI[w, w] <- 1
    for (n in (w + 1):N) {
      num <- 0
      for (i in 1:min(n - w, w)) {
        num <- num + WI[n - w, i]
      }
      WI[n, w] <- num
    }
  }
  for (n in 2:N) {
    for (w in 2:n){
      WI[n, w] <- WI[n, w] + WI[n, w - 1]
    }
  }
  WI[N,N]<-WI[N,N]+1

  generate_partitions <- function(remaining, current_partition) {
    if (remaining == 0) {
      part <<- c(part, list(current_partition))
    } else if (remaining > 0) {
      for (i in 1:min(remaining, ifelse(length(current_partition) == 0, remaining, current_partition[length(current_partition)]))) {
        generate_partitions(remaining - i, c(current_partition, i))
      }
    }
  }

  generate_partitions(N, c())
  return(list(part,WI))
}

#new coms added
