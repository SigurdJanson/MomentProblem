library(profvis)
#source("../SimPdf.R")


Mn <- function(datasample, n=0) {
  N <- length(datasample) #Sample size
  Mn <- 0 #Initialize moment
  for (i in 1:N){
    Mn <- Mn + (datasample[i])^n #Accumulate moment
  }
  Mn <- Mn/N #Averaging over sample size
  return(Mn)
}
MomentSet <- function(datasample, nmax=10) { 
  moments <- c(1) #Zero-th moment value
  n <- c(0) #Zero-th moment order
  for (i in 1:nmax) {
    moments[i+1] <- Mn(datasample, i) #n-th moment value
    n[i+1] <- i #n-th moment order
  }
  MS <- data.frame(n, moments) #Constructing output data frame
  return(MS)
}



Mn.1 <- function(datasample, n = 0) {
  N <- length(datasample) #Sample size
  Mn <- 0 #Initialize moment
  for (i in 1:N){
    Mn <- Mn + (datasample[i])^n #Accumulate moment
  }
  
  Mn <- Mn/N #Averaging over sample size
  return(Mn)
}
MomentSet.1 <- function(datasample, nmax=10) {
  moments <- c(1) #Zero-th moment value
  n <- c(0) #Zero-th moment order
  #for (i in 1:nmax) {
  
    moments[i+1] <- Mn.1(datasample, i) #n-th moment value
    n[i+1] <- i #n-th moment order
  #}
  MS <- data.frame(n, moments) #Constructing output data frame
  return(MS)
}


dt <- rnorm(1e5)

#r <- profvis({
  print(system.time(MomentSet(dt, nmax=8)))

  print(system.time(MomentSet.1(dt, nmax=8)))
#})

#print(r)