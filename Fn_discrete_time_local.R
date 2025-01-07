# model = "Wang": Wang and Loreau, 2016, Ecology Letters
# model = "Kilpatrick": Kilpatrick, 2003, Nature


Fn_discrete_time_local <- function(N0, para, len.t=10000, model="Wang"){
  
  #N0 should be a matrix with 1*S
  K <- para$K # 1*S, carrying capacity
  i.good <- which(K>1)
  
  K <- K[i.good]
  r <- para$r[i.good] # 1*S, growth rate
 
  S <- length(r)  #number of species
  
  
  #c[i,j] means competition coefficient of j on i 
  c <- para$c[i.good, i.good] # a matrix of s*s; notice the diagonal is 1
  
  env <- para$env[,i.good]   #environmental fluctuation, a matrix with t * S
  
  Nt <- N0[i.good]    #to store the temporal abundance
  Nt1 <- Nt   #abundance at the next time step
  N <- matrix(NA, len.t, S)   #to store the abundance of sp.
  N[1,] <- N0[i.good]
  
  if(model=="Wang"){
    for(i in 1:(len.t-1)){
      det.t <- 0.1
      for(t1 in 1:(1/det.t)){
        competition <- c %*% Nt / K
        growth <- r * Nt * (1 - competition ) * det.t  #logistic growth term
        environ <- Nt * env[i,] * sqrt(det.t)   #environmental
        Nt1 <- Nt + growth +  environ
        Nt <- Nt1
        Nt[Nt<1e-10] <- 0
        if(sum(Nt>K)>0){Nt[Nt>K] <- K[Nt>K]} #make sure not higher than K
      }
      N[i+1, ] <- Nt
    }
  }else if(model=="Kilpatrick"){
    for(i in 1:(len.t-1)){
      competition <- c %*% N[i,]/K
      competition[competition<0] <- 0  #make sure growth rate not higher than r
      Nt2 <- N[i,] * exp(r * (1 - competition)+env[i,])
      Nt2[Nt2<1e-10] <- 0
      if(sum(Nt2>K)>0){Nt2[Nt2>K] <- K[Nt2>K]} #make sure not higher than K
      N[i+1,] <- Nt2
    }
  }else if(model=="Kilpatrick.mod"){  #
    for(i in 1:(len.t-1)){
      competition <- c %*% N[i,]/K
      competition[competition<0] <- 0  #make sure growth rate not higher than r
      Nt2 <- N[i,] * exp(r * (1 - competition)+env[i,])
      Nt2[Nt2<1e-10] <- 0
      if(sum(Nt2>K)>0){Nt2[Nt2>K] <- K[Nt2>K]} #make sure not higher than K
      N[i+1,] <- Nt2
    }
  }
  
  N.all <- matrix(0, nrow=len.t, ncol=length(N0))
  N.all[,i.good] <- N
  return(N.all)
}

