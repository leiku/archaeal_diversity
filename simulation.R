

#3 treatments, 18 replicates：
#control: C = 0.3, strength: 0.2 - 0.4
#RC(randomized reduced connectance):  randomly remove some competition links
#RS (randomized reduced strength): randomly decrease competition strength by 0.1-0.2


library(vegan)
library(mvtnorm)
library(dplyr)
library(tidyverse)
source("Fn_discrete_time_local.R")
source("Fn_comp_matrix.R")


set.seed(1)


####### build up a species pool ########
NN <- 1000   #sp number in the pool
S <- 200  # get 750 species from 1000, and 500 are common species
reps <- 18


# set up growth rate
r0 <- runif(NN, 0.5, 1)    #maximum growth rate

# set up original competition interaction
comp0 <- Fn_comp_matrix(S=NN, C=0.3, min.c = 0.2, max.c=0.4)

library(MASS)
#K0 <- runif(NN, 50, 5000)
K0 <- exp(rnorm(NN, 3, 1))
i.common <- sample(1:NN, S/2, prob=K0, replace = F)
#i.common <- order(K0, decreasing = T)[1:(S/2)]
# sample common species from pool (with a half common species) #


i.NN <- 1:NN
I.SP <- matrix(0, reps, S)
for(i in 1:reps){
  i.rare <- sample(i.NN[-i.common], S-length(i.common), replace = F)
  I.SP[i,] <- c(i.common, i.rare)
}



########## simulation ###########

Treatment <- data.frame(treat=rep(c("control","RS","RC"), each=reps),
                        replicate = rep(1:reps, 3))

B <- matrix(NA, nrow = nrow(Treatment), ncol = NN)  
rownames(B) <- paste0(Treatment$treat,"_", Treatment$replicate)
colnames(B) <- paste0("S", 1:NN)

comp.list <- vector("list", nrow(B))
names(comp.list) <- rownames(B)

for(i in 1:nrow(B)){
  i.sp <- I.SP[Treatment$replicate[i],]
  
  r1 <- r0[i.sp]
  K1 <- K0[i.sp]
  
  comp1 <- comp0[i.sp, i.sp]
  
  ##competition
  if(Treatment$treat[i]=="RC"){
    C <- 0.15
    L <- ceiling(S * (S-1) * C )  #number of links
    i.links <- which(comp1>0 & comp1<1)  #index of original links
    comp <- comp1
    comp[sample(i.links, length(i.links) - L)] <- 0 #remove other links to make sure C
  }else if(Treatment$treat[i] == "RS"){
    comp <- comp1
    i.links <- which(comp1>0 & comp1<1)  #index of original links
    comp[i.links] <- comp1[i.links] - runif(length(i.links), 0.1, 0.2)
  }else{comp <- comp1}
  
  comp.list[[i]] <- comp
  
  len.t <- 1000  # time length
  N0 <- rep(1, S) # original density
  
  ### environmental fluctuation
  rho <- 0
  sigma2.e <- 0.006
  sigma2 <- matrix(sigma2.e*rho, S, S)
  diag(sigma2) <- sigma2.e
  eps <- matrix(rnorm(len.t*S,0, sigma2.e), len.t, S)
  
  para <- list(r=r1, K=K1, c=comp, env=eps)
  N1 <- Fn_discrete_time_local(N0, para, len.t=len.t, model="Kilpatrick")
  N.eq <- colMeans(N1[(len.t-99):len.t, ])
  i.bad <- which(N1[len.t,]==0)
  if(length(i.bad)>0){N.eq[i.bad] <- 0} #set species going extinction as 0
  
  B[i, i.sp] <- N.eq
  print(i)
}
B[is.na(B)] <- 0
saveRDS(list(Treatment=Treatment, Biomass=B, I.SP = I.SP, comp.list=comp.list), 
        "Data_simulation.RDS")






