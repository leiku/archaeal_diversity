
# function to create competition matrix
#S: number of species
#C: connectance


Fn_comp_matrix <- function(S, C, min.c=0.1, max.c=0.5){

   cm <- matrix(0, S, S)
   L <- ceiling(S * (S-1) * C )  #number of links
   i.off <- which(upper.tri(cm) | lower.tri(cm))  #index of upper or lower triangle in matrix
   
   cm[sample(i.off, L)] <- runif(L, min.c , max.c)   #undirected
   diag(cm) <- 1

  
 return(cm) 
}