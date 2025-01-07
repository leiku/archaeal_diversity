
#### get population matrix #####  
library(vegan)
library(ggplot2)
library(patchwork)
library(Hmisc)  #for rcorr
library(igraph)

colors <- c(light="deepskyblue", heavy="brown1")


d.OTU <- readRDS("Data_OTU_Archaea.RDS")
m.OTU <- d.OTU$OTU[, -(1:4)]
m.OTU <- t(m.OTU)
colnames(m.OTU) <- d.OTU$OTU$Taxon.name

design <- d.OTU$design
#re-assign index: H for heavy, L for light; A for BH, B for HX, C for YQ
design$Site[design$Site=="BH"] <- "A"
design$Site[design$Site=="HX"] <- "B"
design$Site[design$Site=="YQ"] <- "C"
design$Treat <- paste0(rep(rep(c("H","L"),each=6),3), rep(c("A","B","C"), each=12),
                       design$Replicate)

rownames(m.OTU) <- design$Treat
m.OTU <- m.OTU[order(design$Salt),]
m <- m.OTU / rowSums(m.OTU)
design <- design[order(design$Salt),]
design$Salt

#get favor of genus
d.P <- data.frame(ID=colnames(m), P=rep(0, ncol(m)),
                  avg.light=rep(0, ncol(m)), avg.heavy=rep(0, ncol(m)))
for(i in 1:ncol(m)){
  z <- t.test(m[1:18, i], m[19:36, i])
  d.P$P[i] <- z$p.value
  d.P[i, 3:4] <- z$estimate
}
d.P$Favor <- "none"
d.P$Favor[d.P$P<0.05 & (d.P$avg.light>d.P$avg.heavy)] <- "light"
d.P$Favor[d.P$P<0.05 & (d.P$avg.light<d.P$avg.heavy)] <- "heavy"
table(d.P$Favor)

d.P <- d.P[,c(1,5, 2,3,4)]


#### overall  ##########

i1 <- which(design$Salt == "light")
m1 <- as.matrix(m[i1,])
n.zero <- apply(m1, 2, function(x) sum(x==0))
m1 <- m1[, n.zero<=9]
#m1 <- m1[, colMeans(m1)>0.0001]
tmp <- rcorr(m1, type="spearman")
R1 <- tmp$r
diag(R1) <- 0
P1 <- tmp$P
R1[P1>0.05] <- 0  # That means, R[abs(R) < 0.47] <- 0
sum(R1!=0)/nrow(R1)^2
G1 <- graph_from_adjacency_matrix(R1, mode="undirected", weighted=T, diag=F)
E1 <- as.data.frame(as_edgelist(G1))
E1$weight <- E(G1)$weight
colnames(E1) <- c("Source", "Target", "Weight")
E1$Type="Undirected"
E1$Id <- 0:(nrow(E1)-1)
E1$Cor <- "Positive"
E1$Cor[E1$Weight<0] <- "Negative"
E1 <- E1[,c(1,2,4,5,6,3)]

write.csv(E1, "Network/Correlation_light.csv", row.names = F)


d.P1 <- d.P[d.P$Label %in% colnames(R1),]
write.csv(d.P1, "Network/Node_favor_light.csv", row.names = F)



i2 <- which(design$Salt == "heavy")
m2 <- as.matrix(m[i2,])
n.zero <- apply(m2, 2, function(x) sum(x==0))
m2 <- m2[, n.zero<=9]
#m2 <- m2[, colMeans(m2)>0.0001]
tmp <- rcorr(m2, type="spearman")
R2 <- tmp$r
diag(R2) <- 0
P2 <- tmp$P
R2[P2>0.05] <- 0  # That means, R[abs(R) < 0.47] <- 0
sum(R2!=0)/nrow(R2)^2
G2 <- graph_from_adjacency_matrix(R2, mode="undirected", weighted=T, diag=F)
E2 <- as.data.frame(as_edgelist(G2))
E2$weight <- E(G2)$weight
colnames(E2) <- c("Source", "Target", "Weight")
E2$Type="Undirected"
E2$Id <- 0:(nrow(E2)-1)
E2$Label <- "Positive"
E2$Label[E2$Weight<0] <- "Negative"
E2 <- E2[,c(1,2,4,5,6,3)]
write.csv(E2, "Network/Correlation_heavy.csv", row.names = F)

d.P2 <- d.P[d.P$Label %in% colnames(R2),]
write.csv(d.P2, "Network/Node_favor_heavy.csv", row.names = F)




#### each site #######

sites <- c("A","B","C")
salts <- c("light", "heavy")
res.net <- data.frame(Site=rep(sites, each=2), Salt=rep(salts, 3),
                      P.neg=rep(NA,6), L=rep(NA,6), C=rep(NA,6), S=rep(NA,6))
res.net$Salt <- factor(res.net$Salt, levels=c("light","heavy"))
k <- 1
for(i in 1:3){
  for(j in 1:2){
    ii <- which(design$Site==sites[i] & design$Salt == salts[j])
    m1 <- as.matrix(m[ii,])
    n.zero <- apply(m1, 2, function(x) sum(x==0))
    #m1 <- m1[, colMeans(m1)>0.0001]
    m1 <- m1[, n.zero<=3]
    tmp <- rcorr(m1, type="spearman")
    R <- tmp$r
    diag(R) <- 0
    P <- tmp$P
    R[P>0.05] <- 0  # That means, R[abs(R) < 0.81168] <- 0
    #R[abs(R) < 0.6] <- 0
    R1 <- abs(R)
    R1[R1>0] <- 1
      
    res.net$Site[k] <- sites[i]
    res.net$Salt[k] <- salts[j]
      
    n.neg <- sum(R<0)
    res.net$P.neg[k] <- n.neg / sum(R!=0)
      
    nn <- nrow(R)
    res.net$L[k] <- sum(R1) / 2 / nn
    res.net$C[k] <- sum(R1) / nn / (nn-1)
    res.net$S[k] <- nn
    
    ### edge list for Gephi
    G <- graph_from_adjacency_matrix(R, mode="undirected", weighted=T, diag=F)
    EL <- as.data.frame(as_edgelist(G))
    EL$weight <- E(G)$weight
    colnames(EL) <- c("Source", "Target", "Weight")
    EL$Type="Undirected"
    EL$Id <- 0:(nrow(EL)-1)
    EL$Cor <- "Positive"
    EL$Cor[EL$Weight<0] <- "Negative"
    #write.csv(EL, paste0("Network/Edge_",sites[i],"_",salts[j],".csv"), row.names = F)
    
    d.P2 <- d.P[d.P$ID %in% colnames(R),]
    #write.csv(d.P2, paste0("Network/Node_",sites[i],"_",salts[j],".csv"), row.names = F)
     
    k <- k+1
  }
}

t.test(res.net$P.neg[c(1,3,5)], res.net$P.neg[c(2,4,6)], paired=T)
t.test(res.net$L[c(1,3,5)], res.net$L[c(2,4,6)], paired=T)
t.test(res.net$C[c(1,3,5)], res.net$C[c(2,4,6)], paired=T)
ggplot(res.net, aes(x=Salt, y=P.neg, fill=Salt)) + 
  geom_boxplot() +
  annotate("text", x=1.5, y=0.48, label="NS") +
  labs(x="", y="Proportion of negative links")+
  scale_fill_manual(values=colors) +
  theme_bw()
ggplot(res.net, aes(x=Salt, y=C, fill=Salt)) + 
  geom_boxplot() +
  annotate("text", x=1.5, y=0.19, label="NS") +
  labs(x="", y="Connectance")+
  scale_fill_manual(values=colors) +
  theme_bw()

