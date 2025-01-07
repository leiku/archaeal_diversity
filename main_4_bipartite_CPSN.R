
library(bipartite)
library(ggplot2)
library(patchwork)
library(tidyr)
library(nlme)

colors <- c(light="deepskyblue", heavy="brown1")

### abundance of OTU  ##########
d.OTU <- readRDS("Data_OTU_Archaea.RDS")
m.OTU <- d.OTU$OTU[, -(1:4)]
m.OTU <- t(m.OTU)
colnames(m.OTU) <- d.OTU$OTU$Taxon.name

design <- d.OTU$design
design$Salt <- as.character(design$Salt)
#re-assign index: H for heavy, L for light; A for BH, B for HX, C for YQ
design$Site[design$Site=="BH"] <- "A"
design$Site[design$Site=="HX"] <- "B"
design$Site[design$Site=="YQ"] <- "C"
design$Treat <- paste0(rep(rep(c("H","L"),each=6),3), rep(c("A","B","C"), each=12),
                       design$Replicate)

rownames(m.OTU) <- design$Treat
m.OTU <- m.OTU / rowSums(m.OTU)
m.OTU <- m.OTU[, order(colMeans(m.OTU), decreasing = T)]

###### abundance of gene  ##########
d.gene <- readRDS("Data_gene_taxonomy_Archaea.RDS")

m.gene <- as.matrix(d.gene$gene[,-c(1:6)]) 
info.gene <- d.gene$gene[, 1:6]

#create a new column to tell the difference between the same gene in different groups
info.gene$g <- paste(rep(c("C","N","P","S"), c(272,86,286,240)), info.gene$gene_name, sep="_")  

m1 <- apply(m.gene, 2, function(x) tapply(x, info.gene$g, sum))  #get sum for the same gene
tmp <- strsplit(rownames(m1), "_")
info <- data.frame(group=rep(NA, nrow(m1)), gene=rep(NA, nrow(m1)))
info$group <- unlist(lapply(tmp, "[[", 1))
info$gene <- unlist(lapply(tmp, "[[", 2))
info$name <- rownames(m1)

gene.mat <- t(m1)
m.gene <- gene.mat / rowSums(gene.mat)
rownames(m.gene) <- design$Treat
info <- info[order(colMeans(m.gene), decreasing = T),]
m.gene <- m.gene[,order(colMeans(m.gene), decreasing = T)]

######################  Average across light / heavy  #################

#### get data ready ####

library(readxl)
gene.path <- read_xlsx("Res_Gene_function_Archaea.xlsx", sheet=1)
gene.path <- gene.path[gene.path$path != "NA",]

M.merge <- readRDS("Matrix_OTU_gene_Archaea.RDS")
#res.fit <- read.csv("res_fitting_log-abundance-vs-salt.csv")  #get favor of genus
mycol8 <- c("#7fb961","#356d67","#42465c",
             "#76afda","#ffc556","#e8743c",
             "#b20000","#737373")

level.group <- c("light", "heavy")
m.salt <- vector("list", 2)
names(m.salt) <- level.group
for(i in 1:2){
  i1 <- which(design$Salt==level.group[i])
  x <- M.merge[[i1[1]]]
  for(j1 in 2:length(i1)){
    x <- x+ M.merge[[i1[j1]]]
  }
  x <- x/18
  #x[x>0] <- 1
  
  
  i1 <- match( colnames(m.OTU), rownames(x))  #change the order of OTU
  x <- x[i1,]    #currently the order of OTU is based on the relative abundance
  
  i2 <- match(info$name, colnames(x))
  x <- x[,i2] 
  m.salt[[i]] <- x
}

i.low <- 1:20

###### C P S N  #######
name.cyc <- c("Carbon", "Phosphorus", "Sulfur", "Nitrogen")

for(i in 1:4){
  ## weighted 
  tiff(paste0("Figs/Fig_bipartite_weighted_salt_avg_cyc_",name.cyc[i],"_1.tif"),width=10, height=4, units="in",
       res=300, compression = "lzw")
  name.gene.m <- colnames(m.salt[[1]])
  name.gene.m <- unlist(lapply(strsplit(name.gene.m,"_"),"[[",2))
  
  gene.path1 <- gene.path[gene.path$group==name.cyc[i],]
  gene.path1 <- gene.path1[order(gene.path1$path),]
  gene.path1$path <- as.integer(gene.path1$path)
  
  i.high <- match(gene.path1$gene_name, name.gene.m)
  
  size.low1 <- colMeans(m.OTU[design$Salt=="light", i.low])
  size.low1 <- 1+size.low1 * 20
  names(size.low1) <- 1:20
  size.high1 <- colMeans(m.gene[design$Salt=="light",])
  size.high1 <- 1 + size.high1 * 20
  names(size.high1) <- name.gene.m
  
  m1 <- m.salt[[1]]  #light
  colnames(m1) <- name.gene.m
  m1 <- m1[i.low,i.high]
  rownames(m1) <- 1:20
  plotweb(m1, col.high = mycol8[gene.path1$path],
          col.low="darkgreen",  high.lab.dis=0.06,
          col.interaction = "black", bor.col.interaction = NA, text.rot =0, empty=F,
          method="normal",labsize = 1.3, text.low.col = "black", low.lab.dis = 0.05)  #method="normal"
  
  size.low2 <- colMeans(m.OTU[design$Salt=="heavy", i.low])
  size.low2 <- 1+size.low2 * 20
  names(size.low2) <- 1:20
  size.high2 <- colMeans(m.gene[design$Salt=="heavy",])
  size.high2 <- 1 + size.high2 * 20
  names(size.high2) <- name.gene.m
  dev.off()
  
  tiff(paste0("Figs/Fig_bipartite_weighted_salt_avg_cyc_",name.cyc[i],"_2.tif"),width=10, height=4, units="in",
       res=300, compression = "lzw")
  m2 <- m.salt[[2]]  #heavy
  colnames(m2) <- name.gene.m
  m2 <- m2[i.low,i.high]
  rownames(m2) <- 1:20
  plotweb(m2, col.high = mycol8[gene.path1$path],
          col.low="darkgreen", high.lab.dis=0.06,
          col.interaction = "black", bor.col.interaction = NA, text.rot =0, empty=F,
          method="normal",labsize = 1.3, text.low.col = "black", low.lab.dis = 0.05)  #method="normal"
 
  #legend("right", fill=mycol8, legend=1:8)
  dev.off()
  
}


