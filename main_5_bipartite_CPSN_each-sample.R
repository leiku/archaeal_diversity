
library(bipartite)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(nlme)

colors <- c(low="deepskyblue", high="brown1")

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

design$Salt <- factor(design$Salt, levels=c("light","heavy"), labels=c("low","high"))

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


#### statistic ####

library(readxl)
gene.path <- read_xlsx("Res_Gene_function_Archaea.xlsx", sheet=1)
gene.path <- gene.path[gene.path$path != "NA",]

M.merge <- readRDS("Matrix_OTU_gene_Archaea.RDS")
#res.fit <- read.csv("res_fitting_log-abundance-vs-salt.csv")  #get favor of genus


name.cyc <- c("Carbon", "Nitrogen", "Phosphorus", "Sulfur")
name.salt <- c("low", "high")
i.low <- 1:20
res.degree <- matrix(0, 1, 7)
colnames(res.degree) <- c("site","salt","cyclying","genus","abundance",
                          "degree","degree.weighted")
res.degree <- as.data.frame(res.degree)

res.connectance <- matrix(0, 36*4, 4)
colnames(res.connectance) <- c("site","salt","cyclying","C")
res.connectance <- as.data.frame(res.connectance)
k <- 1

for(ii in 1:4){
  
  for(i in 1:2){
    design1 <- design[design$Salt==name.salt[i],]
    M.merge1 <- M.merge[design$Salt==name.salt[i]]
    m.OTU1 <- m.OTU[design$Salt==name.salt[i], i.low]
    m.gene1 <- m.gene[design$Salt==name.salt[i],]
    
    for(j in 1:18){
      x <- M.merge1[[j]]
      i1 <- match( colnames(m.OTU), rownames(x))  #change the order of OTU
      x <- x[i1,]    #currently the order of OTU is based on the relative abundance
      
      i2 <- match(info$name, colnames(x))
      x <- x[,i2] 
      
      name.gene.m <- colnames(x)
      name.gene.m <- unlist(lapply(strsplit(name.gene.m,"_"),"[[",2))
      
      gene.path1 <- gene.path[gene.path$group==name.cyc[ii],]
      gene.path1 <- gene.path1[order(gene.path1$path),]
      gene.path1$path <- as.integer(gene.path1$path)
      
      i.high <- match(gene.path1$gene_name, name.gene.m)
      
      colnames(x) <- name.gene.m
      x1 <- x[i.low,i.high]
      
      x1.binary <- x1
      x1.binary[x1.binary>0] <- 1
      
      
      
      res1 <- matrix(0, nrow(x1), 7)
      colnames(res1) <- c("site","salt","cyclying","genus","abundance",
                                "degree","degree.weighted")
      res1 <- as.data.frame(res1)
      res1$site <- design1$Treat[j]
      res1$salt <- name.salt[i]
      res1$cyclying <- name.cyc[ii]
      res1$genus <- rownames(x1)
      res1$abundance <- m.OTU1[j,]
      res1$degree <- rowSums(x1.binary)
      res1$degree.weighted <- rowSums(x1)
      
     res.degree <- rbind(res.degree, res1)
     
     res.connectance$site[k] <- design1$Treat[j]
     res.connectance$salt[k] <- name.salt[i]
     res.connectance$cyclying[k] <- name.cyc[ii]
     res.connectance$C[k] <- sum(x1.binary) /nrow(x1) / ncol(x1)
     k <- k+1
    }
  }
}
res.degree <- res.degree[-1,]

res.avg <- res.degree %>% group_by(salt, cyclying, genus) %>%
  summarise(abundance.avg=mean(abundance), abundance.sd=sd(abundance),
            degree.avg=mean(degree), degree.sd=sd(degree),
            degree.weighted.avg=mean(degree.weighted), 
            degree.weighted.sd=sd(degree.weighted))
res.avg$salt <- factor(res.avg$salt, levels=c("low","high"))

for(i in 1:2){
  for(j in 1:4){
    z <- summary(lm(log10(degree.avg+1)~log10(abundance.avg), 
                    data=res.avg[res.avg$cyclying==name.cyc[j] & res.avg$salt==name.salt[i],]))
    print(z$coef[2,4])
  }
}
anova(lm(log10(degree.avg+1)~log10(abundance.avg) * salt, 
       data=res.avg[res.avg$cyclying==name.cyc[1], ]))


myfitting <- function(data=res.avg, cols=colors, legends=F){
  data$x <- log10(data$abundance.avg)
  data$y <- log10(data$degree.avg + 1)
  
  z1 <- lm(y~x, data=data[data$salt=="low",])
  z2 <- lm(y~x, data=data[data$salt=="high",])
  z11 <- summary(z1)
  z21 <- summary(z2)
  text1 <- paste("slope == ",round(coef(z1)[2], 3),"~~R^{2}==", round(z11$r.squared,3))
  text2 <- paste("slope == ",round(coef(z2)[2], 3),"~~R^{2}==", round(z21$r.squared,3))
  
  ymax <- max(data$y)
  ymin <- min(data$y)
  xmin <- min(data$x)
  xmax <- max(data$x)
  g0 <- ggplot(data, aes(x=x, y=y, col=salt)) + 
    geom_point() +
    geom_smooth(method="lm") + 
    labs(x="Relative abundance (log)", y="Number of genes (log)") + 
    ylim(ymin, ymax+0.05*(ymax - ymin))+
    scale_color_manual(values=colors) +
    annotate("text", x=xmin+0*(xmax-xmin), y=ymax, 
             label=text1, parse=T, hjust=0, col=colors[1]) +
    annotate("text", x=xmin+0*(xmax-xmin), y=ymax-0.1*(ymax-ymin), 
             label=text2, parse=T, hjust=0, col=colors[2]) +
    theme_bw()
  if(legends){g0 <- g0 + theme(legend.position=c(0.85,0.1),
                               legend.title = element_blank(),
                                panel.grid=element_blank())}else{
                g0 <- g0 + theme(legend.position="none",
                                panel.grid=element_blank())
                                }
  return(g0)
}


tiff("Figs/Fig_bipartite_degree.tif",width=9, height=8, units="in",
     res=300, compression = "lzw")
g1 <- myfitting(data=res.avg[res.avg$cyclying=="Carbon",], legends=T)
g2 <- myfitting(data=res.avg[res.avg$cyclying=="Nitrogen",])
g3 <- myfitting(data=res.avg[res.avg$cyclying=="Phosphorus",])
g4 <- myfitting(data=res.avg[res.avg$cyclying=="Sulfur",])
g1 + g2 + g3 + g4 + 
  plot_annotation(tag_levels = 'A')+
  plot_layout(ncol=2)
dev.off()



tiff("Figs/Fig_bipartite_weighted_degree.tif",width=9, height=8, units="in",
     res=300, compression = "lzw")
ggplot(data=res.avg, aes(x=log10(abundance.avg), y=log10(degree.weighted.avg), col=salt)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(.~cyclying)+
  labs(x="Relative abundance (Log10)", y="Sum of gene abundance (Log10)")+
  theme_bw()+
  scale_color_manual(values=colors)
dev.off()



### connectance

res.connectance$salt <- factor(res.connectance$salt, levels=c("low","high"))
t.test(C~salt, data=res.connectance[res.connectance$cyclying=="Carbon",])
t.test(C~salt, data=res.connectance[res.connectance$cyclying=="Nitrogen",])
t.test(C~salt, data=res.connectance[res.connectance$cyclying=="Phosphorus",])
t.test(C~salt, data=res.connectance[res.connectance$cyclying=="Sulfur",])

tiff("Figs/Fig_bipartite_connectance.tif",width=5, height=4, units="in",
     res=300, compression = "lzw")
ggplot(data=res.connectance, aes(x=cyclying, y=C, fill=salt)) +
  geom_boxplot() +
  labs(x="", y="Connectance")+
  annotate("segment", x=c(0.8,1.8,2.8,3.8), xend=c(1.2,2.2,3.2,4.2),
           y=c(0.17,0.25, 0.27, 0.25), yend=c(0.17,0.25, 0.27, 0.25)) +
  annotate("text", x=1:4, y=c(0.17,0.25, 0.27, 0.25) + 0.01,
           label=c("***","***","***","***")) + 
  theme_bw()+
  scale_fill_manual(values=colors)+
  theme(legend.position = c(0.1,0.9),
        legend.title = element_blank(),
        panel.grid=element_blank())
dev.off()



#### bipartite ####

library(readxl)
gene.path <- read_xlsx("Res_Gene_function_Archaea.xlsx", sheet=1)
gene.path <- gene.path[gene.path$path != "NA",]

M.merge <- readRDS("Matrix_OTU_gene_Archaea.RDS")
#res.fit <- read.csv("res_fitting_log-abundance-vs-salt.csv")  #get favor of genus
mycol8 <- c("#7fb961","#356d67","#42465c",
            "#76afda","#ffc556","#e8743c",
            "#b20000","#737373")


name.cyc <- c("Carbon", "Phosphorus", "Sulfur", "Nitrogen")
name.salt <- c("light", "heavy")
i.low <- 1:20

for(ii in 1:4){
  pdf(paste0("Figs/Fig_bipartite_weighted_each_",name.cyc[ii],".pdf"))
  op <- par(mfrow=c(3,2), mar=c(4,4,1,1), mgp=c(2,1,0))
  
  for(i in 1:2){
    design1 <- design[design$Salt==name.salt[i],]
    M.merge1 <- M.merge[design$Salt==name.salt[i]]
    m.OTU1 <- m.OTU[design$Salt==name.salt[i], i.low]
    m.gene1 <- m.gene[design$Salt==name.salt[i],]
    
    for(j in 1:18){
      x <- M.merge1[[j]]
      i1 <- match( colnames(m.OTU), rownames(x))  #change the order of OTU
      x <- x[i1,]    #currently the order of OTU is based on the relative abundance
      
      i2 <- match(info$name, colnames(x))
      x <- x[,i2] 
      
      #x[x>0] <- 1
      
      name.gene.m <- colnames(x)
      name.gene.m <- unlist(lapply(strsplit(name.gene.m,"_"),"[[",2))
      
      gene.path1 <- gene.path[gene.path$group==name.cyc[ii],]
      gene.path1 <- gene.path1[order(gene.path1$path),]
      gene.path1$path <- as.integer(gene.path1$path)
      
      i.high <- match(gene.path1$gene_name, name.gene.m)
      
      size.low1 <- m.OTU1[j,]
      size.low1 <- 1+size.low1 * 20
      names(size.low1) <- 1:20
      size.high1 <- m.gene1[j,]
      size.high1 <- 1 + size.high1 * 20
      names(size.high1) <- name.gene.m
      
      colnames(x) <- name.gene.m
      x1 <- x[i.low,i.high]
      rownames(x1) <- 1:20
      plotweb(x1, col.high = mycol8[gene.path1$path],
              col.low="darkgreen", 
              col.interaction = "black", bor.col.interaction = NA, text.rot =0, empty=F,
              low.abun =  size.low1, high.abun = size.high1[i.high],
              low.abun.col = "darkgray", high.abun.col = "darkgray",
              method="normal",labsize = 1.3, text.low.col = "black", low.lab.dis = 0.1)  #method="normal"
      mtext(design1$Treat[j], side=3, line=0)
    }
  }
  
  par(op)
  dev.off()
}



