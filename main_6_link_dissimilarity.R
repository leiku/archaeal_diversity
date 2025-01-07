
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(nlme)
library(vegan)

colors <- c(low="deepskyblue", high="brown1")

### abundance of OTU  ##########
d.OTU <- readRDS("Data_OTU_Archaea.RDS")
m.OTU <- d.OTU$OTU[, -(1:4)]
m.OTU <- t(m.OTU)
colnames(m.OTU) <- d.OTU$OTU$Taxon.name

m.OTU.b <- m.OTU / rowSums(m.OTU)
tmp <- colMeans(m.OTU.b)
m.OTU.b <- m.OTU.b[, order(tmp, decreasing = T)]

name.genus20 <- colnames(m.OTU.b)[1:20]

design <- d.OTU$design
design$Salt <- as.character(design$Salt)
#re-assign index: H for heavy, L for light; A for BH, B for HX, C for YQ
design$Site[design$Site=="BH"] <- "A"
design$Site[design$Site=="HX"] <- "B"
design$Site[design$Site=="YQ"] <- "C"
design$Treat <- paste0(rep(rep(c("H","L"),each=6),3), rep(c("A","B","C"), each=12),
                       design$Replicate)
design$Salt <- factor(design$Salt, levels=c("light","heavy"),labels=c("low","high"))



M.merge <- readRDS("Matrix_OTU_gene_Archaea.RDS")
m <- M.merge[[1]]
ii <- match(rownames(m), name.genus20)
tmp <- data.frame(index=1:nrow(m), position=ii)
tmp <- tmp[!is.na(tmp$position),]
i.position <- tmp[order(tmp$position),]

for(i in 1:length(M.merge)){
  m <- M.merge[[i]]
  m <- m[i.position$index,]
  M.merge[[i]] <- m
}

tmp <- dim(M.merge[[1]])



name.link <- "na"
for(i in 1:tmp[2]){
  name.gene <- colnames(M.merge[[1]])[i]
  name.link <- c(name.link, paste0(1:20,"-", name.gene))
}
name.link <- name.link[-1]

edge.list <- matrix(0, tmp[1] * tmp[2], length(M.merge))
for(i in 1:length(M.merge)){
  m <- M.merge[[i]]
  edge.list[,i] <- as.vector(m)
}
edge.list <- t(edge.list) 

i.good <- which(colSums(edge.list)>0)
edge.list <- edge.list[,i.good]
name.link <- name.link[i.good]

ii <- order(colSums(edge.list), decreasing=T)
edge.list <- edge.list[,ii]
name.link <- name.link[ii]

tmp <- strsplit(name.link, "-")

colnames(edge.list) <- name.link
rownames(edge.list) <- design$Treat



########## beta diversity of links #########

dist.edge  <- vegdist(edge.list, method="bray", binary=F)

mod.edge <- betadisper(dist.edge, group=design$Salt,
                       type="centroid", bias.adjust = F)

permutest(mod.edge)
TukeyHSD(mod.edge)
anosim(edge.list, grouping=design$Salt, distance="bray")

tiff("Figs/Fig_links_beta_disper.tif",width=5, height=7, units="in", res=300, compression = "lzw")
op <- par(mfrow=c(2,1), mar=c(4,4,2.5,1), mgp=c(2,1,0))

prop1 <- round(mod.edge$eig[1]/sum(mod.edge$eig), 3)*100
prop2 <- round(mod.edge$eig[2]/sum(mod.edge$eig), 3)*100
plot(mod.edge, main=NA, sub=NA, ellipse=F, hull =T, segments=T, 
     label=F, pch=rep(20,5), col=colors,cex=1.3, 
     xlab=paste0("PC1 (", prop1,"%)"), ylab=paste0("PC2 (", prop2,"%)"))
legend("topleft", legend=c("low", "high"), col=colors, pch=20)
mtext(LETTERS[1], side=3, line=0.8, adj=-0.15)

boxplot(mod.edge$distances~mod.edge$group,  col=colors, ylim=c(0.2, 0.8),
        xlab=NA, ylab="Distance to centroid")
text(1.5, 0.8, "***")
segments(x0=1.3, y0=0.78, x1=1.7, y1=0.78, col="black")
mtext(LETTERS[2], side=3, line=0.8, adj=-0.15)

par(op)
dev.off()



#########  show the dissimilarity of links  ##########


library (reshape2)#

### light #####
edge.light <- edge.list[6+c(1:6,13:18,25:30),]
edge.light <- edge.light[, order(colSums(edge.light), decreasing = T)]
edge.light1 <- as.data.frame(edge.light[,1:50])
edge.light1 <- edge.light1[,50:1]
e.light <- edge.light1 %>% mutate(Treat=row.names(.)) %>% melt()
e.light <- e.light[e.light$value>0,]
tmp <- strsplit(levels(e.light$variable), "-")

tmp1 <- unlist(lapply(tmp, function(x) paste0(x[1], substr(x[2],2,nchar(x[2])))))
e.light$variable <- factor(e.light$variable, labels=tmp1)

p0<-ggplot(e.light,aes(x=Treat,y=variable)) #热图绘制
p1 <- p0+scale_color_gradientn(colours = c('#6699CC','#FFFF99','#FF822E','#CC3333'))+
  theme_bw()+
  geom_point(aes( size=value, color=value))+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =90,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="gray", size=1, linetype="solid"))+
  geom_vline(xintercept=c(18.5),size=.8) +
  theme(axis.text = element_text(size = 12))


##### heavy ######
edge.heavy <- edge.list[c(1:6,13:18,25:30),]
edge.heavy <- edge.heavy[, order(colSums(edge.heavy), decreasing = T)]
edge.heavy1 <- as.data.frame(edge.heavy[,1:50])
edge.heavy1 <- edge.heavy1[,50:1]
e.heavy <- edge.heavy1 %>% mutate(Treat=row.names(.)) %>% melt()
e.heavy <- e.heavy[e.heavy$value>0,]

tmp <- strsplit(levels(e.heavy$variable), "-")
tmp1 <- unlist(lapply(tmp, function(x) paste0(x[1], substr(x[2],2,nchar(x[2])))))
e.heavy$variable <- factor(e.heavy$variable, labels=tmp1)

p0<-ggplot(e.heavy,aes(x=Treat,y=variable)) #heatmap
p2 <- p0+scale_color_gradientn(colours = c('#6699CC','#FFFF99','#FF822E','#CC3333'))+
  theme_bw()+
  geom_point(aes( size=value, color=value))+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =90,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="gray", size=1, linetype="solid"))+
  geom_vline(xintercept=c(18.5),size=.8) +
  theme(axis.text = element_text(size = 12))


tiff("Figs/Fig_link_diversity.tif",width=9, height=10, units="in",
     res=300, compression = "lzw")
p1 + p2 + 
  plot_annotation(tag_levels = 'A')+
  plot_layout(ncol=2)
dev.off()


