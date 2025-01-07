 
library(vegan)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggpubr)
library(bipartite)


d.OTU <- readRDS("Data_OTU_Archaea.RDS")
m.OTU <- d.OTU$OTU[, -(1:4)]
m.OTU <- t(m.OTU)
colnames(m.OTU) <- d.OTU$OTU$Taxon.name

design <- d.OTU$design
#re-assign index: H for heavy, L for light; A for BH, B for HX, C for YQ
design$Treat <- paste0(rep(rep(c("H","L"),each=6),3), rep(c("A","B","C"), each=12),
                      design$Replicate)

rownames(m.OTU) <- design$Treat
m.OTU <- m.OTU[order(design$Salt),]
design <- design[order(design$Salt),]


##### Archaea #####
m.OTU.b <- m.OTU / rowSums(m.OTU)
tmp <- colMeans(m.OTU.b)
m.OTU.b <- m.OTU.b[, order(tmp, decreasing = T)]

data_m <- melt(m.OTU.b[, 1:20]) 
colnames(data_m) <- c("treat", "taxon", "proportion")

ggplot(data_m,aes(x=taxon,y=treat)) +
  geom_tile(aes(fill=proportion)) +
  labs(x="", y="")+
  theme_bw()+
  theme(axis.text.x = element_text(
    angle = 45,hjust = 1,vjust = 1))  +
  gradient_fill(c("white","red"))



####### cluster and stacked

dis_bray <- vegdist(m.OTU.b, method="bray")  #Jaccard
tree <- hclust(dis_bray)   

tree$order <- rev(tree$order)

group <- data.frame(names=design$Treat, group=design$Salt, group_num=as.numeric(design$Salt))
grp <- group[2]
rownames(grp) <- group$names
group_col <- c('deepskyblue', 'brown1')
names(group_col) <- c('light', 'heavy')
group_name <- c('light', 'heavy')


###### plot ####

tiff("Figs/Fig_S1_stacked_20.tif", height=6, width=11, unit="in", res=300, compression = "lzw")

#
layout(t(c(1, 2, 2, 2, 3)))
par(mar = c(5, 2, 5, 0))

plot(0, type = 'n', xaxt = 'n', yaxt = 'n', frame.plot = FALSE, xlab = '', ylab = '',
     xlim = c(-max(tree$height), 0), ylim = c(0, length(tree$order)))
legend('topleft', legend = c("low","high"), pch = 15, col = group_col, bty = 'n', cex = 1.5)

#clustering
treeline <- function(pos1, pos2, height, col1, col2) {
  meanpos = (pos1[1] + pos2[1]) / 2
  segments(y0 = pos1[1] - 0.4, x0 = -pos1[2], y1 = pos1[1] - 0.4, x1 = -height,  col = col1,lwd = 2)
  segments(y0 = pos1[1] - 0.4, x0 = -height,  y1 = meanpos - 0.4, x1 = -height,  col = col1,lwd = 2)
  segments(y0 = meanpos - 0.4, x0 = -height,  y1 = pos2[1] - 0.4, x1 = -height,  col = col2,lwd = 2)
  segments(y0 = pos2[1] - 0.4, x0 = -height,  y1 = pos2[1] - 0.4, x1 = -pos2[2], col = col2,lwd = 2)
}

meanpos = matrix(rep(0, 2 * length(tree$order)), ncol = 2)
meancol = rep(0, length(tree$order))
for (step in 1:nrow(tree$merge)) {
  if(tree$merge[step, 1] < 0){
    pos1 <- c(which(tree$order == -tree$merge[step, 1]), 0)
    col1 <- group_col[as.character(grp[tree$labels[-tree$merge[step, 1]],1])]
  } else {
    pos1 <- meanpos[tree$merge[step, 1], ]
    col1 <- meancol[tree$merge[step, 1]]
  }
  if (tree$merge[step, 2] < 0) {
    pos2 <- c(which(tree$order == -tree$merge[step, 2]), 0)
    col2 <- group_col[as.character(grp[tree$labels[-tree$merge[step, 2]],1])]
  } else {
    pos2 <- meanpos[tree$merge[step, 2], ]
    col2 <- meancol[tree$merge[step, 2]]
  }
  height <- tree$height[step]
  treeline(pos1, pos2, height, col1, col2)
  meanpos[step, ] <- c((pos1[1] + pos2[1]) / 2, height)
  if (col1 == col2) meancol[step] <- col1 else meancol[step] <- 'grey'
}

##stacked bar
dat <- t(m.OTU.b)
dat <- dat[ ,tree$order]
dat <- dat[1:20,]
dat <- rbind(dat, 1-colSums(dat))
rownames(dat)[21] <- "Others"

phylum_color <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray')
mycol20 <- c("#b0d45d","#7fb961","#4c9568","#356d67","#42465c","#5066a1",
             "#76afda","#fddbc8","#ffc556",
             "#e8743c","#f06152","#b20000","#eb998b","#a14462","#cca69c",
             "#7d4444","#562e3c","#35212e","#737373","black","#bdbdbd")

names(mycol20) <- rownames(dat)

par(mar = c(5, 2, 5, 0))

tt <- barplot(as.matrix(dat), col = mycol20, space = 0.4, width = 0.7, cex.axis = 1, horiz = TRUE, cex.lab = 1.2, 
              xlab = 'Relative Abundance', yaxt = 'n', las = 1, ylim = c(0, ncol(dat)))

#mtext('Top 10 phylums', side = 3, line = 1, cex = 1)
text(x = -0.03, y = tt, labels = colnames(dat), col = group_col[group[tree$order, 2]], xpd = TRUE)

par(mar=c(5,1,5,0))
plot(0, type="n", xaxt="n", yaxt="n", bty="n",xlab="",ylab="")
legend("left", pch=15, col=mycol20, legend=names(mycol20), bty="n", cex=1.2)

dev.off()


