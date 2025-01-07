
#### get population matrix #####  
library(vegan)
library(ggplot2)
library(patchwork)
library(nlme)
library(ade4)
library(ape)   #for pcoa

colors <- c(low="deepskyblue", high="brown1")


d.OTU <- readRDS("Data_OTU_Archaea_all.RDS")
m.OTU <- d.OTU$OTU[, -1]
m.OTU <- t(m.OTU)
colnames(m.OTU) <- d.OTU$taxonomy$lowest_taxon

n.sample <- apply(m.OTU, 2, function(x) sum(x>0))

design <- d.OTU$design
design$Salt <- factor(design$Salt, levels=c("light","heavy"), labels=c("low","high"))

m <- m.OTU / rowSums(m.OTU)

##### alpha diversity ######
res.alpha <- data.frame(Site=design$Site,
                        Salt=design$Salt,
                        Replicate= design$Replicate,
                        Shannon=diversity(m, index="shannon"),
                        Simpson=diversity(m, index="invsimpson"),
                        Richness=specnumber(m),
                        Dominance=rowSums(m^2))  #Simpson dominance
res.alpha$Evenness <- res.alpha$Shannon / log(res.alpha$Richness)
res.alpha$Richness <- log10(res.alpha$Richness)

myplot.alpha <- function(data=res.alpha, x="Salt", y="Richness", ylab="Richness", cols=colors){
  ymax <- tapply(data[,y], data[,x], max)
  ymin <- min(data[,y])
  ystar <- max(ymax)+0.04*(max(ymax) - min(ymin))
  yline <- max(ymax)+0.03*(max(ymax) - min(ymin))
  g0 <- ggplot(data, aes_string(x=x, y=y, fill=x)) + 
    labs(x="", y=ylab) + 
    ylim(ymin, max(ymax)+0.05*(max(ymax) - ymin))+
    scale_fill_manual(values=cols) +
    geom_boxplot(alpha=0.3, outlier.shape = NA) + 
    geom_jitter(aes_string(fill=x), width=0.2, shape=21, size=2) + 
    annotate("text", x=1.5, y=ystar, label=c("***")) +
    annotate("segment", x=1.2, xend=1.8, y=yline, yend=yline)
  return(g0)
}

anova(lme(Richness~Salt, random=~1|Site, data=res.alpha))
g1 <- myplot.alpha(data=res.alpha, x="Salt", y="Richness", ylab="Richness (log)", cols=colors)
anova(lme(Shannon~Salt, random=~1|Site, data=res.alpha))
g2 <- myplot.alpha(data=res.alpha, x="Salt", y="Shannon", ylab="Shannon", cols=colors)
anova(lme(Simpson~Salt, random=~1|Site, data=res.alpha))
g3 <- myplot.alpha(data=res.alpha, x="Salt", y="Simpson", ylab="Simpson", cols=colors)
anova(lme(Dominance~Salt, random=~1|Site, data=res.alpha))
g4 <- myplot.alpha(data=res.alpha, x="Salt", y="Dominance", ylab="Dominance", cols=colors)



tiff("Figs/Fig1_alpha_OTU.tif",width=6, height=5, units="in",
     res=300, compression = "lzw")
g1 + g2 + g3 + g4 + 
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = 'collect', ncol=2)&
  theme_bw()+
  theme(legend.position="none",
        legend.title = element_blank(),
        panel.grid=element_blank()) 
dev.off()



#### pcoa ####

dist1  <- vegdist(m, method="bray", binary=F)
mod1 <- betadisper(dist1, group=design$Salt,
                     type="centroid", bias.adjust = F)

permutest(mod1)
TukeyHSD(mod1)
tmp <- data.frame(distances=mod1$distances, group=mod1$group, Site = design$Site)
anova(lme(distances~group, random=~1|Site, data=tmp))
anosim(m, grouping=design$Salt, distance="bray")

tiff("Figs/Fig2_beta_disper_OTU.tif",width=6, height=2.8, units="in",
     res=300, compression = "lzw")
op <- par(mfrow=c(1,2), mar=c(4,3,1,1), mgp=c(2,1,0))

prop1 <- round(mod1$eig[1]/sum(mod1$eig), 3)*100
prop2 <- round(mod1$eig[2]/sum(mod1$eig), 3)*100
plot(mod1, main=NA, sub=NA, ellipse=F, hull =T, segments=T, 
       label=F, pch=c(15,17), col=colors,cex=1.3, 
       xlab=paste0("PC1 (", prop1,"%)"), ylab=paste0("PC2 (", prop2,"%)"))
legend("bottomright", legend=c("low","high"), col=colors, pch=c(15,17))
text(0, 0.18, "***")
segments(x0=-0.1, y0=0.16, x1=0.1, y1=0.16, col="black")
mtext("E", side=3, line=0, adj=-0.25)
boxplot(mod1$distances~mod1$group, ylim=c(0,0.75), col=alpha(colors, 0.5),
          xlab=NA, ylab="Distance to centroid")
text(1.5, 0.7, "***")
segments(x0=1.3, y0=0.67, x1=1.7, y1=0.67, col="black")
mtext("F", side=3, line=0, adj=-0.25)

par(op)
dev.off()

