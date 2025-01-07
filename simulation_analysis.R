
library(vegan)
library(ggplot2)
library(patchwork)
library(ggsci)
library(nlme)
library(emmeans)

colors <- c("deepskyblue", "coral", "darkmagenta")

d <- readRDS("Data_simulation.RDS")
Treatment <- d$Treatment
B <- d$Biomass
I.SP <- d$I.SP
comp.list <- d$comp.list

######## alpha diversity ########
res.alpha <- Treatment
res.alpha$Shannon=0
res.alpha$Simpson=0
res.alpha$Richness=0
res.alpha$Evenness=0
res.alpha$Dominance=0

for(i in 1:nrow(B)){
  m <- B[i,]
  m <- m[m>0]
  res.alpha$Shannon[i] <- diversity(m, index="shannon", MARGIN = 1)
  res.alpha$Simpson[i] <- diversity(m, index="invsimpson", MARGIN = 1)
  res.alpha$Richness[i] <- specnumber(m)
  p <- m/sum(m)
  res.alpha$Dominance[i] <- sum(p^2)
}
res.alpha$Evenness <- res.alpha$Shannon / log(res.alpha$Richness)

res.alpha$treat <- factor(res.alpha$treat, levels=c("control","RS","RC"))


myplot.alpha <- function(data=res.alpha, x="treat", y="Richness", cols=colors){
  ymax <- tapply(data[,y], data[,x], max)
  ymin <- min(data[,y])
  ysig <- ymax+0.1*(max(ymax) - min(ymin))
  g0 <- ggplot(data, aes_string(x=x, y=y, fill=x)) + 
    labs(x="", y=y) + 
    ylim(ymin, max(ymax)+0.12*(max(ymax) - ymin))+
    scale_fill_manual(values=cols) +
    geom_boxplot(alpha=0.3, outlier.shape = NA) + 
    geom_jitter(aes_string(fill=x), width=0.2, shape=21, size=2) + 
    annotate("text", x=1:3, y=ysig, label=c("a", "b", "c"))
  return(g0)
}


m1 <- lme(Richness~treat, random=~1|replicate, data=res.alpha)
anova(m1)
emmeans(m1, list(pairwise ~ treat), adjust = "tukey")
m2 <- lme(Shannon~treat, random=~1|replicate, data=res.alpha)
anova(m2)
emmeans(m2, list(pairwise ~ treat), adjust = "tukey")
m3 <- lme(Simpson~treat, random=~1|replicate, data=res.alpha)
anova(m3)
emmeans(m3, list(pairwise ~ treat), adjust = "tukey")
m4 <- lme(Dominance~treat, random=~1|replicate, data=res.alpha)
anova(m4)
emmeans(m4, list(pairwise ~ treat), adjust = "tukey")

g1 <- myplot.alpha(data=res.alpha, x="treat", y="Richness", cols=colors)
g2 <- myplot.alpha(data=res.alpha, x="treat", y="Shannon", cols=colors)
g3 <- myplot.alpha(data=res.alpha, x="treat", y="Simpson", cols=colors)
g4 <- myplot.alpha(data=res.alpha, x="treat", y="Dominance", cols=colors)



tiff("Figs/Fig_simulation_alpha.tif",width=8, height=6, units="in", res=300, compression = "lzw")

g1 + g2 +  g3 + g4 +
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = 'collect', ncol=2)&
  theme_bw()+
  theme(legend.position="none",
        #legend.title = element_blank(),
        panel.grid=element_blank())

dev.off()



#### pcoa ####

m <- B
m <- m[, colSums(m)>0]
dist  <- vegdist(m, method="bray", binary=F)

treat <- Treatment$treat
treat <- factor(treat, levels=c("control", "RS", "RC"))
mod <- betadisper(dist, group=treat,
                  type="centroid", bias.adjust = F)

permutest(mod)
TukeyHSD(mod)
anosim(m, grouping=treat, distance="bray")

beta <- data.frame(distances=mod$distances, 
                   group=mod$group,
                   ID=rep(1:18, 3))
m <- lme(distances~group, random=~1|ID, data=beta)
anova(m)
emmeans(m, list(pairwise ~ group), adjust = "tukey")

tiff("Figs/Fig_simulation_beta_disper.tif",width=8, height=3.5, units="in", res=300, compression = "lzw")
op <- par(mfrow=c(1,2), mar=c(4,4,2.5,1), mgp=c(2,1,0))

prop1 <- round(mod$eig[1]/sum(mod$eig), 3)*100
prop2 <- round(mod$eig[2]/sum(mod$eig), 3)*100
plot(mod, main=NA, sub=NA, ellipse=F, hull =T, segments=T, 
     label=F, pch=rep(20,3), col=colors,cex=1.3, 
     xlab=paste0("PC1 (", prop1,"%)"), ylab=paste0("PC2 (", prop2,"%)"))
legend("bottomright", legend=c("control", "RS",  "RC"), col=colors, pch=20)
mtext(LETTERS[5], side=3, line=0.8, adj=-0.15)
boxplot(mod$distances~mod$group,  col=colors, ylim=c(0, 0.4),
        xlab=NA, ylab="Distance to centroid")
text(x=1:3, y=c(0.28, 0.2, 0.39), labels = c("a","a","b"))
mtext(LETTERS[6], side=3, line=0.8, adj=-0.15)

par(op)
dev.off()


##### dissimilarity of links ########
#comp.list中存放的是i.s指示的200个物种间的竞争，需先转回到1000的物种池中，再删掉灭绝的物种，再计算
i.exist <- which(colSums(B) > 0)
nn <- length(i.exist)
edge.list <- matrix(0, nn*nn, nrow(B))

I.SP <- rbind(I.SP, I.SP, I.SP)
for(i in 1:nrow(B)){
  i.sp <- I.SP[i,]   #index of original species
  comp <- comp.list[[i]]
  
  ii <- match(i.sp, i.exist) #position of sp. in existing species
  ii <- ii[!is.na(ii)]
  comp <- comp[ii,ii]
  i1 <- matrix(rep(ii, length(ii)), length(ii))  #index of row
  j1 <- matrix(rep(ii, each=length(ii)), length(ii))  #index of column
  i1 <- as.vector(i1)
  j1 <- as.vector(j1)
  comp1 <- as.vector(comp)
  
  tmp <- matrix(0, nn, nn)
  for(j in 1:length(comp1)){
    tmp[i1[j],j1[j]] <- comp1[j]
  }
  diag(tmp) <- 0
  edge.list[,i] <- as.vector(tmp)
}
edge.list <- edge.list[rowSums(edge.list)>0,]


dist.edge  <- vegdist(t(edge.list), method="bray", binary=F)

treat <- Treatment$treat
treat <- factor(treat, levels=c("control", "RS", "RC"))
mod.edge <- betadisper(dist.edge, group=treat,
                  type="centroid", bias.adjust = F)

permutest(mod.edge)
TukeyHSD(mod.edge)
anosim(t(edge.list), grouping=treat, distance="bray")

tiff("Figs/Fig_simulation_links_beta_disper.tif",width=5, height=7, units="in", res=300, compression = "lzw")
op <- par(mfrow=c(2,1), mar=c(4,4,2.5,1), mgp=c(2,1,0))

prop1 <- round(mod.edge$eig[1]/sum(mod.edge$eig), 3)*100
prop2 <- round(mod.edge$eig[2]/sum(mod.edge$eig), 3)*100
plot(mod.edge, main=NA, sub=NA, ellipse=F, hull =T, segments=T, 
     label=F, pch=rep(20,5), col=colors,cex=1.3, 
     xlab=paste0("PC1 (", prop1,"%)"), ylab=paste0("PC2 (", prop2,"%)"))
legend("bottomright", legend=c("control", "RS",  "RC"), col=colors, pch=20)
mtext(LETTERS[1], side=3, line=0.8, adj=-0.15)

boxplot(mod.edge$distances~mod.edge$group,  col=colors, ylim=c(0.45, 0.63),
        xlab=NA, ylab="Distance to centroid")
text(x=1:3, y=c(0.51, 0.54, 0.62), labels = c("a","b","c"))
mtext(LETTERS[2], side=3, line=0.8, adj=-0.15)

par(op)
dev.off()


