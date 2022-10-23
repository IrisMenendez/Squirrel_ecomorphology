#Project: GM p4 squirrels
#R version: 3.5.1


library(ape)
library(geiger)
library(ggplot2)
library(geomorph)
library(RRPP)
library(motmot)
library(phytools)
library(StereoMorph)
library(mvMORPH)




setwd("/Users/iris/Documents/GitHub/Squirrel_ecomorphology")
set.seed(15)

######READING DATA######
data <- read.table("TableS2_species_info.csv", sep=";", dec=".", header=TRUE)
tree <- read.tree("Sciuridae_tree_Menendez_etal_2017.tree")
tree <- ladderize(tree, right=F)

# subset of data only with species shared between phy and data
rownames(data) <- data$Sp_red
matchdata.phy <- treedata(tree, data[,13:140], sort=T)
phy.inf <- matchdata.phy$phy
mshapes.phy.inf <- matchdata.phy$data
data.phy <- droplevels(data[match(phy.inf$tip.label, data$Sp_red),][,1:12])
mCsize.inf <- data.phy$mean.Centroid.size
names(mCsize.inf)<- data.phy$Sp_red
diet.phy.inf <- data.phy$Diet..this.study.
fundiet.phy.inf <- data.phy$Funtional_diet..this.study.

# subset of data for tree squirrels in phy and data
arb.phy <- drop.tip(tree, as.vector(data.phy[data.phy$Subfamily=="Xerinae",]$Sp_red))
matchdata.phyarb <- treedata(arb.phy, data[,13:140], sort=T)
arb.phy <- matchdata.phyarb$phy
mshapes.phy.inf.arb <- matchdata.phyarb$data
data.phy.inf.arb <- droplevels(data[match(arb.phy$tip.label, data$Sp_red),][,1:12])
mCsize.phy.inf.arb <- data.phy.inf.arb$mean.Centroid.size
names(mCsize.phy.inf.arb) <- data.phy.inf.arb$Sp_red


# subset of data for ground squirrels in phy and data
terr.phy <- drop.tip(tree, as.vector(data.phy[data.phy$Subfamily!="Xerinae",]$Sp_red))
matchdata.phyterr <- treedata(terr.phy, data[,13:140], sort=T)
terr.phy <- matchdata.phyterr$phy
mshapes.phy.inf.terr <- matchdata.phyterr$data
data.phy.inf.terr <- droplevels(data[match(terr.phy$tip.label, data$Sp_red),][,1:12])
mCsize.phy.inf.terr <- data.phy.inf.terr$mean.Centroid.size
names(mCsize.phy.inf.terr) <- data.phy.inf.terr$Sp_red




###### COLORS ######
levels(data$Diet..this.study.)
col.diet.inf <- c("#bb3e7b", "#69ccc7", "#4a7f5c", "#ff595e", "#7dd86f", "Pink", "#f4c44e", "#6340c6")
levels(data$Funtional_diet..this.study.)
col.dietCVA.inf <- c( "#a04a4a", "#add34f", "Pink", "Orange")




###### ANALYSES ######


###### 1. tooth size vs weight  v ######
#data
Weight <- data$mean.Weight
Weight <- log10(Weight^(1/3))
names(Weight) <- data$Sp_red
tooth_size <- log10(data$mean.Centroid.size)

size_dat <- as.data.frame(cbind(Weight, tooth_size))
size_dat <- size_dat[!is.na(Weight),]
diet_cor <- data$Diet..this.study.[!is.na(Weight)]

#correlation
cor.test(Weight, tooth_size, data=size_dat)

#Linear model
lm_size <- lm(tooth_size~Weight, data=size_dat)
summary(lm_size)
predlm <- predict(lm_size, interval = "confidence")

#Plot
weight_tooth_plot <- ggplot(data = size_dat, aes(x=Weight, y=tooth_size))+
  geom_point(col=col.diet.inf[diet_cor], cex=1.8)+
  geom_line(aes(y = predlm[,1]), size = 0.5)+
  geom_smooth(col="#039496", method="lm")+
  annotate("text", x=1.3, y=0.5, label=as.character(summary(lm_size)$r.squared))+
  #annotate("text", x=Weight, y=tooth_size,label=names(Weight), cex=1)+
  coord_fixed()+
  ylim(0.3, 1.5)+
  xlim(0.3, 1.5)+
  theme_classic()

pdf("weight_tooth_lm_plot_cuberoot2_text.pdf")
weight_tooth_plot
dev.off()




###### 2. PGLS  ######

##### + Csize ~ diet ######


size_lm <- procD.lm(data.phy$mean.Centroid.size ~ diet.phy.inf, phy.inf, iter = 10000, SS.type="II")
phy.sig_dietsize <- phylosig(phy.inf, size_lm$residuals, method="lambda", nsim=1000)

phy.lambda_dietsize <- rescale(phy.inf, "lambda",phy.sig_dietsize$lambda)
phylosig(phy.lambda_dietsize,size_lm$residuals, method="lambda", nsim=1000)

pgls_size <- procD.pgls(data.phy$mean.Centroid.size ~ diet.phy.inf, phy.lambda_dietsize, iter = 10000, SS.type="II")



size_lm_arb <- procD.lm(mCsize.phy.inf.arb ~ data.phy.inf.arb$Diet..this.study., arb.phy, iter = 10000, SS.type="II")
phy.sig_dietsize_arb <- phylosig(arb.phy, size_lm_arb$residuals, method="lambda", nsim=1000)
phy.sig_dietsize_arb$lambda
phy.lambda_dietsize_arb <- rescale(arb.phy, "lambda",phy.sig_dietsize_arb$lambda)
phylosig(phy.lambda_dietsize_arb,size_lm_arb$residuals, method="lambda", nsim=1000)

pgls_size.arb <- procD.pgls(mCsize.phy.inf.arb ~ data.phy.inf.arb$Diet..this.study., phy.lambda_dietsize_arb,iter = 10000, SS.type="II")


size_lm_terr <- procD.lm(mCsize.phy.inf.terr ~ data.phy.inf.terr$Diet..this.study., terr.phy, iter = 10000, SS.type="II")
phy.sig_dietsize_terr <- phylosig(terr.phy, size_lm_terr$residuals, method="lambda", nsim=1000)
phy.sig_dietsize_terr$lambda
phy.lambda_dietsize_terr <- rescale(terr.phy, "lambda",phy.sig_dietsize_terr$lambda)
phylosig(phy.lambda_dietsize_terr,size_lm_terr$residuals, method="lambda", nsim=1000)

pgls_size.terr <- procD.pgls(mCsize.phy.inf.terr ~ data.phy.inf.terr$Diet..this.study., phy.lambda_dietsize_terr,iter = 10000, SS.type="II")



write.table(anova(pgls_size)$table, "PGLS_dietsize.csv", dec=".", sep=";")
write.table(anova(pgls_size.arb)$table, "PGLS_dietsize.csv", dec=".", sep=";", append=T)
write.table(anova(pgls_size.terr)$table, "PGLS_dietsize.csv", dec=".", sep=";", append=T)


#pairwise comparisons 
pairw_dietsize <- pairwise(pgls_size, groups = diet.phy.inf)
pairw_dietsize_summary <- summary(pairw_dietsize, confidence = 0.95, test.type = "dist", stat.table =TRUE)

pair_adjust_size <- pairw_dietsize_summary$pairwise.tables$D
pair_adjust_size[upper.tri(pair_adjust_size)] <- p.adjust(pairw_dietsize_summary$pairwise.tables$P[upper.tri(pairw_dietsize_summary$pairwise.tables$P)], method="holm")

write.table(pair_adjust_size, "pairwise_dietsize_PGLS_ADJUSTED.csv", sep=";", dec = ".")





##### + Shape ~ diet*Csize ######

gdf <- geomorph.data.frame(Diet=diet.phy.inf,Shape=mshapes.phy.inf,LCS.diet=data.phy$mean.Centroid.size)
GLS.diet <- procD.lm(Shape~LCS.diet*Diet, data=gdf)
GLS.diet

physignal(phy=phy.inf,GLS.diet$residuals)
#Observed Phylogenetic Signal (K): 0.094


#1) use transformPhylo.ML to obtain lambda

summary(prcomp(GLS.diet$residuals))
# 10 PCs == 0.96
phy.sig <-transformPhylo.ML(prcomp(GLS.diet$residuals)$x[,1:10], phy.inf, lambdaEst=TRUE, model="lambda")
phy.sig$Lambda[1] #0.1976564


#2) rescale the branch lengths for that lambda value

phy.lambda <- rescale(phy.inf, "lambda", phy.sig$Lambda[1])


#3) repeat PGLS with the scaled phylogeny

gdf2=geomorph.data.frame(Diet=diet.phy.inf,Shape=mshapes.phy.inf,LCS.diet=data.phy$mean.Centroid.size,DietTree.scaled=phy.lambda)
Diet.PGLS.II_gdf2 <-procD.pgls(Shape~LCS.diet*Diet,DietTree.scaled,data=gdf2, SS.type="II")

physignal(phy= phy.lambda,GLS.diet$residuals)
#Observed Phylogenetic Signal (K): 0.9804


anova(Diet.PGLS.II_gdf2)



#PGLS FOR TREE AND GROUND SQUIRRELS SEPARATELY

#arboreal
gdf.arb=geomorph.data.frame(Diet= data.phy.inf.arb$Diet..this.study.,Shape= mshapes.phy.inf.arb,LCS.diet=mCsize.phy.inf.arb)

GLS.diet.arb=procD.lm(Shape~LCS.diet*Diet,data=gdf.arb)
GLS.diet.arb

physignal(phy=arb.phy,GLS.diet.arb$residuals)
#Observed Phylogenetic Signal (K): 0.0867


summary(prcomp(GLS.diet.arb$residuals))
# 10 PCs == 0.96
phy.sig_arb <-transformPhylo.ML(prcomp(GLS.diet.arb$residuals)$x[,1:10], arb.phy, lambdaEst=TRUE, model="lambda")
phy.sig_arb$Lambda[1] #0.08355203

phy.lambda_arb <- rescale(arb.phy, "lambda", phy.sig_arb$Lambda[1])


gdf2.arb=geomorph.data.frame(Diet=data.phy.inf.arb$Diet..this.study.,Shape=mshapes.phy.inf.arb,LCS.diet=mCsize.phy.inf.arb,DietTree.scaled=phy.lambda_arb)


Diet.PGLS.II_gdf2arb <-procD.pgls(Shape~LCS.diet*Diet,DietTree.scaled,data=gdf2.arb, SS.type="II")

physignal(phy= phy.lambda_arb,GLS.diet.arb$residuals)
#Observed Phylogenetic Signal (K): 1.0295

anova(Diet.PGLS.II_gdf2arb)





#terrestrial
gdf.terr=geomorph.data.frame(Diet= data.phy.inf.terr$Diet..this.study.,Shape=mshapes.phy.inf.terr,LCS.diet=mCsize.phy.inf.terr)

GLS.diet.terr=procD.lm(Shape~LCS.diet*Diet,data=gdf.terr)
GLS.diet.terr

physignal(phy=terr.phy,GLS.diet.terr$residuals)
#Observed Phylogenetic Signal (K): 0.1261


summary(prcomp(GLS.diet.terr$residuals))
# 8 PCs == 0.96
phy.sig_terr <-transformPhylo.ML(prcomp(GLS.diet.terr$residuals)$x[,1:8], terr.phy, lambdaEst=TRUE, model="lambda")
phy.sig_terr$Lambda[1] #1e-08
phy.lambda_terr <- rescale(terr.phy, "lambda",phy.sig_terr$Lambda[1])


gdf2.terr=geomorph.data.frame(Diet=data.phy.inf.terr$Diet..this.study.,Shape=mshapes.phy.inf.terr,LCS.diet=mCsize.phy.inf.terr,DietTree.scaled=phy.lambda_terr)

Diet.PGLS.II_gdf2terr <-procD.pgls(Shape~LCS.diet*Diet,DietTree.scaled,data=gdf2.terr, SS.type="II")

physignal(phy=phy.lambda_terr,GLS.diet.terr$residuals)
#Observed Phylogenetic Signal (K): 1

anova(Diet.PGLS.II_gdf2terr)




write.table(anova(Diet.PGLS.II_gdf2)$table, "PGLS.csv", dec=".", sep=";")
write.table(anova(Diet.PGLS.II_gdf2arb)$table, "PGLS.csv", dec=".", sep=";", append=T)
write.table(anova(Diet.PGLS.II_gdf2terr)$table, "PGLS.csv", dec=".", sep=";", append=T)



#pairwise comparisons

pairw <- pairwise(Diet.PGLS.II_gdf2, groups=diet.phy.inf)
pairw_summary <- summary(pairw, confidence = 0.95, test.type = "var", stat.table =TRUE)
pairw_summary$pairwise.tables$P
pairw_summary$summary.table$`Pr > d`
paiw_djust <- pairw_summary$pairwise.tables$D

paiw_djust[upper.tri(paiw_djust)] <- p.adjust(pairw_summary$pairwise.tables$P[upper.tri(pairw_summary$pairwise.tables$P)], method="holm")

write.table(paiw_djust, "pairwise_PGLS_ADJUSTED.csv", sep=";", dec = ".")





###### 3. Allometry plots ######


xy <- plot(Diet.PGLS.II_gdf2, type = "regression", predictor = log(data.phy$mean.Centroid.size), reg.type = "RegScore", col=col.diet.inf[diet.phy.inf], pch=16)
text(log(data.phy$mean.Centroid.size), xy$RegScore,labels=data.phy$Sp_red, cex= 0.5)

xy.arb <- plot(Diet.PGLS.II_gdf2arb, type = "regression", predictor = log(mCsize.phy.inf.arb), reg.type = "RegScore", col=col.diet.inf[-c(5,7)][data.phy.inf.arb$Diet..this.study.], pch=16)
text(log(mCsize.phy.inf.arb), xy.arb$RegScore,labels=data.phy.inf.arb$Sp_red, cex= 0.5)

xy.terr <- plot(Diet.PGLS.II_gdf2terr, type = "regression", predictor = log(mCsize.phy.inf.terr), reg.type = "RegScore", col=col.diet.inf[-c(3,4,6)][data.phy.inf.terr$Diet..this.study.], pch=16)
text(log(mCsize.phy.inf.terr), xy.terr$RegScore,labels=data.phy.inf.terr$Sp_red, cex= 0.5)








#### 4. Phylomorphospace ####
PCA.phy.inf <- gm.prcomp(A= arrayspecs(mshapes.phy.inf, 64,2), phy=phy.inf)


#PC1 and PC2
#X11()
pdf("phylomorphospace_PC1vsPC2.pdf", width = 11.69, height = 8.27)
plotGMPhyloMorphoSpace(phy.inf, A= mshapes.phy.inf, tip.labels = F, node.labels = F,
                       ancStates = F, xaxis = 1, yaxis = 2, zaxis = NULL,
                       plot.param=list(t.bg=col.diet.inf[diet.phy.inf], t.cex=data.phy$mean.Centroid.size*0.1, n.cex=0, lwd=0.5, txt.cex=0.6), shadow = FALSE)
legend("bottomleft", legend= levels(diet.phy.inf), fill =col.diet.inf, cex=0.6)
dev.off()

#PC1 and PC3
#X11()
pdf("phylomorphospace_PC1vsPC3.pdf", width = 11.69, height = 8.27)
plotGMPhyloMorphoSpace(phy.inf, A= mshapes.phy.inf, tip.labels = F, node.labels = F,
                       ancStates = F, xaxis = 1, yaxis = 3, zaxis = NULL,
                       plot.param=list(t.bg=col.diet.inf[diet.phy.inf], t.cex=data.phy$mean.Centroid.size*0.1, n.cex=0, lwd=0.5, txt.cex=0.6), shadow = FALSE)
legend("bottomleft", legend= levels(diet.phy.inf), fill =col.diet.inf, cex=0.6)
dev.off()




###with outlines: https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html


# Get generalized Procrustes coordinates
#gpa_array <- gpagen(lm_array)$coords

# Get generalized Procrustes coordinates
gpa_array <- arrayspecs(mshapes.phy.inf, 64,2)
dimnames(gpa_array)[[1]] <- c(1:64)
dimnames(gpa_array)[[2]] <- c("x", "y")



# Convert array to matrix for PCA
gpa_mat <- t(apply(gpa_array, 3, function(y) matrix(t(y), 1)))


# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores
scores <- gpa_mat %*% resEig$vectors
#scores <- scores*-1
# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

sum(per_var)


#define order of landmarks and semilandmarks to draw outlines
order <-c(1,5:19,2,20:34,3,35:49,4,50:64)


# Define function to draw shape
plot_squi_teeth <- function(xy, coor, size=1, col='black'){
  
  # If 3D, rotate points about x-axis using 3D rotation matrix
  if(ncol(coor) == 3){
    coor <- coor %*% matrix(c(1,0,0, 0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), nrow=3, ncol=3)
  }
  
  # Get just x,y coordinates (orthographic projection into xy-plane)
  coor <- coor[, 1:2]
  
  # Get plot aspect ratio
  w <- par('pin')[1]/diff(par('usr')[1:2])
  h <- par('pin')[2]/diff(par('usr')[3:4])
  asp <- w/h
  
  # Correct for plot aspect ratio not necessarily being 1:1
  coor[, 1] <- coor[, 1] * (1/asp)
  
  # Scale points and place back in position
  coor <- coor*size
  
  # Center about zero based on range of coordinates
  coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
                        nrow=nrow(coor), ncol=ncol(coor), byrow=TRUE)
  
  # Move shape to PC score
  coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow=TRUE)
  
  # Set order in which to draw points to create polygon
  polygon_order <- order
  
  # Create filled polygon
  polygon(coor[polygon_order, ], col=col, border=col)
}


# Set PCs to plot
pcs <- c(1,3)
pcs <- 1:2

# Open PDF graphics device
pdf('Backtransform PCA1_text.pdf', width=9, height=6.5)

# Create plot box with axes and axis labels
plot(scores[, pcs], type='n', main='Backtransform morphospace',
     xlab=paste0('PC', pcs[1], ' (', round(per_var[pcs[1]]), '%)'),
     ylab=paste0('PC', pcs[2], ' (', round(per_var[pcs[2]]), '%)'), asp = 1)
# Plot backtransform shapes
btShapes(scores=scores, vectors=resEig$vectors, fcn=plot_squi_teeth, 
         pcs=pcs, n=c(5,5), m=dim(lm_array)[2], row.names=dimnames(gpa_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.08, col=gray(0.7))
# Plot points for each species
points(scores[, pcs], bg=col.diet.inf[diet.phy.inf],col="transparent", cex=data.phy$mean.Centroid.size*0.1, lwd=0.5, pch=21)
# Add text labels
text(scores[, pcs], labels=substr(rownames(scores), 0, 3), cex=0.8, 
     pos=1, offset=0.3)

# Close the PDF graphics device
dev.off()





###### 5. Morphological disparity ######

##by diet
mdiet.inf <- data$Diet..this.study.
mdiet.infCVA <- data$Funtional_diet..this.study.
mshapes.inf <- as.matrix(data[,13:140])

# Morphological disparity of diets to the overall mean
morph.disp.mdiet <- morphol.disparity(mshapes.inf~1, groups=mdiet.inf, partial=FALSE)

morph.disp.mdiet.holm <- morph.disp.mdiet$PV.dist.Pval
morph.disp.mdiet.holm[1:8,1:8] <- p.adjust(morph.disp.mdiet$PV.dist.Pval, method= "holm")

write.table(morph.disp.mdiet.holm, "pairwise_holm_disparity.csv", sep=";", dec=".")


# Morphological disparity of species to the overall mean
morph.disp.mdiet.sp <- morphol.disparity(mshapes.inf~1, groups=rownames(mshapes.inf), partial=FALSE)



dfq <- as.data.frame(morph.disp.mdiet.sp$Procrustes.var)
colnames(dfq) <- "Dispar"
dfq$Diet <- mdiet.inf


q <- ggplot(dfq, aes(x=Diet, y=Dispar))




pdf("disparity_diet_geomorph.pdf", width = 11.69, height = 8.27)
q+
  geom_violin(trim = FALSE, color="transparent",fill="gray80", show.legend = FALSE) + 
  scale_color_manual(values = col.diet.inf)+
  geom_dotplot(
    aes(fill = Diet, color = Diet), trim = TRUE,
    binaxis='y', stackdir='center', dotsize = 0.4,
    position = position_dodge(0.8), show.legend = FALSE, method="dotdensity")+
  scale_fill_manual(values = col.diet.inf)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
               geom = "pointrange", color = "black")+
  theme_classic()
dev.off()




###### +size by dietary categories v ######


summary(aov(mCsize.inf~ mdiet.inf))
bonf <- pairwise.wilcox.test(mCsize.inf, mdiet.inf, p.adjust.method = "holm")
symnum(bonf$p.value, cutpoints = c(0, 0.05, 1), symbols=c("**", " "))


write.table(bonf$p.value, "disparity/pairwise_holm_size.csv", sep=";", dec=".")



e <- ggplot(data.frame(V1=data.phy$mean.Centroid.size, diet=diet.phy.inf), aes(x=diet, y=V1))


pdf("disparity/size_diet1.pdf", width = 11.69, height = 8.27)
e+
  geom_violin(trim = FALSE, fill="gray80", color="transparent", show.legend = FALSE) + 
  scale_color_manual(values = col.diet.inf)+
  geom_dotplot(
    aes(fill = diet, color = diet), trim = FALSE,
    binaxis='y', stackdir='center', dotsize = 0.4,
    position = position_dodge(0.8), show.legend = FALSE)+
  scale_fill_manual(values = col.diet.inf)+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
               geom = "pointrange", color = "black")+
  geom_hline(yintercept=mean(data.phy$mean.Centroid.size))+
  #geom_text(label=rownames(data.phy$mean.Centroid.size))+
  theme_classic()
dev.off()




###### 5. MvMorph ######


  #simpler models
fitBM1 <- mvBM(phy.inf,  PCA.phy.inf$x[,1:6], model = "BM1", param=list(decomp="diagonal"))
saveRDS(fitBM1, file = "p4fitBM1.RData")
fitOU1 <- mvOU(phy.inf, PCA.phy.inf$x[,1:6], model = "OU1", param=list(decomp="diagonal"))
saveRDS(fitOU1, file = "p4fitOU1.RData")
fitEB <- mvEB(phy.inf, PCA.phy.inf$x[,1:6], param=list(decomp="diagonal"))
saveRDS(fitEB, file = "p4fitEB.RData")


  #multi-peak models

#simmap trees reconstructions (reconstruction of diet)


#Diet 8

#geiger fitDiscrete
Diet_ER=fitDiscrete(phy.inf, diet.phy.inf,model="ER")
Diet_SYM=fitDiscrete(phy.inf, diet.phy.inf,model="SYM")
Diet_ARD=fitDiscrete(phy.inf, diet.phy.inf,model="ARD")

#The smaller the better
Diet_ER$opt$aicc
Diet_SYM$opt$aicc
Diet_ARD$opt$aicc

#rec
mtree <- make.simmap(tree=ladderize(phy.inf, right=F), x=diet.phy.inf, model="ER", nsim=1000)

pd <- describe.simmap(mtree)
col.pd <- col.diet.inf
names(col.pd) <- levels(data.phy$Diet..this.study.)


#Diet CVA

col.pdCVA <- col.dietCVA.inf
names(col.pdCVA) <- levels(data.phy$Funtional_diet..this.study.)


#geiger fitDiscrete
Diet_ER_CVA=fitDiscrete(phy.inf, fundiet.phy.inf,model="ER")
Diet_SYM_CVA=fitDiscrete(phy.inf, fundiet.phy.inf,model="SYM")
Diet_ARD_CVA=fitDiscrete(phy.inf, fundiet.phy.inf,model="ARD")

#The smaller the better
Diet_ER_CVA$opt$aicc
Diet_SYM_CVA$opt$aicc
Diet_ARD_CVA$opt$aicc

mtreeCVA <- phytools:::make.simmap(tree=ladderize(phy.inf, right=F), x=fundiet.phy.inf, model="SYM", nsim=1000)
pdCVA <- describe.simmap(mtreeCVA)



###

#Function to transform ancestral states reconstruction to a SIMMAP like tree (to be added to mvMORPH. Please cite this package accordingly)
source("paintAllTree.R") #see https://github.com/EllenJCoombs/Asymmetry-evolution-cetaceans
mv_pd <- paintAllTree(ladderize(phy.inf, right=F), pd, as.character(as.vector(diet.phy.inf)))
mv_pdCVA <- paintAllTree(ladderize(phy.inf, right=F), pdCVA, as.character(as.vector(fundiet.phy.inf)))



#diet
fitBMM.diet <- mvBM(mv_pd, PCA.phy.inf$x[,1:6], model = "BMM", param=list(decomp="diagonal"))
saveRDS(fitBMM.diet, file = "fitBMM.diet.RData")
fitOUM.diet <- mvOU(mv_pd, PCA.phy.inf$x[,1:6], model = "OUM", param=list(decomp="diagonal"))
saveRDS(fitOUM.diet, file = "fitOUM.diet.RData")


#for functional diet

fitBM.cva <- mvBM(mv_pdCVA, PCA.phy.inf$x[,1:6], model = "BMM", param=list(decomp="diagonal"))
saveRDS(fitBM.cva, file = "cva_p4fitBM.RData")
fitOUM.cva <- mvOU(mv_pdCVA, PCA.phy.inf$x[,1:6], model = "OUM", param=list(decomp="diagonal"))
saveRDS(fitOUM.cva, file = "cva_p4fitOUM.RData")


#paint tree with PhyloEM peaks

OU_EM_tree <- ladderize(phy.inf, right = F)
OU_EM_tree <- paintSubTree(OU_EM_tree ,OU_EM_tree$edge[1], 0)
OU_EM_tree <- paintSubTree(OU_EM_tree , getMRCA(OU_EM_tree, c("Tam_dor", "Sci_dav")), 1)
OU_EM_tree <- paintSubTree(OU_EM_tree ,getMRCA(OU_EM_tree, c("Spe_dau", "Xer_spi")), 2)
OU_EM_tree <- paintSubTree(OU_EM_tree , getMRCA(OU_EM_tree, c("Cyn_leu", "Cyn_gun", "Cyn_par", "Cyn_lud") ),3)

plot(OU_EM_tree)
OU_EM_tree$mapped.edge
OU_EM_tree$tip.label
states <- data.frame(tip=OU_EM_tree$tip.label, stringsAsFactors = FALSE)
rownames(states) <- OU_EM_tree$tip.label
state1 <- extract.clade(OU_EM_tree, getMRCA(OU_EM_tree, c("Tam_dor", "Sci_dav")))$tip.label
state2 <- extract.clade(OU_EM_tree, getMRCA(OU_EM_tree, c("Spe_dau", "Xer_spi")))$tip.label
state3 <- extract.clade(OU_EM_tree, getMRCA(OU_EM_tree, c("Cyn_leu", "Cyn_gun", "Cyn_par", "Cyn_lud")))$tip.label
states$tip <- 0
states$tip[rownames(states) %in% state1] <- 1
states$tip[rownames(states) %in% state2] <- 2
states$tip[rownames(states) %in% state3] <- 3




fitOU_EM4 <- mvOU(OU_EM_tree, PCA.phy.inf$x[,1:6], model = "OUM", param=list(decomp="diagonal"))
saveRDS(fitOU_EM4, file = "EM_p4fitOU.RData")
fitBM_EM4 <- mvBM(OU_EM_tree, PCA.phy.inf$x[,1:6], model = "BMM", param=list(decomp="diagonal"))
saveRDS(fitBM_EM4, file = "EM_p4fitBM.RData")

mvMORPH_squirrels.EM <- data.frame(AICc = fitOU_EM4$AICc)
mvMORPH_squirrels.EM [2,] <- fitBM_EM4$AICc
row.names(mvMORPH_squirrels.EM ) <- c("OUM", "BM4")
write.table(mvMORPH_squirrels.EM , file = "EM_mvMORPH_squirrels.csv", append = FALSE, quote = F, sep=";", dec=".", col.names = T, row.names = T)



###RESULTS###
results <- list(fitBM1, fitOU1, fitEB, fitBMM.diet,fitOUM.diet, fitBM.cva, fitOUM.cva, fitBM_EM4, fitOU_EM4)
aicw(results, aicc=TRUE)

results2 <- list(fitBM1, fitOU1, fitEB, fitBMM.diet,fitOUM.diet, fitBM.cva, 
                fitOUM.cva, fitBM_EM4, fitOU_EM4, 
                fitBM.zel, fitOUM.zel,
                fitBM.mclean, fitOUM.mclean)
aicw(results2, aicc=TRUE)


##### +++simulations of optima ####


##THIS SIMULATIONS TAKE LONG TIME###

set.seed(15)
nsim=50

# Alternatively, we can fit using mvSIM rather than a wrapper:
bootSampleOU1 <- mvSIM(param=fitOU1, tree=OU_EM_tree, nsim=nsim)
bootSampleOU4 <- mvSIM(param=fitOU_EM4, tree=OU_EM_tree, nsim=nsim)
bootSampleOU8 <- mvSIM(param=fitOUM.diet, tree=mv_pd, nsim=nsim)


#fitou8_sims <- list()
for (i in 1:50) {
  fitou8_sims[[i]] <- mvOU(mv_pd, bootSampleOU8[[i]], 
                           model="OUM", method="sparse", 
                           echo=FALSE, diagnostic=FALSE, param=list(decomp="diagonal"))
  print(i)
}

saveRDS(fitou8_sims, "fitou8_sims_simulations.RData")

#fitou8_sims <- readRDS("fitou8_sims_simulations.RData")



 #fit_phyOU4_sims <- list()
for (i in 1:50) {
  fit_phyOU4_sims[[i]] <- mvOU(OU_EM_tree, bootSampleOU4[[i]], 
                           model="OUM", method="sparse", 
                           echo=FALSE, diagnostic=FALSE, param=list(decomp="diagonal"))
  print(i)
}

saveRDS(fit_phyOU4_sims, "fit_phyOU4_sims_simulations.RData")




#fitou1_sims <- list()
for (i in 21:50) {
  fitou1_sims[[i]] <- mvOU(OU_EM_tree, bootSampleOU1[[i]], 
                               model="OUM", method="sparse", 
                               echo=FALSE, diagnostic=FALSE, param=list(decomp="diagonal"))
  print(i)
}

saveRDS(fitou1_sims, "fitou1_sims_simulations.RData")






optima_OUM <- data.frame(states= NA, PC1=NA, PC2=NA, PC3=NA)
a <- 1
for (i in 1:10) {
  optima_OUM[a:(a+3),2:4] <- fitoum_sims[[i]]$theta[,1:3,drop=TRUE]
  optima_OUM[a:(a+3),1] <- rownames(fitoum_sims[[i]]$theta[,1:3,drop=TRUE])
  a <- a+4
}

optima_OU1 <- data.frame(states= 1, PC1=NA, PC2=NA, PC3=NA)
a <- 1
for (i in 1:10) {
  optima_OU1[a:(a+1),2:4] <- fitou1_sims[[i]]$theta[,1:3,drop=TRUE]
  optima_OU1[a:(a+1),1] <- 1
  a <- a+2
}

optima_OU8 <- data.frame(states= 1, PC1=NA, PC2=NA, PC3=NA)
a <- 1
for (i in 1:10) {
  optima_OU8[a:(a+7),2:4] <- fitou8_sims[[i]]$theta[,1:3,drop=TRUE]
  optima_OU8[a:(a+7),1] <- rownames(fitou8_sims[[i]]$theta[,1:3,drop=TRUE])
  a <- a+8
}








plot_density1 <- ggplot(data=optimos_OUM, aes(x=PC1,y=PC2))+
  stat_density2d(aes(alpha = ..level.., fill = states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_manual(values=c("black","#800000", "#008000","#000080"))+
  coord_fixed()+
  xlim(-0.11,0.25)+
  ylim(-0.06,0.05)+
  ggtitle("Phy OU4 optima")+
  theme_classic()
plot_density2 <- ggplot(data=optimos_OUM, aes(x=PC1,y=PC3))+
  stat_density2d(aes(alpha = ..level.., fill = states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_manual(values=c("black","#800000", "#008000","#000080"))+
  coord_fixed()+
  xlim(-0.11,0.25)+
  ylim(-0.06,0.05)+
  ggtitle("Phy OU4 optima")+
  theme_classic()


plot_density4 <- ggplot(data=optimos_OU1, aes(x=PC1,y=PC2))+
  stat_density2d(aes(alpha = ..level.., fill=states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_continuous(low = "grey", high = "gray")+
  coord_fixed()+
  xlim(-0.15,0.25)+
  ylim(-0.1,0.08)+
  ggtitle("OU1 optimum")+
  theme_classic()

plot_density5 <- ggplot(data=optimos_OU1, aes(x=PC1,y=PC3))+
  stat_density2d(aes(alpha = ..level.., fill=states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_continuous(low = "grey", high = "gray")+
  coord_fixed()+
  xlim(-0.15,0.25)+
  ylim(-0.1,0.08)+
  ggtitle("OU1 optimum")+
  theme_classic()

plot_density6 <- ggplot(data=optimos_OU1, aes(x=PC2,y=PC3))+
  stat_density2d(aes(alpha = ..level.., fill=states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_continuous(low = "grey", high = "gray")+
  coord_fixed()+
  xlim(-0.15,0.25)+
  ylim(-0.1,0.08)+
  ggtitle("OU1 optimum")+
  theme_classic()

plot_density7 <- ggplot(data=optima_OU8, aes(x=PC1,y=PC2))+
  stat_density2d(aes(alpha = ..level.., fill = states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_manual(values=col.diet.inf)+
  coord_fixed()+
  xlim(-0.11,0.25)+
  ylim(-0.06,0.05)+
  ggtitle("OU8 optima")+
  theme_classic()

plot_density8 <- ggplot(data=optima_OU8, aes(x=PC1,y=PC3))+
  stat_density2d(aes(alpha = ..level.., fill = states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_manual(values=col.diet.inf)+
  coord_fixed()+
  xlim(-0.11,0.25)+
  ylim(-0.06,0.05)+
  ggtitle("OU8 optima")+
  theme_classic()

plot_density9 <- ggplot(data=optima_OU8, aes(x=PC2,y=PC3))+
  stat_density2d(aes(alpha = ..level.., fill = states),geom= "polygon",alpha=0.3, bins=5, show.legend = FALSE)+
  scale_fill_manual(values=col.diet.inf)+
  coord_fixed()+
  xlim(-0.11,0.25)+
  ylim(-0.06,0.05)+
  ggtitle("OU8 optima")+
  theme_classic()












