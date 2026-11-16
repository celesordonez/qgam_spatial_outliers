"
This program uses Quantile GAM (Generalized Additive Models) to detect and correct positive spatial outliers
The lithology is considered as a categorical covariate
Autho: Celestino Ordóñez
"


# install.packages('Rcpp', dependencies = TRUE)
# install.packages('ggplot2', dependencies = TRUE)
# install.packages('data.table', dependencies = TRUE)
# install.packages('plotly', dependencies = TRUE)
# install.packages("scatterplot3d", dependencies = TRUE )
#install.packages('rgl', dependencies = TRUE)
#install.packages("spdep", dependencies = TRUE )

#install.packages("smoothr", dependencies = TRUE )
#install.packages("energy", dependencies = TRUE ) # distance correlation
#install.packages("plot3D", dependencies = TRUE )
#install.packages("visreg", dependencies = TRUE )

library(readxl)
library(qgam)
library(MASS)
library(plotly)
library(scatterplot3d)
library(mgcViz)
library(rgl)
library(utils)
library(spdep)
library(GeoModels)
library(dbscan) # for k nearest neighbors
library(energy)
library(ggplot2)
library(plot3D)
library(visreg)

# Set the working directory
setwd('D:/ARTICULOSrecientes/carlos_sierra')
# Read the csv file containing the data
data_lito <- read.csv('element_lito_id.csv')
which(is.na(data_lito), arr.ind=TRUE)
# Define ID_UNIDAC attribute as factor
data_lito$ID_UNIDADC <- as.factor(data_lito$ID_UNIDADC)
# is.factor(data_lito$ID_UNIDADC)

# Descriptive statistics
summary(data_lito)
# Box plots
boxplot(data_lito[,3:8], outline=FALSE)

# Distance correlation betwween lito and different metals
#dcor(data_lito$ID_UNIDADC, data_lito$Cd_ppm)
#dcor(data_lito$SIMB_LITO, data_lito$Cd_ppm)
#cor(data_lito$Al_percent, data_lito$Cd_ppm)
#cor(data_lito$Al_percent, data_lito$Cd_ppm)
#plot(data_lito$Al_percent, data_lito$Cd_ppm)

# Test if average concentration of a metal is related to lithology
kruskal.test(data_lito$Al_percent~data_lito$ID_UNIDADC) # la hipótesis nula sería que la concentración media de Al no depende de la litología
# kruskal.test(data_lito$Al_percent~data_lito$SIMB_LITO) # la hipótesis nula sería que la concentración media de Al no depende de la litología
x = data_lito$ID_UNIDADC  # independent variable
y = data_lito$Al_percent  # response variable

#p = ggplot(data_lito) + geom_boxplot(aes(x, y), outlier.shape=NA)
# Plot the distribution of some of the pollutants for the different lithologies
p = ggplot(data_lito) + geom_boxplot(aes(x, y))
p1 = p + coord_cartesian(ylim = c(0,15)) + xlab('Lithology') + ylab('Al%')
print(p1)
# Box plot
boxplot(data_lito[,3:8], outline=FALSE)

hist(data_lito$Cu_ppm)  # Cu distribution
# Additve model with the lithology as covariate
model_fac <-gam(Cu_ppm ~ ID_UNIDADC  + s(Al_percent, by=ID_UNIDADC), bs="fs", family=Gamma(link=log), data = data_lito, method = "REML", control = gam.control(keepData = TRUE))
# Additive model without using the lithology a covariate
model<-gam(Cu_ppm ~  s(Al_percent), bs="fs", family=Gamma(link=log), data = data_lito, method = "REML", control = gam.control(keepData = TRUE))

# Compare both models
summary(model_fac)
summary(model)
gam.check(model_fac)
plot(model_fac, ylim=c(-5,5))
anova(model, model_fac)

# Visualize the correlation between an element and Al depending on the lithology
visreg(model_fac, xvar = "Al_percent",
       by = "ID_UNIDADC", data =  data_lito, alpha=0.5, ylim=c(-10,10), breaks=c(73,74,75,76), ylab="Cd")

# Visualize the correlation of an element with Al for a specific lithology
visreg(model_fac, xvar = "Al_percent",
       data =  data_lito, cond=list(ID_UNIDADC=70), alpha=0.5, ylim=c(-1,8), xlab="Al|Litho=70", ylab="Cu")

# Visualize the effect of Al
visreg(model_fac, 'Al_percent', type='conditional')
# Visualize the effects of Al depending on the lithology
visreg(model_fac, 'Al_percent', by='ID_UNIDADC', ylim=c(-1,5))
# 2D representation of the effect
visreg2d(model_fac, 'Al_percent', 'ID_UNIDADC', plot.type="image")

# Represent the effect for a each lithology
b<-getViz(model_fac)
o <- plot( sm(b, 6) )
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

#  All together
listLayers(o)
print(plot(b, allTerms = T), pages = 1)


##################################################################
# Here starts the algorithm to detect outliers from the residuals #
################################################################
res <- residuals(model_fac)
k=6  # number of neihgbors
alpha=0.05  # significance level
coords=data_lito[,1:2]  # coordinates of the observations
outlier=GeoOutlier(data=res, coordx = coords,neighb=k, alpha=alpha)  # detect the outliers
outlier # outliers
quilt.plot(coords, data_lito$Pb_ppm) # represent the data
# Extract the indices of the rows of where there are outliers
ids<-as.integer(row.names(outlier))
# Plot the outliers on top of the observations
for (i in 1:nrow(outlier))  plotrix::draw.circle(outlier[i,1], outlier[i,2],radius=2.5, lwd=6)
# Replace the outliers values for the element under study with the median of its nearest neighbors
nn <- kNN(x=coords, k = k, query=coords[ids,])  # nearest neighbors
#nn <- kNN(x=coords, k = 5, query=coords[c(4,373,377,378,381,382),])
plot(nn,coords)  # plot the k nearest neighbors around the outliers 
nn$id
nn$dist
outliers <- data_lito[ids,]
weights <- 1/nn$dist[,2:k]/sum(1/nn$dist[,2:k])  # weights to correct the outliers values
sum(weights)  # must be 1
# Create a new dataframe to obtain the corrected values
data_lito_corr <- data_lito
# Compare the value of the outlier with median its neighbor. It could be a positive or negative value
conci<-vector()
medt<-vector()
outlier_pos<-vector()
n = NROW(ids)
for (i in 1:n){
  for (j in 1:k-1){
    conci[j]<-data_lito$Zn_ppm[nn$id[i,j+1]]
  }
  medt[i]<-mean(conci)
  if(data_lito$Zn_ppm[ids[i]]>medt[i]){outlier_pos[i]<-ids[i]}
}

outlier_pos<-as.numeric(na.omit(outlier_pos))
outlier_pos  # positive outliers

### Weighted mean of the positive outliers ###
nnp <- kNN(x=coords, k = k, query=coords[outlier_pos,])  # nearest neighbors
m = NROW(outlier_pos)
weightsp<-matrix(data=NA, nrow=m, ncol=k)
weightsp[,1]<-0  # the own observation has weight 0
mat1<-1/nnp$dist[,2:k]
mat2<-as.matrix(rowSums(mat1))
mat3<-rep(mat2, times=k-1)  # transform mat3 to a matrix (mxk-1)
mat5<-mat1/mat4
weightsp[,2:k] <- mat5  # matrix of weigths

# Correct the values of the outliers
for (i in 1:m){
    data_lito_corr$Ni_ppm[outlier_pos[i]]= sum(data_lito$Ni_ppm[nnp$id[i,]]*weightsp[i,])
  }

######################## Quantile GAM #######################################################
# Instead of estimating the mean value of y we could estimate a different quantile using qgam #
#############################################################################################
# Fit the quantile model using only the coordinates as explanatory variabes
fitxy <- qgam(Zn_ppm ~s(XUTM, YUTM), data = data_lito_corr, qu=0.5)
summary(fitxy)
# pred <- predict(fitxy, newdata = data, se=FALSE)
# xSeq <- data.frame(cbind("Zn_ppm" = rep(0, 1e3), "XUTM" = seq(365500, 486350, length.out = 1e3), "YUTM" = seq(4735600, 4817400, length.out = 1e3)))
# df1 <- data.frame(pred, xSeq$XUTM, xSeq$YUTM)
# plot_ly(z = ~xtabs(pred ~ data.XUTM + data.YUTM, data = df1)) %>% add_surface()

# Export model predictions
x <- seq(349500, 487966, length.out = 1e2)
y <- seq(4734618, 4818240, length.out = 1e2)
coords <- expand.grid(x=x, y=y)
#z <- rep(0, 1e4)
#df_grid <- data.frame(XUTM=x,YUTM=y,As_ppm=z)
df_grid <- data.frame(XUTM=coords$x,YUTM=coords$y)
#xSeq <- data.frame(cbind("As_ppm" = rep(0, 1e3), "XUTM" = seq(365500, 486350, length.out = 1e3), "YUTM" = seq(4735600, 4817400, length.out = 1e3)))
predxy <- predict(fitxy, newdata = df_grid, k=20, se=TRUE)
#zf<-matrix(predxy$fit,nrow=100, ncol=100)
# Export to .csv format the prediction and the standar error
el_gbg<-data.frame(XUTM=coords$x,YUTM=coords$y,Zn_ppm=predxy)
write.csv(el_gbg,'Cu_ppm.csv')


# Calculate confidence intervals for the prediction and represent them in 3D perspective
zf <- matrix(predxy$fit, nrow=100, byrow=TRUE)
zf_sup <- matrix(predxy$se.fit, nrow=100, byrow=TRUE)
zf_inf <- - matrix(predxy$se.fit, nrow=100, byrow=TRUE)
# open3d()
# surface3d(x,y,zf, alpha = 0.4, front = "lines", back = "lines")
# fig <-plot_ly(type = 'surface', x=as.matrix(xSeq$XUTM), y=as.matrix(xSeq$YUTM), z=as.matrix(predxy))%>%add_trace()
# fig # No 
# 3D representation of the prediction
persp(x=x,y=y,z=zf*1000,theta=45, phi=20, axes=TRUE, shade = 0.7, expand = 0.9, d = 1, scale=FALSE, box=TRUE, nticks=4, ticktype="detailed", xlab="X", ylab="Y", zlab="Z")

col.pal<-colorRampPalette(c("blue" ,"green", "yellow", "red"))
colors<-col.pal(300)
# Height of facets
z.facet.center <- (zf[-1, -1] + zf[-1, -ncol(zf)] + zf[-nrow(zf), -1] + zf_sup[-nrow(zf), -ncol(zf)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 300)

# Try another vertical scale
persp(x=x,y=y,z=zf*200,theta=45, phi=20, axes=TRUE, shade = 0.7, expand = 0.9, d = 1, scale=FALSE, box=TRUE, nticks=4, col=colors[z.facet.range], ticktype="detailed", xlab="X", ylab="Y", zlab="Z")
image.plot(legend.only=T, zlim=range(z.facet.center), col=colors)

# Represent the upper limit of the confidence interval
col.pal<-colorRampPalette(c("yellow", "red"))
colors<-col.pal(300)
# height of facets
z.facet.center <- (zf_sup[-1, -1] + zf_sup[-1, -ncol(zf_sup)] + zf_sup[-nrow(zf_sup), -1] + zf_sup[-nrow(zf_sup), -ncol(zf_sup)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 300)
persp(x=x,y=y,z=zf*200,theta=45, phi=20, axes=FALSE, shade = 0.7, expand = 0.9, d = 1, scale=FALSE, box=TRUE, nticks=4, col=colors[z.facet.range], ticktype="detailed", xlab="X", ylab="Y", zlab="Z")
image.plot(legend.only=T, zlim=range(z.facet.center), col=colors)

#########################################################################################
### In this part we don´t take into account the lithology
########################################################################################
data <-read_xlsx('Datos_cant.xlsx', 2)  # data without lithology

# scatter3D plot
data2 <- data[, c(1,2,4,3)]
scatterplot3d(data2[,1:3])

# Histogram
hist(data$As_ppm)  
hist(log(data$As_ppm))  

# Regresión simple
ols <- lm(As_ppm ~Al_percent, data=data)
plot(data$As_ppm, rstandard(ols))
abline(h=0)
summary(ols)
summary(ols)$sigma

# Robust linear regression
rols <- rlm(As_ppm ~Al_percent, data=data)
summary(rols)
summary(rols)$sigma


# Effect of Al in As
fit <- qgam(As_ppm ~s(Al_percent), data = data, qu=0.5)
summary(fit)
plot(fit,pages=1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2)
res <- residuals(fit)
pred <- predict(fit, newdata = data, se=FALSE)
res1 <- data$As_ppm - pred

plot(data$As_ppm, res)
cor(data$As_ppm, res)  # There is a strong correlation between the data and the residuals
# Calculate outliers on the residuals
k=5
alpha=0.002
coords=data[,1:2]
outlier=GeoOutlier(data=data$As_ppm, coordx = coords,neighb=k, alpha=alpha)
outlier
quilt.plot(coords, data$As_ppm)
for (i in 1:nrow(outlier))  plotrix::draw.circle(outlier[i,1], outlier[i,2],radius=2.5, lwd=6)
# replace the values of As_ppm on the outliers for the median of their neighbors
nn <- kNN(coords, coords[c(101,141,214,238),], k = 5)
nn$id
nn$dist
outliers <- data[c(101,141,214,238),]
weights <- 1/nn$dist[,2:5]/sum(1/nn$dist[,2:5])
sum(weights)
### Weighted mean ###
id <- c(101, 141, 214, 238)
for (i in 1:4){
  for (j in 2:5){
  data$As_ppm[id[i]]= sum(data$As_ppm[nn$id[i,j]]*weights[i,j])
  }
}

# Fit As using the spatial coordinates as covariates
fit <- qgam(As_ppm ~s(XUTM,YUTM), data = data, qu=0.5)
summary(fit)
plot(fit,pages=1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2)
vis.gam(fit)
xSeq <- data.frame(cbind("As_ppm" = rep(0, 1e3), "XUTM" = seq(365500, 486350, length.out = 1e3), "YUTM" = seq(4735600, 4817400, length.out = 1e3)))
predxy <- predict(fitxy, newdata = xSeq, se=FALSE)

# Add Al as covariate
fit <-qgam(As_ppm ~s(XUTM,YUTM)+ s(Al_percent), data = data, qu=0.5)
summary(fit)

pred <- predict(fit, newdata = data, se=TRUE)

res <- data$As_ppm - pred$fit

# https://www.rdocumentation.org/packages/mgcv/versions/1.9-0/topics/plot.gam
# Plot the effects of both spatial coordinates Al_percent
plot(fit,pages=1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2)

# 2D representation of the effect
fviz = getViz(fit)
# o <- plot( sm(fviz, 1) )
# o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=x), alpha = 0.8) +
#   l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
#   l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

# https://cran.r-project.org/web/packages/mgcViz/vignettes/mgcviz.html
plot(sm(fviz, 1)) + l_fitRaster() + l_fitContour() + l_points()

# Different graphics that help to understand the results
op <- par(mfrow=c(1,2))
plot(fit,all.terms=TRUE)
par(op) 
op <- par(mfrow=c(1,2))
plot(fit,all.terms=TRUE,scheme=1)
par(op)
op <- par(mfrow=c(1,2))
plot(fit,all.terms=TRUE,scheme=c(2,1))
par(op)

# Another way                                                                                                                      
plot(fit,pages=1) 

# 3D representation
vis.gam(fit, theta = 120, n.grid = 50, lwd = .4)

xSeq <- data.frame(cbind("As_ppm" = rep(0, 384), "XUTM" = seq(365500, 486350, length.out = 384), "YUTM" = seq(4735600, 4817400, length.out = 384), 'Al_percent' = seq(0.39, 18, length.out = 384)))
plot(data$Al_percent, data$As_ppm, xlab  = "XUTM", ylab = 'As', ylim = c(3.5, 132))
#pred <- predict(fit, newdata = xSeq, se=TRUE)


# Divide between train y test
# library(caTools)
# 
# ind = sample.split(Y=data$As_ppm, SplitRatio = 0.7)
# 
# train = data[ind,]
# dim(train)
# test = data[!ind,]
# dim(test)
# 
# pred <- predict(fit, newdata=test, se=TRUE)
# plot(test$As_ppm, pred$fit)
# 
# xSeq <- data.frame(cbind("As_ppm" = rep(0, 384), "XUTM" = seq(365500, 486350, length.out = 384), "YUTM" = seq(4735600, 4817400, length.out = 384)))
# pred <- predict(fit, newdata = xSeq, se=TRUE)
# 
# plot(xSeq$XUTM, data$As_ppm, xlab = "XUTM", ylab = "As_ppm", ylim = c(2, 132))
# lines(xSeq$XUTM, pred$fit, lwd = 1)
# lines(xSeq$XUTM, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
# lines(xSeq$XUTM, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)
# 
# 
# plot(test$As_ppm, pred)
# 
# xSeq <- data.frame(cbind("As_ppm" = rep(0, 1e3), "XUTM" = seq(365500, 486350, length.out = 1e3), "YUTM" = seq(4735600, 4817400, length.out = 1e3)))
# pred <- predict(fit, newdata = xSeq, se=TRUE)
# 
# plot(data$XUTM, data$As_ppm, xlab = "YUTM", ylab = "As_ppm", ylim = c(2, 132))
# lines(xSeq$XUTM, pred$fit, lwd = 1)
# lines(xSeq$XUTM, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
# lines(xSeq$XUTM, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)

# fig <-plot_ly(x=xSeq$XUTM, y=xSeq$XUTM, z=pred$fit)
# fig

#### Try a GLM approach #####

ml<-glm(Cu_ppm ~ Al_percent:ID_UNIDADC, data = data_lito, family=gaussian(link='log'))
summary(ml)
plot(ml)
visreg(ml, xvar = "Al_percent",
       by = "ID_UNIDADC", data =  data_lito, ylim=c(-10,10))