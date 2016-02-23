#' Notes:
#' . The main advantage of copulas is that univariate
#' marginal distributions can be defined independently
#' from the joint behaviour of the variables involved. Hence, a
#' copula allows for modelling the dependence structure of random
#' variables regardless the family that the marginal distributions
#' belong to. Besides, joint return periods can be easily
#' estimated from copulas, which represents an additional benefit
#' as the study of joint return periods is essential to flood
#' frequency analysis.
#' The theory of copulas is based on the Sklar’s theorem
#' (Sklar, 1959), which in the case of a bivariate case can be
#' written in the form:
#'  H (x, y) = C{F (x), G(y)}, x, y ∈ R, (1)
#' where H (x, y) is the joint cumulative distribution function
#' of the random variables X and Y , F (x) and G(y) are the
#' marginal distribution functions of X and Y , respectively,
#' and the mapping function C : [0, 1]2 →[0, 1] is the copula
#' function.
#' Further details about copulas can be found in Joe (1997),
#' Nelsen (1999) and Salvadori et al. (2007).
#' Although copula models have been extensively applied in
#' other fields such as finance, they have been only recently applied
#' to model hydrological events such as floods, storms
#' and droughts. Overall, the Archimedean and extreme value
#' copula families are the most used in modelling flood variables.
#' The Archimedean copulas can be constructed easily
#' and, as a great deal of copulas belongs to this family, a broad
#' kind of dependence can be considered. Some authors used
#' Archimedean copulas such as the Frank copula (Favre et al.,
#' 2004) or the Clayton copula (Shiau et al., 2006) to characterise
#' the dependence structure between peak and volume variables. 
#' Meanwhile, extreme value copulas have the advantage
#' that they are able to connect the extreme values of the
#' studied variables, which is very important in flood frequency
#' analysis. A lot of authors considered the Gumbel copula as
#' the copula that best represents the relation between peak and
#' volume (Zhang and Singh, 2006, among others).

#' Flood peak = Q, Volume = V, Discharge = R
library(copula)
library(spcopula)
library(VineCopula)
library(evd)
library(MASS)
library(hexbin)
library(rgl)

load("~/R/ASTDProject/simTriHydro.RData")
head(all_QVD)
str(all_QVD)
simTriHydro <- all_QVD[,-1]
head(simTriHydro)
colnames(simTriHydro) <- c("Peak","Volume", "Duration")
head(simTriHydro)
plot(simTriHydro$Peak,simTriHydro$Volume)
abline(lm(simTriHydro$Volume~simTriHydro$Peak), col = "blue", lwd = 3)
plot3d(simTriHydro)

# Function to show marginal distribution with histograms along with scatter plot

scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

scatterhist(simTriHydro$Peak,simTriHydro$Volume,xlab = "Peak", ylab = "Volume")
scatterhist(simTriHydro$Volume,simTriHydro$Duration,xlab = "Volume", ylab = "Duration")
scatterhist(simTriHydro$Peak,simTriHydro$Duration,xlab = "Peak", ylab = "Duration")


# Compute the pseudo-observations for simTriHydro and check scatter plot to find correlations
hydro <- pobs(simTriHydro[])
plot(cbind(hydro[,1], hydro[,2]), cex=0.5, asp=1, xlab = "Peak", ylab = "Volume")

# Fit different multivariate copulas to find their correlation and dependencies
gaussFit <- fitCopula(normalCopula(dim=3), hydro, method = "ml")
gaussFit

u1 = rCopula(200,normalCopula(coef(gaussFit),dim=3))
points(u1[,1],u1[,2],col = "grey")

gumbFit <- fitCopula(gumbelCopula(dim=3), hydro)
gumbFit
plot(cbind(hydro[,1], hydro[,2]), cex=0.5, asp=1, xlab = "Peak", ylab = "Volume")
u1 = rCopula(500,gumbelCopula(coef(gumbFit),dim=3))
points(u1[,1],u1[,2],col = "grey")


clayFit <- fitCopula(claytonCopula(dim=3), hydro)
clayFit
plot(cbind(hydro[,1], hydro[,2]), cex=0.5, asp=1, xlab = "Peak", ylab = "Volume")
u1 = rCopula(500,claytonCopula(coef(clayFit),dim=3))
points(u1[,1],u1[,2],col = "grey")

franFit <- fitCopula(frankCopula(dim=3), hydro)
franFit
plot(cbind(hydro[,1], hydro[,2]), cex=0.5, asp=1, xlab = "Peak", ylab = "Volume")
u1 = rCopula(500,frankCopula(coef(franFit),dim=3))
points(u1[,1],u1[,2],col = "grey")

# Add further taildep functions from the spcopula package
rtPair <- 1-as.matrix(rankTransform(simTriHydro[,c(1,3)]))
tdfEmp <- empTailDepFun(rtPair)
plot(tdfEmp,ylim=c(0,1), ylab="tail index", xlab="u")
abline(v=0.5, col="grey")

gaussCop <- fitCopula(normalCopula(0), rtPair)@copula
tdfGauss <- tailDepFun(gaussCop)
curve(tdfGauss, add=T,col="green",n=500)

gumbelCop <- fitCopula(gumbelCopula(2),rtPair)@copula
tdfGumbel <- tailDepFun(gumbelCop)
curve(tdfGumbel,add=T, col="blue",n=500)

frankCop <- fitCopula(frankCopula(), rtPair)@copula
tdffrank <- tailDepFun(frankCop)
curve(tdffrank, add=T,col="red",n=500)

clayCop <- fitCopula(claytonCopula(), rtPair)@copula
tdfclay <- tailDepFun(clayCop)
curve(tdfclay, add=T,col="purple",n=500)

legend("bottom",
       c("empirical", "Gaussian", "Gumbel", "Frank", "Clayton"),
       col=c("black", "green", "blue", "red", "purple"), lwd = 2, lty=1, y.intersp = 0.45, seg.len = 1)


# Fit marginal distributions

peakGev <- fgev(simTriHydro$Peak)
volumeGev <- fgev(simTriHydro$Volume)
durationGev <- fgev(simTriHydro$Duration)

dpeak <- function(x)
  dgev(x, peakGev$estimate[1], peakGev$estimate[2], peakGev$estimate[3])
# dexp(x, peakExp$estimate)

hist(simTriHydro$Peak, freq=F, n=20, ylim=c(0,0.01), main = "Histogram of Peak")
curve(dpeak, add=T, col="red")

ppeak <- function(x)
  pgev(x, peakGev$estimate[1], peakGev$estimate[2], peakGev$estimate[3])
# pexp(x, peakExp$estimate)

plot(ecdf(simTriHydro$Peak),main = "Empirical cumulative distribution function: Peak" , lwd =2)
curve(ppeak, add=T, col="red")

dvolume <- function(x)
  dexp(x, volumeExp$estimate)

hist(simTriHydro$Volume, freq=F, n=20, ylim=c(0,1e-6), main = "Histogram of Volume")
curve(dvolume, add=T, col="red")

pvolume <- function(x)
  pexp(x, volumeExp$estimate)

plot(ecdf(simTriHydro$Volume), main = "Empirical cumulative distribution function:Volume", lwd = 2)
curve(pvolume, add=T, col="red")

dduration <- function(x)
  dgev(x, durationGev$estimate[1], durationGev$estimate[2], durationGev$estimate[3]-0.15)

hist(simTriHydro$Duration, freq=F, n=20, main = "Histogram of Duration")#, ylim=c(0,1e-6))

curve(dduration, add=T, col="red")

pduration <- function(x)
  pgev(x, durationGev$estimate[1], durationGev$estimate[2], durationGev$estimate[3]-0.15)

plot(ecdf(simTriHydro$Duration), main = "Empirical cumulative distribution function:Volume", lwd = 2)
curve(pduration, add=T, col="red")

# Selects an appropriate bivariate copula family for peak & volume distribution

BiCopSelect(hydro[,1], hydro[,2])
RVC <- RVineStructureSelect(hydro)
RVineLogLik(hydro, RVC)$loglik

RVineTreePlot(hydro, RVC, edge.labels= c("family", "theotau"), tree = 1:3)

