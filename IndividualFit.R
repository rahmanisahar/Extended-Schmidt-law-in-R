# Modified Script to perform Baysian linear regression using three parameters with measurment error i.e. The Extended Schmidt law
# Modified Nov 2013 by Sahar Rahmani & Dr. Pauline Barmby 
# Script to perform Bayesian mutiple linear regression, with measurement error.
# Orginal Script was from Rahul Shetty (2012) fot fitting Schmidt-Kennicutt law  
# Uses other scripts from:
# Kruschke, J. K. (2011). Doing Bayesian Data Analysis

graphics.off()
if (1) {
#rm(list=ls(all=TRUE))
if ( .Platform$OS.type != "windows" ) {
  windows <- function( ... ) X11( ... )
}
require(rjags)         
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
    for( i in 1 : Ndata ) {

        x[i] ~ dnorm( mux[i] , 1.0/(xerr[i]^2) )
        mux[i] ~ dnorm( xbar, taux )

        x2[i] ~ dnorm( mux2[i] , 1.0/(x2err[i]^2) )
        mux2[i] ~ dnorm( x2bar, taux2 )

        y[i] ~ dnorm( muy[i] , 1.0/(yerr[i]^2) )

        muy[i] ~ dnorm( mu[i] , tauy )

        mu[i] <- beta0 + beta1 * mux[i] + beta2 * mux2[i] + muscat[i]
        muscat[i] ~ dnorm( 0 , tauscat )

    }

    beta0 ~ dnorm( 0 , 1.0E-12 )
    beta1 ~ dnorm( 0 , 1.0E-12 )
    beta2 ~ dnorm( 0 , 1.0E-12 )
    xbar ~ dnorm( 0 , 1.0E-12 )
    x2bar ~ dnorm( 0 , 1.0E-12 )

    tauy ~ dgamma( 0.001 , 0.001 )
    taux ~ dgamma( 0.001 , 0.001 )
    taux2 ~ dgamma( 0.001 , 0.001 )
    tauscat ~ dgamma( 0.001 , 0.001 )

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

alldat=read.csv("M31data.csv")

#SurfD: Gas surface Density ; StarD: Stellar mass surface density
# Columns are "subjID" , "SurfD", "errSurfD", "StarD", "errStarD", "SFR", "errSFR"
colnames(alldat)=c("subjID" , "SurfD" , "errSurfD", "StarD", "errStarD","SFR" , "errSFR")

# set strings for plotting labels and file names:  
yColName = "SFR" ; yPlotLab = "log SFR"
xColName = "SurfD" ; xPlotLab = "Surf Density"
x2ColName = "StarfD" ; x2PlotLab = "Stellar mass surface density"
xerrColName="errSurfD"
x2errCollname="errStarD"
yerrColName="errSFR"
subjColName = "subjID" ; subjPlotLab = "Galaxy"
xName=expression(Sigma[CO])
x2Name=expression(Sigma[Stellar mass])
yName=expression(Sigma[SFR])

fileNameRoot = "m31Fit"

# Extract data info:
Ndata = NROW(alldat)
# To make sure that subj has same order of subjects as dataMat, must use
# levels=unique() argument in factor statement...

x=alldat[,2]
x2=alldat[,4]
y=alldat[,6]
xerr = alldat[,3]
x2err = alldat[,5]
yerr = alldat[,7]

nSubj=length(x)
pdf("first_plot.pdf")
scatterplot3d(x,x2,y)   # needs special package from http://cran.r-project.org/web/packages/scatterplot3d/index.html
dev.off()
#------------------------------------------------------------------------------
# First inspect results from Multiple Linear Regression
# result: lmfit$coef[1] + lmfit$coef[2]*x + lmfit$coef[3]*x2 
lmfit <- lm(y~x+x2) # from http://www.statmethods.net/stats/regression.html


#show(summary(lmfit))
dev.copy(png)
#------------------------------------------------------------------------------

standardize=FALSE
#standardize=TRUE

# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.

if (standardize) {
    xM = mean( x ) ; xSD = sd( x )
    x2M = mean( x2 ) ; x2SD = sd ( x2 )
    yM = mean( y ) ; ySD = sd( y )
    zx = ( x - xM ) / xSD
    zx2 = ( x2 - x2M ) / x2SD
    zy = ( y - yM ) / ySD

    xerr=xerr*zx
    x2err=x2err*zx2
    yerr=yerr*zy
} else {
    zx=x
    zx2=x2
    zy=y
}

minx=min(x)
maxx=max(x)
minx2=min(x2)
maxx2=max(x2)
# Specify data, as a list.
dataList = list(
  x = zx ,
  x2 = zx2 ,
  y = zy ,
  Ndata = nSubj,
  xerr = xerr,
  x2err = x2err,
  yerr = yerr
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
r = lmfit$coef[2]*xSD/ySD # This is calculated here http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/MultiLinRegressInterJags.R as I unedrstand
r2 = lmfit$coef[3]*x2SD/ySD
initsList = list(
    beta0 = 0 ,    # because data are standardized
    beta1 = r ,        # because data are standardized
    beta2 = r2 ,       # beacuase data are standardized 
    tauy = (length(y)*sd(y)^2)/sum(lmfit$res^2),  # this 1 / sqrt((1-r^2)^2 + (1-r2^2)^2)? or this (length(y)*ySD^2)/sum(lmfit$res^2)
    taux = 1 / (1-r^2),  # because data are standardized
    taux2 = 1 / (1-r2^2),  # because data are standardized
    tauscat = 10.0,
    xbar=mean(x)
    x2bar=mean(x2)
)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# RUN THE CHAINS

if (1) {
parameters = c("beta0" , "beta1" , "beta2", "tauy", "taux", "taux2","tauscat", "xbar","x2bar")  # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters ,
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices:
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

checkConvergence = F
show( summary( codaSamples ) )
if ( checkConvergence ) {
  #windows()
  plot( codaSamples , ask=T )
  #windows()
  autocorr.plot( codaSamples , ask=T )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract chain values:
z0 = matrix( mcmcChain[, "beta0" ] )
z1 = matrix( mcmcChain[, "beta1" ] )
z2 = matrix( mcmcChain[, "beta2" ] )
zTauy = matrix( mcmcChain[, "tauy" ] )
zTaux = matrix( mcmcChain[, "taux" ] )
zTaux2 = matrix( mcmcChain[, "taux2" ] )
zTauscat = matrix( mcmcChain[, "tauscat" ] )
zxbar = matrix( mcmcChain[, "xbar" ] )
zx2bar = matrix( mcmcChain[, "x2bar" ] )

} # if on line 2

# Convert to original scale:
if (standardize) {
    b1 = z1 * ySD / xSD # cause there is no correlation between x and x1
    b2 = z2 * ySD / x2SD
    b0 = ( z0 * ySD + yM - z1 * ySD * xM / xSD - z2 * ySD * x2M / x2SD )

    sigmascat = 1/sqrt(zTauscat * ySD)  # Convert precision to SD
    sigmay = 1 / sqrt( zTauy )

} else {

    b1 = z1
    b2 = z2
    b0 = z0
    xbar = zxbar
    x2bar = zx2bar

# Convert precision to SD
    sigmay = 1 / sqrt( zTauy )
    sigmascat = 1 / sqrt( zTauscat ) # Convert precision to SD

    sigmax = 1 / sqrt( zTaux ) # Convert precision to SD, but sigmax not really needed
    sigmax2 = 1 / sqrt( zTaux2 ) # Convert precision to SD, but sigmax not really needed


}

#windows()
pdf("first_hist_b0.pdf")
hist(b0, prob=TRUE)
dev.off()
#dev.copy2pdf(file=pase("first_hist_b0.pdf"))
#windows()
pdf("first_hist_b1.pdf")
hist(b1, prob=TRUE)
dev.off()

pdf("first_hist_b2.pdf")
hist(b2, prob=TRUE)
dev.off()
# Plots and Posterior Prediction? (set showplot=TRUE)
makeplots=TRUE

if (makeplots) {
# Posterior prediction:
# Specify x values for which predicted y's are needed:
#xPostPred = seq(min(x),max(x),length=20)


numPostPred=20
xPostPred = seq(min(x),max(x),length=numPostPred)
postSampSize = length(b1)
xHDIlim = matrix( 0 , nrow=numPostPred , ncol=2 )

# Define matrix for recording posterior predicted y values at each x value.
# One row per x value, with each row holding random predicted y values.

yPostPred = matrix( 0 , nrow=numPostPred , ncol=postSampSize )
expression(sigma[v] (km/s))
# Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim = matrix( 0 , nrow=numPostPred , ncol=2 )

# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.

for ( chainIdx in 1:postSampSize ) {
    yPostPred[,chainIdx] = rnorm( numPostPred ,
                           mean = b0[chainIdx] + b1[chainIdx] * xPostPred + sigmascat[chainIdx],
                           sd = rep( sigmay[chainIdx] , numPostPred ))
}

source("HDIofMCMC.R")
for ( xIdx in 1:numPostPred ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
}

# Display believable beta0 and b1 values
#quartz()
pdf("M31fitSlopeIntercept.pdf")
thinIdx = seq(1,length(b0),length=700)
plot( z1[thinIdx] , z0[thinIdx] , cex.lab=1.75 ,
      ylab="Standardized Intercept" , xlab="X" )
plot( b1[thinIdx] , b0[thinIdx] , cex.lab=1.75 ,
      ylab="Intercept (ht when wt=0)" , xlab="X" )
#dev.copy2pdf( file = paste(fileNameRoot,"SlopeIntercept.pdf",sep="") )
dev.off()
# Display the posterior of the b1:
source("plotPost.R")
#quartz()
pdf("M31fitPostSlope.pdf")
histInfo = plotPost( z1 , xlab="Standardized slope" , compVal=0.0 ,
                     breaks=30  )
histInfo = plotPost( b1 , xlab="Slope" , compVal=0.0 ,
                     breaks=30  )
#dev.copy2pdf( file = paste(fileNameRoot,"PostSlope.pdf",sep="") )
dev.off()
#quartz()
pdf("M31fitPostIntercept.pdf")
histInfo = plotPost( z0 , xlab="Standardized Int" , compVal=0.0 ,
                     breaks=30  )
histInfo = plotPost( b0 , xlab="Intercept" , compVal=0.0 ,
                     breaks=30  )
#dev.copy2pdf( file = paste(fileNameRoot,"PostIntercept.pdf",sep="") )
dev.off()
# Display data with believable regression lines and posterior predictions.
#windows()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
xRang = max(x)-min(x)
yRang = max(y)-min(y)
limMult = 0.25
xLim= c( min(x)-limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y)-limMult*yRang , max(y)+limMult*yRang )
pdf("M31fitDataLines.pdf")
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="X" , ylab="Y", cex.lab=1.5 ,
      main="Data with credible regression lines" , cex.main=1.33  )
# Superimpose a smattering of believable regression lines:
for ( i in seq(from=1,to=length(b0),length=50) ) {
    abline( b0[i] , b1[i] , col="grey" )
}
#dev.copy2pdf( file = paste(fileNameRoot,"DataLines.pdf",sep="") )
dev.off()
# Display data with HDIs of posterior predictions.
#quartz()
# Plot data values:
pdf("M31fitDataPred.pdf")
yLim= c( min(yHDIlim) , max(yHDIlim) )
plot( x , y , cex=1.5 , lwd=2 , col="black" ,
      xlab="X" , ylab="Y", cex.lab=1.5 ,
      main="Data with 95% HDI & Mean of Posterior Predictions" , cex.main=1.33  )

# Superimpose posterior predicted 95% HDIs:
for ( i in 1:numPostPred ) {
    segments( xPostPred, yHDIlim[,1] , xPostPred, yHDIlim[,2] , lwd=3, col="grey" )
    points( xPostPred , rowMeans( yPostPred ) , pch="+" , cex=2 , col="grey" )
}
#dev.copy2pdf( file = paste(fileNameRoot,"DataPred.pdf",sep="") )
dev.off()

#quartz()
pdf("M31fitySigma.pdf")
histInfo = plotPost( sigmay , xlab=expression(sigma[Y]) , compVal=0.0 ,
                     breaks=30  )
#dev.copy2pdf( file = paste(fileNameRoot,"ySigma.pdf",sep="") )
dev.off()
}
}
#------------------------------------------------------------------------------
