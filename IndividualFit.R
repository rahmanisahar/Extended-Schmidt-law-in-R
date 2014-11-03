# Script to perform Bayesian linear regression, with measurement error.
# Created May, 2012 - Rahul Shetty
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

        y[i] ~ dnorm( muy[i] , 1.0/(yerr[i]^2) )

        muy[i] ~ dnorm( mu[i] , tauy )

        mu[i] <- beta0 + beta1 * mux[i] + muscat[i]
        muscat[i] ~ dnorm( 0 , tauscat )

    }

    beta0 ~ dnorm( 0 , 1.0E-12 )
    beta1 ~ dnorm( 0 , 1.0E-12 )
    xbar ~ dnorm( 0 , 1.0E-12 )

    tauy ~ dgamma( 0.001 , 0.001 )
    taux ~ dgamma( 0.001 , 0.001 )
    tauscat ~ dgamma( 0.001 , 0.001 )

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

alldat=read.csv("M51data.csv")

# Columns are "subjID" , "SurfD" , "errSurfD", "SFR" , "errSFR"
colnames(alldat)=c("subjID" , "SurfD" , "errSurfD", "SFR" , "errSFR")

# set strings for plotting labels and file names:  
yColName = "SFR" ; yPlotLab = "log SFR"
xColName = "SurfD" ; xPlotLab = "Surf Density"
xerrColName="errSurfD"
yerrColName="errSFR"
subjColName = "subjID" ; subjPlotLab = "Galaxy"
xName=expression(Sigma[CO])
yName=expression(Sigma[SFR])

fileNameRoot = "M51fit"

# Extract data info:
Ndata = NROW(alldat)
# To make sure that subj has same order of subjects as dataMat, must use
# levels=unique() argument in factor statement...

x=alldat[,2]
y=alldat[,4]
xerr = alldat[,3]
yerr = alldat[,5]

nSubj=length(x)

plot(x,y)

#------------------------------------------------------------------------------
# First inspect results from Simple Linear Regression

lmfit <- lm(y~x)
show(summary(lmfit))
#------------------------------------------------------------------------------

standardize=FALSE
#standardize=TRUE

# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.

if (standardize) {
    xM = mean( x ) ; xSD = sd( x )
    yM = mean( y ) ; ySD = sd( y )
    zx = ( x - xM ) / xSD
    zy = ( y - yM ) / ySD

    xerr=xerr*zx
    yerr=yerr*zy
} else {
    zx=x
    zy=y
}

minx=min(x)
maxx=max(x)

# Specify data, as a list.
dataList = list(
  x = zx ,
  y = zy ,
  Ndata = nSubj,
  xerr = xerr,
  yerr = yerr
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

r = cor(x,y)
initsList = list(
    beta0 = 0 ,    # because data are standardized
    beta1 = r ,        # because data are standardized
    tauy = 1 / (1-r^2),  # because data are standardized
    taux = 1 / (1-r^2),  # because data are standardized
    tauscat = 10.0,
    xbar=mean(x)
)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# RUN THE CHAINS

if (1) {
parameters = c("beta0" , "beta1" , "tauy", "taux", "tauscat", "xbar")  # The parameter(s) to be monitored.
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
  windows()
  plot( codaSamples , ask=T )
  windows()
  autocorr.plot( codaSamples , ask=T )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract chain values:
z0 = mcmcChain[, "beta0" ]
z1 = mcmcChain[, "beta1" ]
zTauy = mcmcChain[, "tauy" ]
zTaux = mcmcChain[, "taux" ]
zTauscat = mcmcChain[, "tauscat" ]
zxbar = mcmcChain[, "xbar" ]

} # if on line 2

# Convert to original scale:
if (standardize) {
    b1 = z1 * ySD / xSD
    b0 = ( z0 * ySD + yM - z1 * ySD * xM / xSD )

    sigmascat = 1/sqrt(zTauscat * ySD)  # Convert precision to SD
    sigmay = 1 / sqrt( zTauy )

} else {

    b1 = z1
    b0 = z0
    xbar = zxbar

# Convert precision to SD
    sigmay = 1 / sqrt( zTauy )
    sigmascat = 1 / sqrt( zTauscat ) # Convert precision to SD

    sigmax = 1 / sqrt( zTaux ) # Convert precision to SD, but sigmax not really needed

}

windows()
hist(b0, prob=TRUE)
windows()
hist(b1, prob=TRUE)

# Plots and Posterior Prediction? (set showplot=TRUE)
makeplots=FALSE

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
quartz()
thinIdx = seq(1,length(b0),length=700)
plot( z1[thinIdx] , z0[thinIdx] , cex.lab=1.75 ,
      ylab="Standardized Intercept" , xlab="X" )
plot( b1[thinIdx] , b0[thinIdx] , cex.lab=1.75 ,
      ylab="Intercept (ht when wt=0)" , xlab="X" )
dev.copy2eps( file = paste(fileNameRoot,"SlopeIntercept.eps",sep="") )

# Display the posterior of the b1:
source("plotPost.R")
quartz()
histInfo = plotPost( z1 , xlab="Standardized slope" , compVal=0.0 ,
                     breaks=30  )
histInfo = plotPost( b1 , xlab="Slope" , compVal=0.0 ,
                     breaks=30  )
dev.copy2eps( file = paste(fileNameRoot,"PostSlope.eps",sep="") )

quartz()
histInfo = plotPost( z0 , xlab="Standardized Int" , compVal=0.0 ,
                     breaks=30  )
histInfo = plotPost( b0 , xlab="Intercept" , compVal=0.0 ,
                     breaks=30  )
dev.copy2eps( file = paste(fileNameRoot,"PostIntercept.eps",sep="") )

# Display data with believable regression lines and posterior predictions.
windows()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
xRang = max(x)-min(x)
yRang = max(y)-min(y)
limMult = 0.25
xLim= c( min(x)-limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y)-limMult*yRang , max(y)+limMult*yRang )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="X" , ylab="Y", cex.lab=1.5 ,
      main="Data with credible regression lines" , cex.main=1.33  )
# Superimpose a smattering of believable regression lines:
for ( i in seq(from=1,to=length(b0),length=50) ) {
    abline( b0[i] , b1[i] , col="grey" )
}
dev.copy2eps( file = paste(fileNameRoot,"DataLines.eps",sep="") )

# Display data with HDIs of posterior predictions.
quartz()
# Plot data values:
yLim= c( min(yHDIlim) , max(yHDIlim) )
plot( x , y , cex=1.5 , lwd=2 , col="black" ,
      xlab="X" , ylab="Y", cex.lab=1.5 ,
      main="Data with 95% HDI & Mean of Posterior Predictions" , cex.main=1.33  )

# Superimpose posterior predicted 95% HDIs:
for ( i in 1:numPostPred ) {
    segments( xPostPred, yHDIlim[,1] , xPostPred, yHDIlim[,2] , lwd=3, col="grey" )
    points( xPostPred , rowMeans( yPostPred ) , pch="+" , cex=2 , col="grey" )
}
dev.copy2eps( file = paste(fileNameRoot,"DataPred.eps",sep="") )


quartz()
histInfo = plotPost( sigmay , xlab=expression(sigma[Y]) , compVal=0.0 ,
                     breaks=30  )
dev.copy2eps( file = paste(fileNameRoot,"ySigma.eps",sep="") )

}
}
#------------------------------------------------------------------------------
