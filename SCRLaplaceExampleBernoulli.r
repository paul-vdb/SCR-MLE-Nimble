library(nimble)
nimbleOptions(buildModelDerivs = TRUE)

source("functions.r")
load('test-data.Rdata')

capt <- test.data$bin.capt
traps <- test.data$trapsT
mask <- as.matrix(test.data$mask)
area <- attr(test.data$mask, "area")

d2mask <- as.matrix(t(apply(mask, 1, FUN = function(x){(x[1]-traps[,1])^2 + (x[2]-traps[,2])^2})))

## Conditional SCR model
SCR_model <- nimbleCode({
	g0 ~ dbeta(1, 1) # Detection rate at distance 0
    sigma ~ dgamma(1, 1)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	D <- K0/logESA

	## Doesn't Work
	ESA <- getESABernoulli(sigma = sigma, g0 = g0, 
		d2mask = d2mask[1:nmask, 1:J], nmask = nmask, 
		ntrap = J, area = A)
	logESA <- log(ESA)
	
    for(k in 1:K0) {
        X[k, 1] ~ dnorm(0,1)
        X[k, 2] ~ dnorm(0,1)
		## Distances
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
        ## hazard rate
		pkj[k, 1:J] <- g0*exp(-d2[k,1:J]*tau2)
		
		zeros[k] ~ dpois(logESA)	## 1/ESA contribution to likelihood
		## Bernoulli Detections
		for(j in 1:J) y[k, j] ~ dbinom(prob=pkj[k, j], size = 1)
		## Zero-truncation part of the likelihood cancels with the conditional 
		## activity centre distribution, leaving just 1/ESA.
	}
})

K <- nrow(capt)
data <- list(y = capt, zeros = rep(0, K))	#, K I added it as data instead of the prior and that allows me to change it but not influence the Log Likelihood.
constants <- list(traps = traps, d2mask = d2mask, J = nrow(traps), nmask = nrow(mask), 
	A = area, area = area*nrow(mask), K0 = nrow(capt))
inits <- list(sigma = 0.5, g0 = 0.6) ##D = 50,  add D back when I get the dist function working.
Rmodel <- nimbleModel(SCR_model, data=data, constants=constants, inits = inits, buildDerivs = TRUE)
Cmodel <- compileNimble(Rmodel)
scr_laplace <-  buildLaplace(Rmodel, paramNodes = c('sigma', 'g0'), randomEffectsNodes = 'X')
Cscr_laplace <- compileNimble(scr_laplace, project = Rmodel)
mle <- Cscr_laplace$findMLE()

## HT-like estimator for Density:
Cmodel[['sigma']] <- mle$par[1]
Cmodel[['lambda']] <- mle$par[2]
Cmodel$calculate()
Dhat <- K/Cmodel$ESA 	# Density Estimate
Nhat <- Dhat*area