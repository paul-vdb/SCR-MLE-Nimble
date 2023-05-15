# remotes::install_github("nimble-dev/nimble", ref = "ADoak", subdir = "packages/nimble")
library(nimble)
nimbleOptions(buildModelDerivs = TRUE)

## Conditional SCR model
SCR_model <- nimbleCode({
	lambda ~ dgamma(1, 1) # Detection rate at distance 0
    sigma ~ dgamma(1, 1)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	# D ~ dgamma(1, 1)

	# logESA <- getLogESA(sigma = sigma, lambda = lambda, d2mask[1:nmask, 1:J], nmask, J, Time, A)
	
	for( i in 1:nmask ) {
		Hazk[i] <-  sum(lambda*Time*exp(-d2mask[i,1:J]*tau2))
		pDetect[i] <- (1-exp(-Hazk[i])
	}
	ESA <- sum(pDetect[1:nmask]))*A
	logESA <- log(ESA)
	
    for(k in 1:K0) {
        X[k, 1] ~ dnorm(0,1)
        X[k, 2] ~ dnorm(0,1)
		## Distances
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
        ## hazard rate
		Hkj[k, 1:J] <- lambda*exp(-d2[k,1:J]*tau2)*Time
		
		zeros[k] ~ dpois(logESA)	## 1/ESA contribution to likelihood
		## Count process
		for(j in 1:J) y[k, j] ~ dpois(Hkj[k,j])
		## Zero-truncation part of the likelihood cancels with the conditional 
		## activity centre distribution, leaving just 1/ESA.
	}
})

simSCR <- function(N = 50, sigma = 0.5, lambda = 0.5, StudyPeriod = 25, traps, xlim, ylim)
{
    locs <- cbind(x = runif(N, xlim[1], xlim[2]), 
                  y = runif(N, ylim[1], ylim[2]))
    J <- nrow(traps)
    capthist <- NULL
	ids <- NULL
    for(i in 1:N)
    {
        d2 <- (locs[i,1] - traps[,1])^2 + (locs[i,2] - traps[,2])^2
        Hkj <- lambda*exp(-d2/(2*sigma^2))*StudyPeriod
		nkj <- rpois(J, Hkj)
		if(sum(nkj) > 0) capthist <- rbind(capthist, nkj)
    }
	as.matrix(capthist)
}

## Start with an easy problem.
traps <- expand.grid(x = 1:5,y = 1:5)
N <- 50
sigma <- 1.5
lambda <- 0.4
StudyPeriod <- 25
xlim <- range(traps[,1]) + c(-3, 3)
ylim <- range(traps[,2]) + c(-3, 3)
mask <- expand.grid(x = seq(xlim[1], xlim[2], 0.25), y = seq(ylim[1], ylim[2], 0.25)  )
nmask <- nrow(mask)
J <- nrow(traps)
A <- 0.25^2
area <- nrow(mask) * A
d2mask <- as.matrix(t(apply(mask, 1, FUN = function(x){(x[1]-traps[,1])^2 + (x[2]-traps[,2])^2})))

y  <- simSCR(N, sigma, lambda, StudyPeriod, traps, xlim, ylim)

K <- nrow(y)
data <- list(y = y, zeros = rep(0, K))	#, K I added it as data instead of the prior and that allows me to change it but not influence the Log Likelihood.
constants <- list(traps = traps, d2mask = d2mask, J = nrow(traps), nmask = nrow(mask), 
	A = A, area = area, K0 = nrow(y), Time = StudyPeriod)
inits <- list(D = 50, sigma = 0.5, lambda = 1)	
Rmodel <- nimbleModel(SCR_model, data=data, constants=constants, inits = inits, buildDerivs = TRUE)
Cmodel <- compileNimble(Rmodel)
scr_laplace <-  buildLaplace(Rmodel, paramNodes = c('sigma', 'lambda'), randomEffectsNodes = 'X')
Cscr_laplace <- compileNimble(scr_laplace, project = Rmodel)
mle <- Cscr_laplace$findMLE()

## HT-like estimator for Density:
Cmodel[['sigma']] <- mle$par[1]
Cmodel[['lambda']] <- mle$par[2]
Cmodel$calculate()
Dhat <- K/Cmodel$ESA 	# Density Estimate
Nhat <- Dhat*area