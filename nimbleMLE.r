library(nimble)
library(spatstat)
library(microbenchmark)
nimbleOptions(buildModelDerivs = TRUE)

# load("test-data.Rdata")
# str(test.data)
source('scr-ll.r') ## #evalulate Ben's likelihood that is slow to compare.

capt <- test.data$bin.capt
traps <- test.data$traps
mask <- as.matrix(test.data$mask)
area <- attr(test.data$mask, "area")

## Step 1: Just compare Nimble and Plain R.
nimESA <- nimbleFunction(
	setup = function(traps, mask, area){
		d2mask <- apply(traps, 1, FUN  = function(x){(x[1] - mask[,1])^2 + (x[2] - mask[,2])^2})
		ntrap <- nrow(traps)
		nmask <- nrow(mask)
	},
	run = function(sigma = double(), lambda = double()){
		sigma2 <- sigma*sigma*2
		pMask <- lambda*exp(-d2mask/sigma2)
		pAvoid = numeric(value = 1, length = nmask)
		for(j in 1:ntrap)
		{
			pAvoid <- pAvoid * (1 - pMask[,j])
		}
		return(area*sum(1-pAvoid))
		returnType(double())
	}
)

calcESA <- nimESA(traps, mask, area)
calcESA.c <- compileNimble(calcESA)

benchmark(calcESA(100, 5))


SCRBernoulli <- nimbleFunction(
  setup = function(capt = double(2), traps = double(2), 
					area = double(), mask = double(2)){
    n <- nrow(capt)
    n.traps <- nrow(traps)
    n.mask <- nrow(mask)
    a <- area

    mask.dists <- crossdist(mask[, 1], mask [, 2],
                            traps[, 1], traps[, 2])
    mask.dists.sq <- mask.dists^2
	mask.prob <- matrix(0, nrow = n.mask, ncol = n.traps)
	
	machinePrec <- .Machine$double.xmin
  },
  methods = list(
    calcESA = function(g0 = double(), sigma = double()){
	  ESA <- 0.0
	  sigma2 <- sigma*sigma*2
	  # mask.prob <- matrix(value = 0, nrow = n.mask, ncol = n.traps)
	  for( i in 1:n.mask )
	  {
	    pAvoid <- 1.0
	    for( j in 1:n.traps ) {
	      pDetect <- g0*exp(-mask.dists.sq[i,j]/sigma2)
		  mask.prob[i,j] <<- ADbreak(pDetect)
	      pAvoid <- pAvoid * (1 - pDetect)
		}
        ESA <- ESA + a*(1 - pAvoid)
	  }
	return(ESA)
	returnType(double())
    },
	gr_ESA = function(g0 = double(), sigma = double()){
  	  ans <- nimDerivs(calcESA(g0, sigma), wrt = 1:2, order = 1)
	  return(ans$jacobian[1,])
	  returnType(double(1))
	},
	negLogLikelihood = function(pars = double(1)){
      g0 <- 1/(1+exp(-pars[1]))
      sigma <- exp(pars[2])
      D <- exp(pars[3])
	  
	  ## Calc ESA and mask probabilities
	  ESA <- 0.0
	  sigma2 <- sigma*sigma*2
	  maskP <- nimMatrix(value = 0, nrow = n.mask, ncol = n.traps)
	  
	  for( i in 1:n.mask )
	  {
	    pAvoid <- 1.0
	    for( j in 1:n.traps ) {
	      # pDetect <- g0*exp(-mask.dists.sq[i,j]^2/sigma2)
		  maskP[i,j] <- g0*exp(-mask.dists.sq[i,j]/sigma2)
	      pAvoid <- pAvoid * (1 - maskP[i,j])
		}
        ESA <- ESA + a*(1 - pAvoid)
	  }

	  ll <- 0.0
	  for( k in 1:n ){
		ll.k <- 0.0
		for( i in 1:n.mask){
		  fcapt.i <- 0.0
	      for( j in 1:n.traps ){
	        fcapt.i <- fcapt.i + capt[k,j]*log(maskP[i,j] + machinePrec) + 
				(1-capt[k,j])*log(1-maskP[i,j] + machinePrec) 
	      }
		  ll.k <- ll.k + exp(fcapt.i)*a
		}
	    ll <- ll + log(ll.k + machinePrec)
	  }
	ll <- ll + n*log(D) - D*ESA
	return(-ll)
	returnType(double())
  }, 
  gr_negLogLikelihood = function(pars = double(1)){
	  ans <- nimDerivs(negLogLikelihood(pars), wrt = 1:3, order = 1)
	  return(ans$jacobian[1,])
	  returnType(double(1))
	}
  ),
  buildDerivs = list(calcESA = list(ignore = c('i', 'j')),
					negLogLikelihood = list(ignore = c('i', 'j', 'k')))
)

test <- SCRBernoulli(capt = capt, traps = traps, 
					area = area, mask = mask)
test$calcESA(0.5, 100)
test$negLogLikelihood(c(qlogis(0.5), log(100), log(20/area)))

# debug(test$negLogLikelihood)

cTest <- compileNimble(test)
cTest$negLogLikelihood(c(qlogis(0.5), log(100), log(20/area)))
cTest$gr_negLogLikelihood(c(qlogis(0.5), log(100), log(20/area)))
cTest$gr_ESA(0.5, 100)

time1 <- system.time(
	MLE <- optim(par = c(0, log(100), log(50)), 
		fn = cTest$negLogLikelihood, gr = cTest$gr_negLogLikelihood, 
		method = "BFGS")
)

## Fitting the model.
time2 <- system.time(
	fit <- optim(par = c(0, log(100), log(50)), scr.nll, capt = test.data$bin.capt, 
		traps = test.data$traps, mask = test.data$mask)
)	

MLE <- nlminb(start = c(0, log(100), log(50)), 
	objective = cTest$negLogLikelihood, gradient = cTest$gr_negLogLikelihood)

plogis(MLE$par[1])
exp(MLE$par[2])
exp(MLE$par[3])

exp(fit$par[3])
exp(fit$par[1])
plogis(fit$par[2])
