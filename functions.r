
## This function compiles using compileNimble directly, but not in the model compile...
getESA <- nimbleFunction(
    run = function(sigma = double(), lambda = double(), d2mask = double(2), 
		nmask = integer(), ntrap = integer(), Time = double(), area = double()) { # type declarations
		ESA <- 0.0
		sigma2 <- sigma*sigma*2
		for( i in 1:nmask ){
		  Hk <- 0.0
		  for( j in 1:ntrap) {
			  d2mask_noAD <- ADbreak(d2mask[i, j])
			  Hk <- Hk + lambda*Time*exp(-d2mask_noAD/sigma2)
			}
			ESA <- ESA + (1-exp(-Hk))*area
		}
        return(ESA)
        returnType(double())  # return type declaration
    },  
	buildDerivs = list(run = list(ignore = c('nmask', 'ntrap', 'Time', 'area', 'i', 'j')))
)

## Testing this one. Same as getLogESA but has setup code and derivsRun so 
## that we can check that it compiles locally.
getLogESATest <- nimbleFunction(
	setup = function(){}, 
    run = function(sigma = double(), lambda = double(), d2mask = double(2), 
		nmask = integer(), ntrap = integer(), Time = double(), area = double()) { # type declarations
		ESA <- 0.0
		sigma2 <- sigma*sigma*2
		for( i in 1:nmask ){
		  Hk <- 0.0
		  for( j in 1:ntrap) {
			Hk <- Hk + lambda*Time*exp(-d2mask[i,j]/sigma2)
		  }
		  ESA <- ESA + (1-exp(-Hk))*area
		}
		lESA <- log(ESA)
        return(lESA)
        returnType(double())  # return type declaration
    },  
	methods = list(
      derivsRun = function(sigma = double(), lambda = double(), d2mask = double(2), 
		nmask = integer(), ntrap = integer(), Time = double(), area = double()) {
		return(derivs(run(sigma, lambda, d2mask, nmask, ntrap,
			Time, area),order = 0:2))
		returnType(ADNimbleList())
    }
  ),
	buildDerivs = list(run = list(ignore = c('d2mask', 'nmask', 'ntrap', 'Time', 'area', 'i', 'j')))
)

## Trying to get it to work as a distribution... Will revist once the other version is working.
dCount_cond <- nimbleFunction(
    run = function(x = integer(0), D = double(), sigma = double(), lambda = double(), 
			d2mask = double(2), Time = double(), area = double(), nmask = integer(), ntrap = integer(), 
			log = integer(0, default = 0)) {
		ESA <- 0.0
		sigma2 <- sigma*sigma*2
		for( i in 1:nmask ){
		  Hk <- 0.0
		  for( j in 1:ntrap) {
			  d2mask_noAD <- ADbreak(d2mask[i, j])
			  Hk <- Hk + lambda*Time*exp(-d2mask_noAD/sigma2)
			}
			ESA <- ESA + (1-exp(-Hk))*area
		}
		
        ll <- -D*ESA + x*log(D)
        returnType(double(0))
		if(log) return(ll) else return(exp(ll))
		},
	buildDerivs = list(run = list(ignore = c('nmask', 'ntrap', 'Time', 'area', 'i', 'j')))
)

