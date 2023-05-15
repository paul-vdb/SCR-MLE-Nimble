
## This function compiles using compileNimble directly, but not in the model compile...
getLogESA <- nimbleFunction(
	# setup = function(){}, 
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
	# methods = list(
      # derivsRun = function(sigma = double(), lambda = double(), d2mask = double(2), 
		# nmask = integer(), ntrap = integer(), Time = double(), area = double()) {
		# wrt <- 1:2
		# return(derivs(run(sigma, lambda, d2mask, nmask, ntrap,
			# Time, area), wrt = wrt, order = 0:2))
		# returnType(ADNimbleList())
    # }
  #),
	buildDerivs = list(run = list(ignore = c('d2mask', 'nmask', 'ntrap', 'Time', 'area', 'i', 'j')))
)

# Rget <- getLogESA()
# Rget$run(0.5, 0.2, d2mask, nmask, J, StudyPeriod, A)
# Rget$derivsRun(0.5, 0.2, d2mask, nmask, J, StudyPeriod, A)
# cgetESA <- compileNimble(Rget)
# cgetESA$run(0.5, 0.2, d2mask, nrow(d2mask), ncol(d2mask), StudyPeriod, A)
# cgetESA$derivsRun(0.5, 0.2, d2mask, nrow(d2mask), ncol(d2mask), StudyPeriod, A)

dCount_cond <- nimbleFunction(
    run = function(x = integer(0), D = double(), sigma = double(), lambda = double(), 
			d2mask = double(2), Time = double(), area = double(), nmask = integer(), ntrap = integer(), 
			log = integer(0, default = 0)) {
		ESA <- 0.0
		sigma2 <- sigma*sigma*2
		for( i in 1:nmask ){
		  Hk <- 0.0
		  for( j in 1:ntrap) {
			Hk <- Hk + lambda*Time*exp(-d2mask[i,j]/sigma2)
		  }
		  ESA <- ESA + (1-exp(-Hk))*area
		} 
		
        ll <- -D*ESA + x*log(D)
        returnType(double(0))
		if(log) return(ll) else return(exp(ll))
		},
	buildDerivs = list(run = list(ignore = c('d2mask', 'nmask', 'ntrap', 'Time', 'area', 'i', 'j')))
)

