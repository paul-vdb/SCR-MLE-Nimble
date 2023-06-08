library(nimble)
nimbleOptions(buildModelDerivs = TRUE)
wombats <- read.csv('Wombats.csv')

model_code <- nimbleCode({
  # priors (will be ignored for MLE)
  alpha ~ dflat()
  beta ~ dflat()
  # random effects and data  
  for(i in 1:n) {
    # data
    y[i] ~ dpois(exp(alpha + beta*x[i]))
  }
})


model <- nimbleModel(model_code, 
	constants = list(n = nrow(wombats), x = wombats$nBurrows),
	data = list(y = wombats$nWombats), buildDerivs = TRUE) # Create nimble model object

compileNimble(model)

logLikelihood_nf <- nimbleFunction(
  setup = function(model, paramNodes) {
    # Determine nodes for calculating the log likelihood for parameters given by
    # paramNodes, ignoring any priors.
    calcNodes <- model$getDependencies(paramNodes, self = FALSE)
    # Set up the additional arguments for nimDerivs involving model$calculate
    derivsInfo <- makeModelDerivsInfo(model, paramNodes, calcNodes)
    updateNodes <- derivsInfo$updateNodes
    constantNodes <- derivsInfo$constantNodes
    # Create a parameter transformation between original and unconstrained
    # parameter spaces.
    transformer <- parameterTransform(model, paramNodes)
  },
  methods = list(
    neg_logLikelihood_p = function(p = double(1)) {
      # Put values in model and calculate negative log likelihood.
      values(model, paramNodes) <<- p
      return(-model$calculate(calcNodes))
      returnType(double())
    },
    neg_logLikelihood = function(ptrans = double(1)) {
      # Objective function for optim,
      # using transformed parameter space.
      p <- transformer$inverseTransform(ptrans)
      return(neg_logLikelihood_p(p))
      returnType(double())
    },
    gr_neg_logLikelihood = function(ptrans = double(1)) {
      # Gradient of neg log likelihood
      p <- transformer$inverseTransform(ptrans)
      d <- derivs(neg_logLikelihood_p(p), wrt = 1:length(p), order = 1,
                  model = model, updateNodes = updateNodes,
                  constantNodes = constantNodes)
      return(d$jacobian[1,])
      returnType(double(1))
    },
    transform = function(p = double(1)) {
      # Give user access to the transformation ...
      return(transformer$transform(p))
      returnType(double(1))
    },
    inverse = function(ptrans = double(1)) { # ... and its inverse.
      return(transformer$inverseTransform(ptrans))
      returnType(double(1))
    }
  ),
  buildDerivs = 'neg_logLikelihood_p'
)

wombatGLM <- logLikelihood_nf(model, c('alpha', 'beta'))
CwombatGLM <- compileNimble(wombatGLM, project = model)

MLE <- optim(c(0,0), # initial values
         fn = wombatGLM$neg_logLikelihood, # function to be minimized
         gr = wombatGLM$gr_neg_logLikelihood, # gradient of function to be minimized
         method = "BFGS") # optimization method to use
MLE$par