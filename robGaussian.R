## Id: robGaussian.R
## Author: Felipe Osorio
## Last update: 11-02-2017.

robGaussian <-
function(x, data, subset, na.action, qpar = 1, maxiter = 200, tol = 1e-6)
{ # robust estimation of center and Scatter from a multivariate Gaussian distribution
  Call <- match.call()
  if (missing(x))
    stop("'x' is not supplied")
  if (inherits(x, "formula")) {
    mt <- terms(x, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "x"] <- "formula"
    mf$qpar <- mf$maxiter <- mf$tol <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  }
  else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("gaussianMLq applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  iter <- 0
  eps <- 1e-2

  ## initial estimates
  wts <- rep(1, n)
  obj <- cov.wt(z, wt = wts, center = TRUE, method = "ML")
  old.center <- obj$center
  old.Sigma <- obj$cov
  old.theta <- c(old.center, old.Sigma[lower.tri(old.Sigma, diag = TRUE)])

  ## iterative loop
  now <- proc.time()
  repeat {
    # weigths computation
    distances <- mahalanobis(z, old.center, old.Sigma)
    val <- det(old.Sigma)^(-.5 * (1 - qpar))
    val <- val * (2 * pi)^(-p * (1 - qpar))
    wts <- val * exp(-.5 * (1 - qpar) * distances)

    # updating estimates
    obj <- cov.wt(z, wt = wts, center = TRUE, method = "ML")
    center <- obj$center
    Sigma <- obj$cov

    # eval convergence
    iter <- iter + 1
    theta <- c(center, Sigma[lower.tri(Sigma, diag = TRUE)])
    conv <- sum((theta - old.theta)^2) / (sum(theta^2) + eps)
    if (conv < tol) # successful completion
      break
    if (iter >= maxiter) # maximum number of iterations exceeded
      break

    # save estimates
    old.center <- center
    old.Sigma <- Sigma
    old.theta <- theta
  }
  speed <- proc.time() - now

  ## computing the Lq-likelihood
  distances <- mahalanobis(z, center, Sigma)
  if (qpar != 1) {
    val <- -.5 * (1 - qpar) * (p * log(2 * pi) + log(det(Sigma)))
    LqLik <- exp(val) * sum(exp(-.5 * (1 - qpar) * distances))
    LqLik <- (LqLik - n) / (1 - qpar)
  } else {
    val <- -.5 * p * log(2 * pi) - .5 * log(det(Sigma))
    LqLik <- n * val - .5 * sum(distances)
  }

  ## creating the output object
  out <- list(call = Call, x = z, dims = dz, center = center, Sigma = qpar * Sigma,
              qpar = qpar, weights = wts, distances = distances, LqLik = LqLik,
              iterations = iter, speed = speed)
  class(out) <- "robGaussian"
  out
}

print.robGaussian <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Call:\n")
  dput(x$call, control = NULL)
  cat("\nCenter:\n ")
  print(format(round(x$center, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  if (x$dims[2] <= 5)
    print.symmetric(x$Sigma, digits = digits)
  else {
    print.symmetric(x$Sigma[1:5,1:5], digits = digits)
    cat("...")
  }
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs)
  p <- x$dims[2]
  cat("\nLq-likelihood:", format(x$LqLik), "on", p * (p + 3) / 2, "degrees of freedom\n")
  invisible(x)
}
