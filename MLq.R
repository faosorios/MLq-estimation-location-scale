MLq <-
function(x, data, subset, na.action, qpar = 1, maxiter = 500, tol = 1e-6)
{
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
    stop("MLq applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  iter <- 0
  eps <- 1e-2

  ## initial estimates
  wts <- rep(1, n)
  obj <- cov.wt(z, wt = wts, center = TRUE, method = "ML")
  center <- obj$center
  Sigma <- obj$cov
  theta <- c(center, Sigma[lower.tri(Sigma, diag = TRUE)])

  ## iterative loop
  now <- proc.time()
  repeat {
    # weigths computation
    distances <- mahalanobis(z, center, Sigma)
    val <- det(Sigma)^(-.5 * (1 - qpar))
    val <- val * (2 * pi)^(-p * (1 - qpar))
    wts <- val * exp(-.5 * (1 - qpar) * distances)

    # updating estimates
    obj <- cov.wt(z, wt = wts, center = TRUE, method = "ML")
    new.center <- obj$center
    new.Sigma <- obj$cov

    # eval convergence
    iter <- iter + 1
    new.theta <- c(new.center, new.Sigma[lower.tri(new.Sigma, diag = TRUE)])
    conv <- sum((new.theta - theta)^2) / (sum(new.theta^2) + eps)
    if (conv < tol) # successful completion
      break
    if (iter >= maxiter) # maximum number of iterations exceeded
      break

    # save estimates
    center <- new.center
    Sigma <- new.Sigma
    theta <- new.theta
  }
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              x = z,
              dims = dz,
              center = new.center,
              Sigma = qpar * new.Sigma,
              qpar = qpar,
              weights = wts,
              distances = distances,
              iterations = iter,
              speed = speed)
  class(out) <- "MLq"
  out
}

print.MLq <-
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
  cat("\nNumber of Observations:", nobs, "\n")
  invisible(x)
}
