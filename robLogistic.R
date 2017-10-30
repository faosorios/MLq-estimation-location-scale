
robLogistic <-
function(formula, data, subset, na.action, qpar = 1, control = glm.control(),
  model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...)
{
  # support functions
  linkfun <- function(mu) log(mu / (1 - mu))
  linkinv <- function(eta) {
    thresh <- -log(.Machine$double.eps)
    eta <- pmin(thresh, pmax(eta, -thresh))
    exp(eta)/(1 + exp(eta))
  }
  kernel <- function(eta) log(1 + exp(eta))
  variance <- function(mu) mu * (1 - mu)

  # extract model matrix and response variable
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$qpar <- mf$control <- mf$model <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$... <- NULL
  # mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)
  ynames <- names(y)
  xnames <- dimnames(x)[[2]]
  dims <- dim(x)
  nobs <- dims[1]

  # additional control parameters
  eps <- 1e-2

  # checking arguments
  if (any(y < 0 | y > 1))
    stop("y values must be 0 <= y <= 1")
  if (is.null(control))
    control <- glm.control()

  now <- proc.time()
  ## initial estimates
  z   <- log((y + 0.5) / (1 - y + .5))
  fit <- lsfit(x, z, intercept = FALSE)
  oldcoef <- fit$coefficients
  eta <- c(x %*% oldcoef)

  # IRLS iterations
  iter <- 0
  repeat {
    mu   <- linkinv(eta)
    rob.wts <- exp((1 - qpar) * (y * eta - kernel(eta)))
    k    <- exp(qpar * kernel(eta) - kernel(qpar * eta))
    wts  <- k * variance(mu)
    z    <- eta + rob.wts * (y - mu) / wts
    fit  <- lsfit(x, z, wt = wts, intercept = FALSE)
    coef <- fit$coef

    iter <- iter + 1
    diff <- coef - oldcoef
    conv <- sum(diff^2) / (sum(oldcoef^2) + eps)

    if (conv < control$epsilon) break
    if (iter >= control$maxit)  break

    eta <- c(x %*% coef)
    oldcoef <- coef
  }

  ## calculate df
  nulldf <- nobs - sum(wts == 0)
  rank   <- fit$qr$rank
  resdf  <- nobs - rank

  speed <- proc.time() - now
  ## creating the output object
  mu <- linkinv(eta)
  out <- list(call = Call,
              dims = dims,
              qpar = qpar,
              coefficients = qpar * coef,
              fitted.values = mu,
              linear.predictor = eta,
              residuals = y - mu,
              numIter = iter,
              control = control,
              rob.weights = rob.wts,
              rank = rank,
              df.null = nulldf,
              df.residual = resdf,
              speed = speed)
  class(out) <- "logistic"
  out
}

print.logistic <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  dput(x$call, control = NULL)
  cat("Converged in", x$numIter, "iterations\n\n")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  }
  else
    cat("No coefficients\n\n")
  cat("\nDistortion parameter:", x$qpar)
  cat("\nDegrees of Freedom:", x$df.null, "Total;", x$df.residual, "Residual\n")
  invisible(x)
}
