# Functions for use with Faraway's book
#
"maxadjr" <-
  function (l,best=3)
# Display the best (3) models from a leaps() object
{
  i <- rev(order(l$a))
  nopreds <- max(l$size)-1
  labels <- apply(l$which,1,function(x) paste(as.character((1:nopreds)[x]),collapse=","))
  m <- round(l$a[i[1:best]],3)
  names(m) <- labels[i[1:best]]
#  m <- cbind(round(l$a[i[1:best]],3),labels[i[1:best]])
#  dimnames(m) <- list(NULL,c("Adj R^2","Model"))
  m
}
"qqnorml" <-
  function(y,main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles",...)
# labeled Q-Q plot
  {
    n <- length(y)
    u <- qnorm((1:n)/(n+1))
    i <- order(y)
    plot(u,y[i],xlab=xlab,ylab=ylab,main=main,type="n")
    text(u,y[i],as.character(1:n)[i])
  }
"Cpplot" <-
  function (cp)
# Construct a Cp plot
{
  p <- max(cp$size)
  i <- (cp$Cp < (p+1.5))
  plot(cp$size[i],cp$Cp[i],xlab="p",ylab="Cp",type="n")
  labels <- apply(cp$which,1,function(x) paste(as.character((1:(p-1))[x]),collapse=""))
  text(cp$size[i],cp$Cp[i],labels[i])
  abline(0,1)
}
vif <- function(object)
UseMethod("vif")

vif.default <- function(object) {
  if(!is.data.frame(object) & !is.matrix(object)) stop("Not matrix or data frame")
  if(is.data.frame(object)) object <- as.matrix(object)
  ncols <- dim(object)[2]
  v <- numeric(ncols)
  names(v) <- dimnames(object)[[2]]
  for(i in 1:ncols) v[i] <- 1/(1-summary(lm(object[,i]~object[,-i]))$r.squared)
  v
}

# function from Bill Venables post on R-digest
vif.lm <- function(object) {
  V <- summary(object)$cov.unscaled
  Vi <- crossprod(model.matrix(object))
        nam <- names(coef(object))
  if(k <- match("(Intercept)", nam, nomatch = FALSE)) {
                v1 <- diag(V)[-k]
                v2 <- (diag(Vi)[-k] - Vi[k, -k]^2/Vi[k,k])
                nam <- nam[-k]
        } else {
                v1 <- diag(V)
                v2 <- diag(Vi)
                warning("No intercept term detected.  Results may surprise.")
        }
        structure(v1*v2, names = nam)
}

prplot <- function(g,i)
{
# Partial residuals plot for predictor i
  xl <- attributes(g$terms)$term.labels[i]
  yl <- paste("beta*",xl,"+res",sep="")
  x <- model.matrix(g)[,i+1]
  plot(x,g$coeff[i+1]*x+g$res,xlab=xl,ylab=yl)
  abline(0,g$coeff[i+1])
  invisible()
}
"halfnorm" <-
function (x, nlab = 2, labs = as.character(1:length(x)), ylab = "Sorted Data",
            ...)
{
  x <- abs(x)
  labord <- order(x)
  x <- sort(x)
  i <- order(x)
  n <- length(x)
  ui <- qnorm((n + 1:n)/(2 * n + 1))
  plot(ui, x[i], xlab = "Half-normal quantiles", ylab = ylab, ylim=c(0,max(x)),
       type = "n", ...)
  if(nlab < n)
    points(ui[1:(n - nlab)], x[i][1:(n - nlab)])
  text(ui[(n - nlab + 1):n], x[i][(n - nlab + 1):n], labs[labord][(n -
                                                              nlab + 1):n])
}
# logit and inverse logit
logit <- function(x){
  if(any(omit <- (is.na(x) | x <=0 | x >= 1))){
    lv <- x
    lv[omit] <- NA
    if(any(!omit))
      lv[!omit] <- Recall(x[!omit])
    return(lv)
  }
  log(x/(1-x))
}
ilogit <- function(x){
  if(any(omit <- is.na(x))){
    lv <- x
    lv[omit] <- NA
    if(any(!omit))
      lv[!omit] <- Recall(x[!omit])
    return(lv)
  }
  exp(x)/(1 + exp(x))
}

# Essential regression summary (idea from Gelman and Hill)
sumary <- function(object){
  digits <- options()$digits
  summ <- summary (object)
  sigma.hat <- summ$sigma
  r.squared <- summ$r.squared
  coef <- summ$coef[,,drop=FALSE]
  n <- summ$df[1] + summ$df[2]
  p <- summ$df[1]
  if (nsingular <- summ$df[3] - summ$df[1]) cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
  printCoefmat(coef,signif.stars=FALSE)
  cat("\n")
  cat (paste ("n = ", n, ", p = ", p,
    ", Residual SE = ", format(round(sigma.hat, digits-2),nsmall=digits-2),
    ", R-Squared = ", format(round(r.squared, 2)), "\n", sep=""))
  invisible(summ)
}
