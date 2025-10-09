# ======================================================================
# R Code for the manuscript 
# "Nonparametric quantile regression captures regional variability
# and scaling deviations in Atlantic surfclam length–weight relationships"
#
# Author: Marta Sestelo - sestelo@uvigo.gal 
# ======================================================================




# --- Load required packages ---
library(sfsmisc)
library(mgcv)
library(ggplot2)

plotgg <- function(x, y, muhat, lci, uci, pcol = "grey80", col = "black",
                   ylim = NULL, main = NULL, CIlinecol = "transparent", 
                   CIcol = "black", points = TRUE, cex = 1, lty = 1, 
                   CIlty = 2, lwd = 0.5, CIlwd = 1, alpha = 0.2, 
                   ylab = NULL, xlab = NULL, add = FALSE, 
                   ggprevious = NULL, ...) {
  
  dat <- data.frame(
    x = as.numeric(x),
    y = y,
    muhat = muhat,
    lci = lci,
    uci = uci, 
    col = rep(col, length(x))
  )
  
  if (points == TRUE) {
    points_layer <- ggplot2::geom_point(
      data = dat,
      ggplot2::aes_string(x = "x", y = "y"),
      colour = pcol,
      size = cex
    )
  } else {
    points_layer <- NULL
  }
  
  ggplot2::ggplot() +
    points_layer + 
    ggplot2::geom_ribbon(
      data = dat,
      ggplot2::aes_string(x = "x", ymin = "lci", ymax = "uci"), 
      alpha = alpha,
      fill = CIcol,
      linetype = lty,
      size = CIlwd,
      col = CIlinecol
    ) +
    ggplot2::geom_line(
      data = dat,
      ggplot2::aes(x = x, y = muhat, colour = col, group = col), 
      size = lwd,
      linetype = lty,
      na.rm = TRUE
    ) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::ggtitle(main)
}


# --- Quantile derivative estimation (with bootstrapping) ---


# Derivative estimation using npregfast (interpolation with der1frfast)
quantder <- function(x, y, quant = 0.5, alpha = 0.05, 
                     nboot = 100, p = 1, kbin = 100, h = 0, hres = 0) {
  
  # Mean estimation
  mu <- der0(y = y, x = x, newdat = NULL, p = p, kbin = kbin, h = h) 
  res <- log((y - mu)^2)
  sd <- der0(y = res, x = x, newdat = NULL, p = p, kbin = kbin, h = hres)
  sd <- sqrt(exp(sd))
  resf <- (y - mu) / sd
  
  q <- mu + sd * quantile(resf, prob = quant)
  
  # First derivative
  q1 <- der1frfast(y = q, x = x)
  
  # Bootstrap samples
  yboot <- replicate(
    nboot,
    mu + sd * sample(resf, size = length(y), replace = TRUE)
  )
  
  # Bootstrap estimation
  muboot <- apply(yboot, 2, der0, x = x, newdat = NULL, 
                  p = p, kbin = kbin, h = h)
  resboot <- log((yboot - muboot)^2)
  sdboot <- apply(resboot, 2, der0, x = x, newdat = NULL, 
                  p = p, kbin = kbin, h = hres)
  sdboot <- sqrt(exp(sdboot))
  resbootf <- (yboot - muboot) / sdboot
  
  aux <- apply(resbootf, 2, quantile, prob = quant)
  qboot <- muboot + sweep(sdboot, MARGIN = 2, aux, "*") 
  
  # Recenter
  sesgo <- q - rowMeans(qboot)
  qboot <- qboot + sesgo
  
  # Bootstrap first derivative
  q1boot <- apply(qboot, 2, der1frfast, x = x)
  
  ic <- apply(qboot, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2))
  ic1 <- apply(q1boot, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2))
  
  res <- list(q = cbind(q, t(ic)), q1 = cbind(q1, t(ic1)))
}



# --- First and second derivative estimation ---


# First derivative at original data or at new data points
der0 <- function(y, x, newdat = NULL, ...) {
  if (is.null(newdat)) {
    predict(
      frfast(y ~ x, data = data.frame(x, y), nboot = 0, ...), 
      newdata = data.frame(x)
    )$Est[, 1]
  } else {
    predict(
      frfast(y ~ x, data = data.frame(x, y), nboot = 0, ...), 
      newdata = newdat
    )$Est[, 1]
  }
}

der1 <- function(y, x, newdat = NULL, ...) {
  if (is.null(newdat)) {
    predict(
      frfast(y ~ x, data = data.frame(x, y), nboot = 0, ...), 
      newdata = data.frame(x)
    )$First_deriv[, 1]
  } else {
    predict(
      frfast(y ~ x, data = data.frame(x, y), nboot = 0, ...), 
      newdata = newdat
    )$First_deriv[, 1]
  }
}

# First derivative at nodes
der0_nodes <- function(y, x) {
  frfast(y ~ x, data = data.frame(x, y), nboot = 0)$p[, 1, ]
}

# First derivative using linear model
der0_lm <- function(y, x) {
  as.numeric(predict(lm(y ~ x), type = "response"))
}

# First derivative using frfast with fixed parameters
der1frfast <- function(y, x, newdat = NULL, ...) {
  if (is.null(newdat)) {
    predict(
      frfast(y ~ x, data = data.frame(x, y), h = 0, p = 2, 
             kbin = 100, nboot = 0, ...),
      newdata = data.frame(x)
    )$First_deriv[, 1]
  } else {
    predict(
      frfast(y ~ x, data = data.frame(x, y), h = 0, p = 2, 
             kbin = 100, nboot = 0, ...),
      newdata = newdat
    )$First_deriv[, 1]
  }
}




# --- Analytical derivative from log-log quantile regression ---



d1 <- function(x, y, quant, kbin, grid) {
  xlog <- log(x)
  ylog <- log(y)
  
  model <- rq(ylog ~ xlog, data = data.frame(xlog, ylog), tau = quant)
  a <- exp(coef(model)[1])
  b <- coef(model)[2]
  
  if (is.null(grid)) {
    grid <- seq(min(x), max(x), length.out = kbin)
  }
  
  der <- a * b * grid^(b - 1)
  return(der)
}





  
  
 # --- Bootstrap confidence intervals for analytical derivative --- 

der1Qallo <- function(x, y, quant, nboot = 100, alpha = 0.05, kbin = 100) {
  
  grid <- seq(min(x), max(x), length.out = kbin)
  if (alpha == 0.10) {
    grid <- seq(3, 155, length.out = kbin)
  }
  
  mu1 <- d1(x, y, quant = quant, kbin, grid = grid)
  mu1boot <- matrix(NA, ncol = nboot, nrow = kbin)
  
  for (iboot in 1:nboot) {
    ii <- sample.int(length(y), size = length(y), replace = TRUE)
    mu1boot[, iboot] <- d1(x = x[ii], y[ii], quant = quant, kbin, grid = grid)
  }
  
  ic <- apply(mu1boot, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2))
  
  # Recenter
  sesgo <- mu1 - apply(ic, 2, mean)
  ic[1, ] <- ic[1, ] + sesgo
  ic[2, ] <- ic[2, ] + sesgo
  
  res <- list(q = cbind(mu1, t(ic)), grid = grid)
}

       




