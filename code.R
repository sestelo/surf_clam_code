# ======================================================================
# R Code for the manuscript 
# "Nonparametric quantile regression captures regional variability
# and scaling deviations in Atlantic surfclam length–weight relationships"
# 
# Author: Marta Sestelo - sestelo@uvigo.gal 
# ======================================================================

# --- Load required packages ---
library(npregfast)
library(quantreg)
library(Qtools)
library(patchwork)

# --- Load auxiliary functions ---
source("aux_functions.R")  # Ensure "functions.R" exists in your working directory

# --- Load the data ---
# Replace this with your actual data loading command
# Example: data <- read.csv("your_data.csv")



# ====================================================
# 1. Estimate Mean Regression Functions
# ====================================================

# --- Allometric model estimation ---
allo <- npregfast::frfast(
  formula = wt ~ lt,
  data = data,
  kbin = 100,
  model = "allo",
  seed = 300716
)

# --- Summary of the allometric model ---
summary(allo)  # Prints estimated coefficients and diagnostics

# --- Plot regression and first derivative for allometric model ---
ders_allo <- lapply(0:1, function(der_order) {
  npregfast::autoplot(
    object = allo,
    der = der_order,
    xlab = "LT (mm)",
    ylab = "WT (g)",
    lwd = 0.5
  ) +
    xlim(35, 190) +
    ylim(0, 355)
})


# --- Nonparametric model estimation ---
np <- npregfast::frfast(
  formula = wt ~ lt,
  data = data,
  kbin = 100,
  model = "np",
  p = 2,
  seed = 300716
)

# --- Plot regression and first derivative for nonparametric model ---
ders_np_all <- lapply(0:1, function(der_order) {
  npregfast::autoplot(
    object = np,
    der = der_order,
    xlab = "LT (mm)",
    ylab = "WT (g)",
    lwd = 0.5
  ) +
    xlim(35, 190) +
    ylim(0, 355)
})



# --- Predictions for specific length values (both models) ---

aux <- seq(50, 170, 20)
xvalues <- data.frame(lt = c(aux))
mu_allo <- predict(allo, newdata = xvalues, seed = 300716)
mu_np <- predict(np, newdata = xvalues, seed = 300716)





# --- Goodness-of-Fit Test for Allometric Model ---

npregfast::allotest(
  formula = wt ~ lt,
  data = data,
  nboot = 1000,
  seed = 300716,
  test = "lrt"
)





# ====================================================
# 2. Estimate Quantile Regression Functions
# ====================================================

# --- Prepare data for quantile regression ---

logwt <- log(data$wt)
loglt <- log(data$lt)
logdata <- data.frame(loglt, logwt)


# --- Allometric quantile regression estimation ---

q50 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.5),
  newdata = data.frame(loglt = loglt),
  interval = "confidence"
)

q10 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.1),
  newdata = data.frame(loglt = loglt),
  interval = "confidence"
)

q90 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.9),
  newdata = data.frame(loglt = loglt),
  interval = "confidence"
)

colors <- wesanderson::wes_palette("GrandBudapest1", 3)
ii <- order(exp(loglt))

p50 <- plotgg(
  x = exp(loglt)[ii],
  y = exp(logwt)[ii],
  muhat = exp(q50[ii, 1]),
  lci = exp(q50[ii, 2]),
  uci = exp(q50[ii, 3]),
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[1]
)

p10 <- plotgg(
  x = exp(loglt)[ii],
  y = exp(logwt)[ii],
  muhat = exp(q10[ii, 1]),
  lci = exp(q10[ii, 2]),
  uci = exp(q10[ii, 3]),
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[2]
)

p90 <- plotgg(
  x = exp(loglt)[ii],
  y = exp(logwt)[ii],
  muhat = exp(q90[ii, 1]),
  lci = exp(q90[ii, 2]),
  uci = exp(q90[ii, 3]),
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[3]
)

der_allo_q_all <- p50 +
  p10$layers[[2]] + p10$layers[[3]] +
  p90$layers[[2]] + p90$layers[[3]] +
  scale_color_manual(values = colors, guide = "none") +
  xlim(35, 190) +
  ylim(NA, 355)






# --- First derivative estimation for allometric quantile regression ---

d50 <- der1Qallo(
  x = data$lt, y = data$wt, quant = 0.5,
  nboot = 1000, kbin = 500
)

d10 <- der1Qallo(
  x = data$lt, y = data$wt, quant = 0.1,
  nboot = 1000, kbin = 500
)

d90 <- der1Qallo(
  x = data$lt, y = data$wt, quant = 0.9,
  nboot = 1000, kbin = 500
)

pd50 <- plotgg(
  x = d50$grid, y = NA, muhat = d50$q[, 1],
  lci = d50$q[, 2], uci = d50$q[, 3],
  xlab = "LT (mm)", ylab = "First derivative",
  col = colors[1], points = FALSE
)

pd10 <- plotgg(
  x = d10$grid, y = NA, muhat = d10$q[, 1],
  lci = d10$q[, 2], uci = d10$q[, 3],
  xlab = "LT (mm)", ylab = "First derivative",
  col = colors[2], points = FALSE
)

pd90 <- plotgg(
  x = d90$grid, y = NA, muhat = d90$q[, 1],
  lci = d90$q[, 2], uci = d90$q[, 3],
  xlab = "LT (mm)", ylab = "First derivative",
  col = colors[3], points = FALSE
)

der1_allo_q_all <- pd50 +
  pd10$layers[[2]] + pd10$layers[[1]] +
  pd90$layers[[2]] + pd90$layers[[1]] +
  scale_color_manual(
    name = "", values = colors,
    labels = c("Q90", "Q50", "Q10")
  ) +
  xlim(35, 190) +
  coord_cartesian(ylim = c(NA, 6.5))










# --- Nonparametric quantile regression estimation (function and derivative) ---


q50 <- quantder(
  x = data$lt, y = data$wt, quant = 0.5,
  nboot = 100, p = 2, kbin = 100,
  hres = -1, h = -1
)

q10 <- quantder(
  x = data$lt, y = data$wt, quant = 0.1,
  nboot = 100, p = 2, kbin = 100,
  hres = -1, h = -1
)

q90 <- quantder(
  x = data$lt, y = data$wt, quant = 0.9,
  nboot = 100, p = 2, kbin = 100,
  hres = -1, h = -1
)

colors <- wesanderson::wes_palette("GrandBudapest1", 3)
ii <- order(data$lt)

p50 <- plotgg(
  x = data$lt[ii],
  y = data$wt[ii],
  muhat = q50$q[ii, 1],
  lci = q50$q[ii, 2],
  uci = q50$q[ii, 3],
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[1]
)

p10 <- plotgg(
  x = data$lt[ii],
  y = data$wt[ii],
  muhat = q10$q[ii, 1],
  lci = q10$q[ii, 2],
  uci = q10$q[ii, 3],
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[2]
)

p90 <- plotgg(
  x = data$lt[ii],
  y = data$wt[ii],
  muhat = q90$q[ii, 1],
  lci = q90$q[ii, 2],
  uci = q90$q[ii, 3],
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[3]
)

der_np_q_all <- p50 +
  p10$layers[[2]] + p10$layers[[3]] +
  p90$layers[[2]] + p90$layers[[3]] +
  scale_color_manual(values = colors, guide = "none") +
  xlim(35, 190) +
  ylim(NA, 355)

# --- first derivative (nonparametric QR) ---

ylim <- c(NA, 5.5)

pd50 <- plotgg(
  x = data$lt[ii],
  y = data$wt[ii],
  muhat = q50$q1[ii, 1],
  lci = q50$q1[ii, 2],
  uci = q50$q1[ii, 3],
  points = FALSE,
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[1],
  ylim = ylim
)

pd10 <- plotgg(
  x = data$lt[ii],
  y = data$wt[ii],
  muhat = q10$q1[ii, 1],
  lci = q10$q1[ii, 2],
  uci = q10$q1[ii, 3],
  points = FALSE,
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[2],
  ylim = ylim
)


pd90 <- plotgg(
  x = data$lt[ii],
  y = data$wt[ii],
  muhat = q90$q1[ii, 1],
  lci = q90$q1[ii, 2],
  uci = q90$q1[ii, 3],
  points = FALSE,
  xlab = "LT (mm)",
  ylab = "WT (g)",
  col = colors[3],
  ylim = ylim
)

der1_np_q_all <- pd50 +
  pd10$layers[[2]] + pd10$layers[[1]] +
  pd90$layers[[2]] + pd90$layers[[1]] +
  scale_color_manual(
    name = "",
    values = colors,
    labels = c("Q90", "Q50", "Q10")
  ) +
  xlim(35, 190) +
  coord_cartesian(ylim = c(NA, 6.5))





# --- Predictions for specific length values (both models) ---

ss <- seq(50, 170, 20)
xvalues <- data.frame(loglt = c(ss))

nq10 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.1),
  newdata = log(xvalues),
  interval = "confidence"
)

nq25 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.25),
  newdata = log(xvalues),
  interval = "confidence"
)

nq50 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.5),
  newdata = log(xvalues),
  interval = "confidence"
)

nq75 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.75),
  newdata = log(xvalues),
  interval = "confidence"
)

nq90 <- predict(
  rq(logwt ~ loglt, data = logdata, tau = 0.9),
  newdata = log(xvalues),
  interval = "confidence"
)

q25 <- quantder(
  x = data$lt, y = data$wt, quant = 0.25,
  nboot = 500, p = 3, kbin = 100,
  hres = -1, h = -1
)

q75 <- quantder(
  x = data$lt, y = data$wt, quant = 0.75,
  nboot = 500, p = 3, kbin = 100,
  hres = -1, h = -1
)

npq10 <- predict(smooth.spline(x = data$lt, y = q10$q[, 1]), x = xvalues)$y
npq25 <- predict(smooth.spline(x = data$lt, y = q25$q[, 1]), x = xvalues)$y
npq50 <- predict(smooth.spline(x = data$lt, y = q50$q[, 1]), x = xvalues)$y
npq75 <- predict(smooth.spline(x = data$lt, y = q75$q[, 1]), x = xvalues)$y
npq90 <- predict(smooth.spline(x = data$lt, y = q90$q[, 1]), x = xvalues)$y

res <- data.frame(
  xvalues,
  exp(nq10[, 1]), exp(nq25[, 1]),
  exp(nq50[, 1]), exp(nq75[, 1]), exp(nq90[, 1]),
  npq10, npq25, npq50, npq75, npq90
)



# --- first derivative--- 

np1q10 <- predict(smooth.spline(x = data$lt, y = q10$q1[, 1]), x = xvalues)$y
np1q25 <- predict(smooth.spline(x = data$lt, y = q25$q1[, 1]), x = xvalues)$y
np1q50 <- predict(smooth.spline(x = data$lt, y = q50$q1[, 1]), x = xvalues)$y
np1q75 <- predict(smooth.spline(x = data$lt, y = q75$q1[, 1]), x = xvalues)$y
np1q90 <- predict(smooth.spline(x = data$lt, y = q90$q1[, 1]), x = xvalues)$y

d25 <- der1Qallo(x = data$lt, y = data$wt, quant = 0.25)
d75 <- der1Qallo(x = data$lt, y = data$wt, quant = 0.75)

allo1q10 <- predict(smooth.spline(x = d10$grid, y = d10$q[, 1]), x = xvalues)$y
allo1q25 <- predict(smooth.spline(x = d25$grid, y = d25$q[, 1]), x = xvalues)$y
allo1q50 <- predict(smooth.spline(x = d50$grid, y = d50$q[, 1]), x = xvalues)$y
allo1q75 <- predict(smooth.spline(x = d75$grid, y = d75$q[, 1]), x = xvalues)$y
allo1q90 <- predict(smooth.spline(x = d90$grid, y = d90$q[, 1]), x = xvalues)$y

res <- data.frame(
  xvalues,
  allo1q10, allo1q25,
  allo1q50, allo1q75, allo1q90,
  np1q10, np1q25, np1q50, np1q75, np1q90
)




# --- Goodness-of-Fit Test for Allometric Model (quantile) ---

model_quant <- rq(
  logwt ~ loglt,
  data = logdata,
  tau = c(0.1, 0.5, 0.9)
)

Qtools::GOFTest(
  model_quant,
  B = 1000,
  seed = 300716
)

