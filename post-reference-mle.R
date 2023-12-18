require(coda)
require(maxLik)
require(VGAM)
require(MCMCpack)
require(numDeriv)
require(xtable)

# Preview -----------------------------------------------------------------

set.seed(2023)
eulerg <- 0.57721566490153
burnin <- 500
jump <- 5

month <- c("January", "February", "March", "April", "May", "June", "July",
           "August", "September", "October", "November", "December")

Table_mu <- 
  data.frame(
    Months = month,
    Reference = NA,
    MLE = NA
  )

Table_varp <- 
  data.frame(
    Months = month,
    Reference = NA,
    MLE = NA
  )


Table_mu_ci <- 
  data.frame(
    Months = month,
    Reference = NA,
    MLE = NA
  )

Table_varp_ci <- 
  data.frame(
    Months = month,
    Reference = NA,
    MLE = NA
  )

df <- rio::import("df_1947_2023.xls")

# Go! ---------------------------------------------------------------------

Generator <- function (index) {
  temp <- df[, index + 1]
  x <- as.numeric(temp[which(temp != "---")])
  n <- length(x)
  media <- mean(x)
  
  b <- 1
  m <- 1
  R <- 10000
  mu <- length(R + 1)
  varp <- length(R + 1)
  inimu <- media - (eulerg * sd(x) * sqrt(6) / pi)
  inivarp <- sd(x) * sqrt(6) / pi
  mu[1] <- inimu
  varp[1] <- inivarp
  
  xp <- c()
  c3 <- rep(0, times = R) 
  c4 <- rep(0, times = R)
  a3 <- 0
  a4 <- 0
  i <- 1
  c10 <- 0
  
  posterior3 <- function (v) {
    p <- (-n * media / v) - (n) * log(v) - n * log(sum(exp(-x / v)))
    return(p)
  }
  
  posterior4 <- function (v, mu) {
    p <- -(n + 1) * log(v) - (n * (media - mu) / v) - sum(exp(-(x - mu) / v))
    return(p)
  }
  
  loglike <- function (theta) {
    mu <- theta[1]
    v <- theta[2]
    aux <- -n * log(v) - (n * (media - mu) / v) - sum(exp(-(x - mu) / v))
    return(aux)
  }
  
  fisher <- function(mu, var) {
    I2 <- matrix(nrow = 2, ncol = 2)
    I2[1, 1] <- 1 / (var ^ 2)
    I2[1, 2] <- -(1 / (var ^ 2)) * (1 + digamma(1))
    I2[2, 1] <- I2[1, 2]
    I2[2, 2] <- (1 / (var ^ 2)) * (1 + trigamma(2) + digamma(2) ^ 2)
    return(I2)
  }
  
  try(while (i <= R) {
    if (i < 1)
      i <- 2
    prop3 <- rgamma(1, shape = b * varp[i], rate = b)
    ratio3 <-
      posterior3(prop3) - posterior3(varp[i]) + 
      dgamma(varp[i], shape = b * prop3, rate = b, log = T) - 
      dgamma(prop3, shape = b * varp[i], rate = b, log = T)
    has <- min(1, exp(ratio3))
    u3 <- runif(1)
    if (u3 < has &
        is.double(has) &
        has != "NaN" &
        has != "Inf" &
        has != "-Inf" &
        prop3 > 0.01) {
      varp[i + 1] <- prop3
      c3[i] <- 0
      a3 <- 0
    } else {
      varp[i + 1] <- varp[i]
      a3 <- a3 + 3
      c3[i] <- 1
    }
    
    prop4 <- rnorm(1, mean = mu[i], sd = m)
    ratio4 <-
      posterior4(varp[i + 1], prop4) - posterior4(varp[i + 1], mu[i]) + 
      dnorm(mu[i], mean = prop4, sd = m, log = T) - 
      dnorm(prop4, mean = mu[i], sd = m, log = T)
    alpha4 <- min(1, exp(ratio4))
    u4 <- runif(1)
    
    if (u4 < alpha4 &
        alpha4 != "NaN" &
        alpha4 != "Inf" &
        alpha4 != "-Inf") {
      mu[i + 1] <- prop4
      c4[i] <- 0
      a4 <- 0
    } else {
      mu[i + 1] <- mu[i]
      a4 <- a4 + 1
      c4[i] <- 1
    }
    if (a3 == 50 | a4 == 50) {
      i <- i - 50
      a3 = 0
      a4 = 0
    }
    i <- i + 1
    c10 <- c10 + 1
    if (c10 == 50000) {
      i = R + 1
      mu <- rep(0, times = R)
      varp <- rep(0, times = R)
    }
  })
  
  try(vvarp <- varp[seq(burnin, R, jump)])
  try(vmu <- mu[seq(burnin, R, jump)])
  
  fit <- 
    maxLik(loglike, start = c(1, 1), 
           constraints = list(ineqA = matrix(c(0, 1), ncol = 2), ineqB = 0))
  
  mle_est <- fit$estimate
  mle_se <- sqrt(diag(solve(n * fisher(mle_est[1], mle_est[2]))))
  
  old <- options(digits = 5)
  ace3 <- (1 - sum(c3) / length(c3))
  atc3 <- mean(acf(vvarp, plot = F)$acf)
  atc4 <- mean(acf(vmu, plot = F)$acf)
  ge3 <- abs(geweke.diag(vvarp)$z[1])
  ge4 <- abs(geweke.diag(vmu)$z[1])
  
  out <-
    matrix(c(
      paste0(round(mean(vmu), 4), " (", round(sd(vmu), 4), ")"),
      paste0(round(mle_est[1], 4), " (", round(mle_se[1], 4), ")"),
      paste0(round(mean(vvarp), 4), " (", round(sd(vvarp), 4), ")"),
      paste0(round(mle_est[2], 4), " (", round(mle_se[2], 4), ")"),
      paste0("(", round(quantile(vmu, 0.025), 4), "; ", round(quantile(vmu, 0.975), 4), ")"),
      paste0("(", round(mle_est[1] - 1.96 * mle_se[1], 4), "; ", round(mle_se[1] + 1.96 * mle_se[1], 4), ")"),
      paste0("(", round(quantile(vvarp, 0.025), 4), "; ", round(quantile(vvarp, 0.975), 4), ")"),
      paste0("(", round(mle_est[2] - 1.96 * mle_se[2], 4), "; ", round(mle_se[2] + 1.96 * mle_se[2], 4), ")")
    ),
    ncol = 2, byrow = T)
  
  return(out)
}

for (i in 1:12) {
  aux <- Generator(index = i)
  Table_mu[i, 2:3] <- aux[1,]
  Table_varp[i, 2:3] <- aux[2,]
  Table_mu_ci[i, 2:3] <- aux[3, ]
  Table_varp_ci[i, 2:3] <- aux[4, ]
}

print(xtable(Table_mu), include.rownames = F)
print(xtable(Table_varp), include.rownames = F)

print(xtable(Table_mu_ci), include.rownames = F)
print(xtable(Table_varp_ci), include.rownames = F)

