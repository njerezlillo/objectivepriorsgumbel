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
    Jeffreys = NA
  )

Table_varp <- 
  data.frame(
    Months = month,
    Jeffreys = NA
  )

Table_mu_ci <- 
  data.frame(
    Months = month,
    Jeffreys = NA
  )

Table_varp_ci <- 
  data.frame(
    Months = month,
    Jeffreys = NA
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
  c1 <- rep(0, times = R)
  c2 <- rep(0, times = R)
  a1 <- 0
  a2 <- 0
  i <- 1
  c10 <- 0
  
  posterior1 <- function (v) {
    p <- (-n * media / v) - (n + 1) * log(v) - n * log(sum(exp(-x / v)))
    return(p)
  }
  
  posterior2 <- function (v, mu) {
    p <- -(n + 2) * log(v) - (n * (media - mu) / v) - sum(exp(-(x - mu) / v))
    return(p)
  }
  
  loglike <- function (theta) {
    mu <- theta[1]
    v <- theta[2]
    aux <- -n * log(v) - (n * (media - mu) / v) - sum(exp(-(x - mu) / v))
    return(aux)
  }
  
  try(while (i <= R) {
    if (i < 1)
      i <- 2
    prop1 <- rgamma(1, shape = b * varp[i], rate = b)
    ratio1 <-
      posterior1(prop1) - posterior1(varp[i]) + 
      dgamma(varp[i], shape = b * prop1, rate = b, log = T) - 
      dgamma(prop1, shape = b * varp[i], rate = b, log = T)
    has <- min(1, exp(ratio1))
    u1 <- runif(1)
    if (u1 < has & is.double(has) &
        has != "NaN" &
        has != "Inf" &
        has != "-Inf" &
        prop1 > 0.01) {
      varp[i + 1] <- prop1
      c1[i] <- 0
      a1 <- 0
    } else {
      varp[i + 1] <- varp[i]
      a1 <- a1 + 1
      c1[i] <- 1
    }
    
    prop2 <- rnorm(1, mean = mu[i], sd = m)
    ratio2 <-
      posterior2(varp[i + 1], prop2) - posterior2(varp[i + 1], mu[i]) + 
      dnorm(mu[i], mean = prop2, sd = m, log = T) - 
      dnorm(prop2, mean = mu[i], sd = m, log = T)
    alpha2 <- min(1, exp(ratio2))
    u2 <- runif(1)
    if (u2 < alpha2 &
        alpha2 != "NaN" &
        alpha2 != "Inf" &
        alpha2 != "-Inf") {
      mu[i + 1] <- prop2
      c2[i] <- 0
      a2 <- 0
    } else {
      mu[i + 1] <- mu[i]
      a2 <- a2 + 1
      c2[i] <- 1
    }
    if (a1 == 50 | a2 == 50) {
      i <- i - 50
      a1 = 0
      a2 = 0
    }
    i <- i + 1
    c10 <- c10 + 1
  })
  
  try(vvarp <- varp[seq(burnin, R, jump)])
  try(vmu <- mu[seq(burnin, R, jump)])
  
  old <- options(digits = 5)
  ace1 <- (1 - sum(c1) / length(c1))
  atc1 <- mean(acf(vvarp, plot = F)$acf)
  atc2 <- mean(acf(vmu, plot = F)$acf)
  ge1 <- abs(geweke.diag(vvarp)$z[1])
  ge2 <- abs(geweke.diag(vmu)$z[1])
  
  xp <- sapply(c(0.50, 0.75, 0.90, 0.95, 0.975), function(z) vmu - vvarp * log(-log(z)))

  x_grid <- c(2, 4, 10, 20, 40)
  plot(x_grid, apply(xp, 2, mean), type = "l", main = month[index],
       xlab = "Return Period (Years)", ylab = "Return Level (mm)")
  lines(x_grid, apply(xp, 2, function(z) quantile(z, 0.975)), lty = 2)
  lines(x_grid, apply(xp, 2, function(z) quantile(z, 0.025)), lty = 2)
  legend("bottomright", c("Jeffreys prior", "Credible Intervals"),
         lty = c(1, 2), bty = "n")
  
  out <- 
    matrix(c(paste0(round(mean(vmu), 4), " (", round(sd(vmu), 4), ")"),
             paste0(round(mean(vvarp), 4), " (", round(sd(vvarp), 4), ")"),
             paste0("(", round(quantile(vmu, 0.025), 4), "; ", round(quantile(vmu, 0.975), 4), ")"),
             paste0("(", round(quantile(vvarp, 0.025), 4), "; ", round(quantile(vvarp, 0.975), 4), ")")
             ), ncol = 2, byrow = T)
  
  return(out)
}

pdf("Fig3.pdf", width = 12, height = 10)
par(mfrow = c(3, 4))
for (i in 1:12) {
  aux <- Generator(index = i)
  Table_mu$Jeffreys[i] <- aux[1, 1]
  Table_varp$Jeffreys[i] <- aux[1, 2]
  Table_mu_ci$Jeffreys[i] <- aux[2, 1]
  Table_varp_ci$Jeffreys[i] <- aux[2, 2]
}
par(mfrow = c(1, 1))
dev.off()

print(xtable(Table_mu), include.rownames = F)
print(xtable(Table_varp), include.rownames = F)

print(xtable(Table_mu_ci), include.rownames = F)
print(xtable(Table_varp_ci), include.rownames = F)
