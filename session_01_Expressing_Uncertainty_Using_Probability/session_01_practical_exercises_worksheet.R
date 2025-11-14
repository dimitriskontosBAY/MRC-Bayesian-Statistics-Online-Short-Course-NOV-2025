# 1 Practical session (1): probability judgements ----
## 1 ----
n_sds <- qnorm(0.975) - qnorm(0.025)
n_sds

n_sds <- (qnorm(0.975, sd=1.234) - qnorm(0.025, sd=1.234))/1.234
n_sds

## 2 ----
upper <- 0.2
lower <- 0.02
mu <- 0.1

# for a 95% CI
sigma <- (upper - lower) / n_sds
sigma

(a <- (mu*(1 - mu)/sigma^2 - 1)*mu)
(b <- (mu*(1 - mu)/sigma^2 - 1)*(1 - mu))

a / (a+b)
qbeta(c(0.025, 0.5, 0.975), a, b)

## 3 ----
upper <- 0.2
lower <- 0.02
mu <- 0.1

# for a 50% CI
n_sds50 <- qnorm(0.75) - qnorm(0.25)
n_sds50
sigma50 <- (upper - lower) / n_sds50
sigma50

# (a <- (mu*(1 - mu)/sigma^2 - 1)*mu)
(a50 <- (mu*(1 - mu)/sigma50^2 - 1)*mu)

# (b <- (mu*(1 - mu)/sigma^2 - 1)*(1 - mu))
(b50 <- (mu*(1 - mu)/sigma50^2 - 1)*(1 - mu))


a50/(a50+b50)
qbeta(c(0.25, 0.5, 0.75), a50, b50)

## 4 ----
x <- ppoints(100)
plot(x, dbeta(x, a, b), type="l")
lines(x, dbeta(x, a50, b50), col="blue")

beta_dist <- function(mu, lower, upper , ci_width = 0.95){
  tail_prob <- (1 - ci_width)/2
  n_sds <- qnorm(1 - tail_prob) - qnorm(tail_prob)
  sigma <- (upper - lower) / n_sds
  # defining quantities that are used more than once 
  # (e.g. n here. This quantity is actually the “effective sample size” — as explained in the lecture
  n <- mu*(1 - mu)/sigma^2 - 1
  a <- n*mu
  b <- n*(1 - mu)
  c(a=a, b=b)
}

beta_dist(0.1, 0.02, 0.2) 
beta_dist(0.1, 0.02, 0.2, ci_width=0.5)


# 2 Practical session (2): predictions and decisions ----
install.packages("BiocManager")
BiocManager::install("Biobase")
install.packages("TailRank")
library(TailRank)

## 1 ----
### (a) ----
10

### (b) ----
1 - pbb(20, 100, a, b)

### (c) ----
qbb(c(0.025, 0.975), 100, a, b)

## 2 ----
R <- 200000
theta <- rbeta(R, a, b)
y_pred <- rbinom(R, 100, theta)
mean(y_pred > 20)
mean(c(18, 19, 21) > 20)
mean(c(FALSE, FALSE, TRUE))
mean(c(0, 0, 1))
(mcse <- sd(y_pred > 20) / sqrt(R))
quantile(y_pred, c(0.025, 0.975))

## 3 -----

nsim <- 100000
p_old <- rbeta(nsim, a, b)
odds_old <- p_old / (1 - p_old)

log_or <- rnorm(nsim, 1, 1)
or <- exp(log_or)

odds_new <- odds_old * or
p_new <- odds_new / (1 + odds_new)

mean(p_old)
mean(p_new)
quantile(p_new, c(0.025, 0.975))

## 4 ----
### (a) ----
1 - pnorm(145, mean=120, sd = sqrt(10^2 + 8^2))

### (b) ----
qgamma(c(0.025, 0.975), 16, 2)

### (c) ----
R <- 100000
sigma_rep <- rgamma(R, 16, 2)
X_rep <- rnorm(R, mean=120, sd=sqrt(10^2 + sigma_rep^2))
mean(X_rep > 145)

c(
  sd_fixed = sqrt(10^2 + 8^2),
  sd_uncertain = sd(X_rep)
)
# 3 Practical session (3): introduction to JAGS ----
# install.packages("rjags")
library(rjags)

## 1 Defining and supplying a model ----

mod <- "
model {
  pop_size <- 100
  theta ~ dbeta(4.168, 37.515)
  cases_pred ~ dbinom(theta, pop_size)
}
"
mod_jag <- jags.model(textConnection(mod), data=NULL)

# or another way
mod_file <- tempfile()
mod_file

cat(mod, file=mod_file)
mod_jag <- jags.model(mod_file, data=NULL)

## 2 Sampling from a model ----
### (a) ----

### (b) ----
pars <- c("theta", "cases_pred")
nsim <- 10000
sam_c <- coda.samples(mod_jag, pars, n.iter=nsim)
# This is a matrix with rows for samples, and columns for different variables.
sam_m <- sam_c[[1]]

## 3 Summarizing samples from a model ----
cases_pred <- sam_m[,"cases_pred"]
quantile(cases_pred, c(0.025, 0.975))
quantile(sam_m[,"theta"], c(0.025, 0.975))
