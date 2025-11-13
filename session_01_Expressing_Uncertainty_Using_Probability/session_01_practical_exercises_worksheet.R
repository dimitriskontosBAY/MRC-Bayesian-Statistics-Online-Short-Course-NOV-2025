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
