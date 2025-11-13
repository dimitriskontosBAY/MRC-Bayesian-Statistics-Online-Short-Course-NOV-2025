# from the TailRank package (Kevin R. Coombes)
# https://CRAN.R-project.org/package=TailRank


pbb <- function (q, N, u, v) 
{
  sapply(q, function(xx) sum(dbb(0:xx, N, u, v)))
}

dbb <-function (x, N, u, v, log = FALSE) 
{
  logval <- lbeta(x + u, N - x + v) - lbeta(u, v) + lchoose(N, 
                                                            x)
  if (log) {
    ret <- logval
  }
  else {
    ret <- exp(logval)
  }
  ret
}

qbb <- function (p, N, u, v) 
{
  pp <- cumsum(dbb(0:N, N, u, v))
  sapply(p, function(x) sum(pp < x))
}

rbb <- function (n, N, u, v) 
{
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}
