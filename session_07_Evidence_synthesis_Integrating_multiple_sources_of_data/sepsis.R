# Sepsis example, Ohlsson & Lacy 2013
library(rjags)
library(tidyverse)
library(posterior)

# control data is in the first column, treatment in the second
sepsis_dat <- list(
  Ns = 10,
  Na = 2,
  y = matrix(c(23, 8, 5, 14, 209, 5, 13, 13, 8, 39, 20, 2, 0, 8, 186, 4,
               10, 19, 3, 40),
             nrow = 10,
             ncol = 2),
  n = matrix(c(65, 43, 59, 32, 1212, 50, 34, 41, 40, 381, 61, 43, 56, 34,
               1204, 100, 68, 40, 40, 372),
             nrow = 10,
             ncol = 2)
)

# for use in plotting
StudyNames <- c("Bussel (1990)","Chirico (1987)","Clapp (1989)",
                "Conway (1990)","Fanaroff (1994)","Haque (1986)",
                "Ratrisawadi (1991)","Sandberg (2000)","Tanzer (1997)",
                "Weisman (1994)","Mean")

# independent effects model
sepsis_ie_model <- "
model
{
  # For each study, Ns = total number of studies
  for(i in 1:Ns)
  {
    # for each of the two arms
    for(k in 1:Na)
    {
      # Binomial likelihood
      y[i,k] ~ dbin(p[i,k], n[i,k])
    }
      
    # on logit scale, proportion is probability of success in terms of 
    # study baselines mu and study-specific treatment contrasts delta 
    # (log odds ratio, relative to study baseline, =0 for arm 1)
    logit(p[i,1]) <- mu[i]
    logit(p[i,2]) <- mu[i] + delta[i]
    
    # for the treatment arm, contrasts (log odds ratios) are independent effects
    # each with a vague Normal prior
    delta[i] ~ dnorm(0, 0.01)
  
    # study-specific baselines, vague priors
    mu[i] ~ dnorm(0, 0.01)
  }
}
"

# Inits
sepsis_ie_inits <- list(
  list(mu = rep(0,sepsis_dat$Ns),
       delta = rep(0,sepsis_dat$Ns),
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(7195)
  ),
  list(mu = c(1,-1,-2,0,0,-2,1,0,2,2),
       delta = rep(1,sepsis_dat$Ns),
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(168422)
  )
)

# numbers of chains, burn-in iterations and iterations to keep
nChains <- 2
nBurn <- 1000
nIter <- 5000

# Initialise model
sepsis_ie_jm <- jags.model(textConnection(sepsis_ie_model),
                           data = sepsis_dat,
                           inits = sepsis_ie_inits,
                           n.chains = nChains)

# burn-in
update(sepsis_ie_jm, n.iter = nBurn)

# Parameters to monitor
sepsis_ie_params <- c("p","mu","delta")

# samples to keep
sepsis_ie_out <- coda.samples(sepsis_ie_jm,
                              variable.names = sepsis_ie_params,
                              n.iter = nIter,
                              n.thin = 1)


# common effect model
sepsis_ce_model <- "
model
{
  # For each study, Ns = total number of studies
  for(i in 1:Ns)
  {
    # for each of the two arms
    for(k in 1:Na)
    {
      # Binomial likelihood
      y[i,k] ~ dbin(p[i,k], n[i,k])
    }
    
    # now only have a single treatment effect parameter
    # so delta not indexed by study i anymore
    logit(p[i,1]) <- mu[i]
    logit(p[i,2]) <- mu[i] + delta
  
    # study-specific baselines, vague priors
    mu[i] ~ dnorm(0, 0.01)
  }
    
  # vague prior for the common effect
  delta ~ dnorm(0, 0.01)
}
"

# Inits
sepsis_ce_inits <- list(
  list(mu = rep(0,sepsis_dat$Ns),
       delta = 0,
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(7195)
  ),
  list(mu = c(1,-1,-2,0,0,-2,1,0,2,2),
       delta = 1,
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(168422)
  )
)


# Initialise model
sepsis_ce_jm <- jags.model(textConnection(sepsis_ce_model),
                           data = sepsis_dat,
                           inits = sepsis_ce_inits,
                           n.chains = nChains)

# burn-in
update(sepsis_ce_jm, n.iter = nBurn)

# samples to keep
sepsis_ce_out <- coda.samples(sepsis_ce_jm,
                              variable.names = sepsis_ie_params,
                              n.iter = nIter,
                              n.thin = 1)





# random effects model
sepsis_re_model <- "
model
{
  # For each study, Ns = total number of studies
  for(i in 1:Ns)
  {
    # for each of the two arms
    for(k in 1:Na)
    {
      # Binomial likelihood
      y[i,k] ~ dbin(p[i,k], n[i,k])
    }
    
    # on logit scale, proportion is probability of success in terms of 
    # study baselines mu and study-specific treatment contrasts delta 
    # (log odds ratio, relative to study baseline, =0 for arm 1)
    logit(p[i,1]) <- mu[i]
    logit(p[i,2]) <- mu[i] + delta[i]
    
    # for the treatment arm, contrasts (log odds ratios) are random effects 
    # with a common mean d
    delta[i] ~ dnorm(d, prec.d)
  
    # study-specific baselines, vague priors
    mu[i] ~ dnorm(0, 0.01)
  }

  # Priors for basic parameters:
  
  # mean log odds ratio of treatment vs control
  d ~ dnorm(0, 0.01)
  
  # sd of study-specific log odds ratios
  prec.d <- 1 / (sd.d * sd.d)
  sd.d ~ dunif(0,10)
}
"

# Inits
sepsis_re_inits <- list(
  list(d = 0,
       sd.d = 5,
       mu = rep(0,sepsis_dat$Ns),
       delta = rep(0,sepsis_dat$Ns),
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(7195)
  ),
  list(d = 0.1,
       sd.d = 1,
       mu = c(1,-1,-2,0,0,-2,1,0,2,2),
       delta = rep(1,sepsis_dat$Ns),
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(168422)
  )
)

# Initialise model
sepsis_re_jm <- jags.model(textConnection(sepsis_re_model),
                           data = sepsis_dat,
                           inits = sepsis_re_inits,
                           n.chains = nChains)

# burn-in
update(sepsis_re_jm, n.iter = nBurn)

# Parameters to monitor
sepsis_re_params <- c("p","mu","delta","d","sd.d")

# samples to keep
sepsis_re_out <- coda.samples(sepsis_re_jm,
                              variable.names = sepsis_re_params,
                              n.iter = nIter,
                              n.thin = 1)


# Combine treatment effect (odds ratio) estimates from all three models 
# in a single tibble for plotting
sepsis_post <- as_tibble(as_draws_matrix(sepsis_re_out), rownames = "Iteration") %>%
  select("Iteration",contains("delta"), "d") %>%
  pivot_longer(
    cols = contains("d"),
    names_to = "Study",
    values_to = "LOR"
  ) %>%
  mutate(
    Model = "Random effects"
  ) %>%
  bind_rows(
    as_tibble(as_draws_matrix(sepsis_ie_out), rownames = "Iteration") %>%
      select("Iteration",contains("delta")) %>%
      pivot_longer(
        cols = contains("delta"),
        names_to = "Study",
        values_to = "LOR"
      ) %>%
      mutate(
        Model = "Independent effects"
      )
  ) %>%
  bind_rows(
    as_tibble(as_draws_matrix(sepsis_ce_out), rownames = "Iteration") %>%
      select("Iteration",contains("delta")) %>%
      rename(d = delta) %>%
      pivot_longer(
        cols = contains("d"),
        names_to = "Study",
        values_to = "LOR"
      ) %>%
      mutate(
        Model = "Common effect"
      )
  ) %>%
  mutate(
    Study = factor(Study,
                   levels = c(paste0("delta[",1:10,"]"), "d"),
                   labels = StudyNames),
    OR = as.double(exp(LOR)),
    Model = factor(Model,
                   levels = c("Common effect","Independent effects","Random effects"))
  )

sepsis_summary <- sepsis_post %>%
  group_by(Model, Study) %>%
  summarise(
    Median = median(OR),
    Lower = quantile(OR, probs = 0.025),
    Upper = quantile(OR, probs = 0.975),
    .groups = "keep"
  ) %>%
  ungroup()


# Create the repeated datasets for the cross-validation
ycv <- as_tibble(sepsis_dat$y) %>%
  slice(rep(row_number(), 10)) %>%
  data.matrix %>% 
  array(dim = c(10,10,2)) %>%
  aperm(c(2,1,3))
ncv <- as_tibble(sepsis_dat$n) %>%
  slice(rep(row_number(), 10)) %>%
  data.matrix %>% 
  array(dim = c(10,10,2)) %>%
  aperm(c(2,1,3))
sepsis_cv_dat <- list(
  Ns = 10,
  Na = 2,
  ycv = ycv,
  ncv = ncv
)

# Initial values for sepsis cross-validation
sepsis_cv_inits <- list(
  list(d = rep(0,sepsis_cv_dat$Ns),
       sd.d = rep(5,sepsis_cv_dat$Ns),
       mu = t(structure(.Data = rep(rep(0,sepsis_cv_dat$Ns),sepsis_cv_dat$Ns),
                        .Dim = c(sepsis_cv_dat$Ns,sepsis_cv_dat$Ns))),
       delta = t(structure(.Data = rep(rep(0,sepsis_cv_dat$Ns),sepsis_cv_dat$Ns),
                           .Dim = c(sepsis_cv_dat$Ns, sepsis_cv_dat$Ns))),
       yrep = sepsis_dat$y,
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(7195)
  ),
  list(d = rep(0.1,sepsis_cv_dat$Ns),
       sd.d = rep(1,sepsis_cv_dat$Ns),
       mu = t(structure(.Data = rep(c(1,-1,-2,0,0,-2,1,0,2,2),sepsis_cv_dat$Ns),
                        .Dim = c(sepsis_cv_dat$Ns,sepsis_cv_dat$Ns))),
       delta = t(structure(.Data = rep(rep(1,sepsis_cv_dat$Ns),sepsis_cv_dat$Ns),
                           .Dim = c(sepsis_cv_dat$Ns, sepsis_cv_dat$Ns))),
       yrep = sepsis_dat$y+1,
       .RNG.name = c("base::Mersenne-Twister"),
       .RNG.seed = c(168422)
  )
)


