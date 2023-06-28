
library(multcomp)
library(tidyverse)
library(bayesnec)
library(glmmTMB)
library(jagsNEC)

require(cmdstanr)
set_cmdstan_path("C:/cmdstan")
options(brms.backend = "cmdstanr")

load("caste_study_results.RData")

dat <- Urchin.out2$mod.fits$NEC3param$mod.dat |> 
  bind_cols() |> 
  mutate(suc=as.integer(round(N*y)),
         Tr=as.factor(x)) |> 
  data.frame()

fit <- bnec(suc | trials(N) ~ crf(x, "ecxwb2"), data=dat, family=beta_binomial2)
fit2 <- bnec(y ~ crf(x, "ecxwb2"), data=dat)
fit3 <- bnec(suc | trials(N) ~ crf(x, "ecxwb2"), data=dat)

fit.glmer<-glmmTMB(y~Tr, family=Beta(link = "logit", link_phi = "log"), data=dat) 

summary(fit.glmer)
drop1(fit.glmer)
summary(glht(fit.glmer, linfct=mcp(Tr="Dunnett")))         


out1 <- fit.jagsNEC(
  data = dat,
  x.var = "x",
  y.var = "suc",
  trials.var = "N",
  model = "ECxWeibull2")

plot(out1)
out1$over.disp

out2 <- fit.jagsNEC(
  data = dat,
  x.var = "x",
  y.var = "y",
  model = "ECxWeibull2")
out2$over.disp

extract_ECx(out1)
extract_ECx(out2)

ec10 <- extract_ECx(Urchin.out2)


exp(Urchin.out2$NEC)
