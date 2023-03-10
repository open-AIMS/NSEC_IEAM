
require(jagsNEC)
require(tidyverse)
library(ggplot2)
library(ggthemes) 
library(ggnewscale)
require(kableExtra)
source('functions.R')

load(file="data_raw/binomial_sim_results.RData")

x.range <- c(0.1, 10)
x.vec <- seq(x.range[1], x.range[2], length=50)

# set a seed
set.seed <- 1000

# create the binomial examples
sim.dat.binom <- create_sim_dat(scenario.dat)

sim.params <- list()
for(i in 1:length(sim.dat.binom)){
  name.fit <- names(sim.dat.binom)[i]
  if(length(grep("NEC", name.fit))==1){
    dat.i <- sim.dat.binom[[i]] %>%
      data.frame() %>%
      mutate(trials=as.integer(trials),
             top=as.numeric(top),
             NEC=as.numeric(NEC),
             beta=as.numeric(beta),
             suc=as.integer(suc))
    
    dat.params <- dat.i %>%
      select(top, beta, NEC) %>%
      unique() 
    
    y.pred <- NEC3param(x=x.vec, top=dat.params$top, beta=dat.params$beta, NEC=dat.params$NEC)     
  }
  if(length(grep("Sigmoidal", name.fit))==1){
    dat.i <- sim.dat.binom[[i]] %>%
    data.frame() %>%
    mutate(trials=as.integer(trials),
           top=as.numeric(top),
           d=as.numeric(d),
           beta=as.numeric(beta),
           suc=as.integer(suc))
  
    dat.params <- dat.i %>%
     select(top, beta, d) %>%
     unique() 
  
    y.pred <- ECxsigmoidal(x=x.vec, top=dat.params$top, beta=dat.params$beta, d=dat.params$d)    
    
  }
  sim.params <- c(sim.params, list(list(y.pred=y.pred, x.vec=x.vec, dat.i=dat.i, dat.params=dat.params)))
}



# make the new complete set
fitted.binom <- sapply(names(sim.dat.binom),function(x) NULL)

# fit the models
for(i in 1:length(sim.dat.binom)){
  name.fit <- names(sim.dat.binom)[i]
  dat.i <- sim.dat.binom[[i]] %>%
    data.frame() %>%
    mutate(trials=as.integer(trials),
           suc=as.integer(suc))
  

  fit.i <- try(fit.jagsMANEC(data = dat.i, x.var = "x", y.var = "suc", trials.var = "trials",
                             model.set = c("NEC3param", "ECxsigmoidal"), 
                             burnin = 10000))
  if(class(fit.i)!="try-error"){
    out.i <- fit.i
  }else{
    out.i <- "error"
  }
  
  fitted.binom[[name.fit]] <- out.i

}

names(fitted.binom) <- names(sim.dat.binom)
names(sim.params) <- names(sim.dat.binom)
save(sim.dat.binom, fitted.binom, sim.params, file="fitted_binom_examples.RData")


