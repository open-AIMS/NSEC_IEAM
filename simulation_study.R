require(jagsNEC)
require(tidyverse)
require(doParallel)
source('functions.R')


## SIMULATION STUDY  --------
# Plot the scenarios ---------
x.range <- c(0.1, 10)
scenarios.list <- list("Sigmoidal 1"= c(top=0.9, beta=5e-03, NEC = NA, d=3.5), 
                       "NEC 1" = c(top=0.9, beta=0.75, NEC = 2.5, d=NA), 
                       "Sigmoidal 2" = c(top=0.9, beta=1e-08, NEC = NA, d=9), 
                       "NEC 2" = c(top=0.9, beta=0.75, NEC = 6, d=NA))

scenarios <- names(scenarios.list)
scenarios.dat <- do.call("rbind", scenarios.list) %>%
  data.frame() %>%
  rownames_to_column(var = "scenarios")

x.vec <- seq(x.range[1], x.range[2], length=50)

# Build scenario data -------
scenario.list <- list(
  n.treatments = c(8, 12),
  n = c(5, 10),
  top = 0.9,
  trials = c(10, 20),
  scenarios = scenarios
) 

scenario.dat <- expand.grid(scenario.list) %>%
  data.frame() 

scenario.dat$key <- unlist(unite(scenario.dat, col="key"))

scenario.dat <- scenario.dat %>%
  left_join(scenarios.dat) %>%
  mutate(model=ifelse(grepl("Sigmoidal", scenarios),"Sigmoidal", "NEC"))

# try several times in case there is an issue
for(r in 1:25){
  ## make a complete list of fitted sims ---------
  fitted.binom.allsims <- list()
  
  # run new simulations ----------------
  
  cl <- makePSOCKcluster(6)
  registerDoParallel(cl)
  
  n.sims <- 6*15
  
  out <- foreach(n=1:n.sims, .packages=c("jagsNEC", "tidyverse")) %dopar% {
    # create new simulated data
    sim.dat.binom <- create_sim_dat(scenario.dat)
    # make the new complete set
    fitted.binom <- sapply(names(sim.dat.binom),function(x) NULL)
    # fit the models
    for(i in 1:length(sim.dat.binom)){
      name.fit <- names(sim.dat.binom)[i]
      dat.i <- sim.dat.binom[[i]]
      fit.i <- try(fit.jagsMANEC(data = dat.i, x.var = "x", y.var = "suc", trials.var = "trials",
                                 model.set = c("NEC3param", "ECxsigmoidal"), 
                                 n.tries=10, burnin = 10000))
      if(class(fit.i)!="try-error"){
        out.i <- extract_output(fit=fit.i)
      }else{
        out.i <- "error"
      }
      fitted.binom[[name.fit]] <- out.i
    }
    return(fitted.binom)
    
  }
  stopCluster(cl)
  
  fitted.binom.allsims <- c(fitted.binom.allsims, out)
  save(fitted.binom.allsims, file=paste("binomial_allsims_", Sys.Date(), ".RData", sep=""))
  
}

# join the old and new output --------------------
# Note files can be downloaded from zip folder binomial_allsims_original.zip
load(file="data_raw/binomial_allsims_original.RData")
length(fitted.binom.allsims)

all.fitted.temp <- fitted.binom.allsims
load(file="data_raw/binomial_allsims_new_2020-08-14.RData")
length(fitted.binom.allsims)
fitted.binom.allsims <- c(fitted.binom.allsims, all.fitted.temp)
length(fitted.binom.allsims)

all.fitted.temp <- fitted.binom.allsims
load(file="data_raw/binomial_allsims_new_2020-08-13.RData")
length(fitted.binom.allsims)
fitted.binom.allsims <- c(fitted.binom.allsims, all.fitted.temp)
length(fitted.binom.allsims)

# process the output ---------------------
# weight dat      
weight.dat <- list()
for(l in 1:length(fitted.binom.allsims)){
  h=fitted.binom.allsims[[l]]
  tt1 <- lapply(h, FUN=function(x){
    if(x[1]!="error"){
      out = x$mod.stats %>%
        rownames_to_column(var="model") %>%
        select(model, wi)
    }else{
      out <- data.frame("model"=NA, "wi"=NA)
    }
    return(out)
  })
  tt1 <- bind_rows(tt1, .id="key") %>%
    pivot_wider(names_from = model, values_from = wi) %>%
    select(key, NEC3param, ECxsigmoidal) %>%
    mutate(sim=l) %>%
    data.frame()
  weight.dat <- c(weight.dat, list(tt1))
}

weight.dat <- do.call("rbind", weight.dat) %>%
  na.omit()

temp.list <- list()
for(n in 1:length(fitted.binom.allsims)){
  temp.list <- c(temp.list, list(extract_sim_dat(fitted.binom.allsims=fitted.binom.allsims, n=n)))
}

# collate output with the scenario data ---------------
head(scenario.dat)
# calculate theoretical ECx values for each simulation

x.vec <- seq(x.range[1], x.range[2], length=1000)

ECx.vals.out <- list()

for(l in 1:nrow(scenario.dat)){
  scen <- scenario.dat[l,]
  if(scen[,"model"]=="Sigmoidal"){
    y.vec <- ECxsigmoidal(x=x.vec, top=scen[,"top"], beta=scen[,"beta"], d=scen[,"d"])}
  
  if(scen[,"model"]=="NEC"){   
    y.vec <- NEC3param(x=x.vec, NEC=scen[,"NEC"], top=scen[,"top"], beta=scen[,"beta"])}
  
  ECx.vals <- list(EC10 = 10,EC5 = 5,EC1 = 1)
  
  range.y <- c(0, scen[,"top"])
  
  ECx.vals.out <- c(ECx.vals.out, list(do.call("cbind", ECx.out <- lapply(ECx.vals, FUN=function(ecx){
    ECx.y <- max(range.y) - diff(range.y) * (ecx / 100)
    ECx.x <- x.vec[which.min(abs(y.vec - ECx.y))]
    return(ECx.x)
  }))))
  
}

ECx.dat <- do.call("rbind", ECx.vals.out)

out <- cbind(scenario.dat, ECx.dat)  

head(out)

all.sim.out <- do.call("rbind", temp.list) %>%
  data.frame() %>%
  right_join(weight.dat) %>% 
  pivot_wider(names_from = rowname, values_from = value) %>%
  data.frame() %>%
  left_join(out) %>%
  mutate(total.reps=trials*n*n.treatments,
         n.treatments=factor(n.treatments))
colnames(all.sim.out)
str(all.sim.out)
dim(all.sim.out)

true.endpoints <- list(
  actual.NEC = all.sim.out %>%
    mutate(y=NEC) %>%
    select(scenarios, y) %>%
    unique() %>%
    mutate(endpoint="NEC"),
  actual.EC1 = all.sim.out %>%
    mutate(y=EC1) %>%
    select(scenarios, y) %>%
    unique() %>%
    mutate(endpoint="EC1"),
  actual.EC5 = all.sim.out %>%
    mutate(y=EC5) %>%
    select(scenarios, y) %>%
    unique() %>%
    mutate(endpoint="EC5"),
  actual.EC10 = all.sim.out %>%
    mutate(y=EC10) %>%
    select(scenarios, y) %>%
    unique() %>%
    mutate(endpoint="EC10")
)

save(scenarios.list, scenario.dat, all.sim.out, true.endpoints, 
     file="data_raw/binomial_sim_results.RData")