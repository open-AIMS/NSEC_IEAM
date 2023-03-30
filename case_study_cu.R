library(jagsNEC)
library(tidyverse)

dat <- read.csv("ignore/C. armigera combined.csv") |> 
  dplyr::filter(Metal == "Cu") |> 
  dplyr::mutate(log.x = log(Conc + 1))
head(dat)

jagsfits <- fit.jagsMANEC(data = dat,
                          x.var = "log.x",
                          y.var = "Pcon")
jagsfits <- modify_jagsMANEC(jagsfits, drop.models = "NECHormesis") 

plot(jagsfits)
plot(jagsfits, all_models = TRUE)
jagsfits$mod.stats |> 
  dplyr::mutate(wi=round(wi,3),
                DIC=round(DIC,3)) |> 
  dplyr::select(DIC, wi) |> 
  write.csv("cs_cu_weights.csv")

i=1
modnames <- names(jagsfits$mod.fits)
all_fit_dat <- lapply(jagsfits$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model") |> 
  rbind(data.frame(Model="averaged",jagsfits$pred.vals[c("x", "y", "up", "lw")]))

head(all_fit_dat)

all_endpoint_dat <- lapply(jagsfits$mod.fits, FUN = function(x){
  data.frame(
    EC10=extract_ECx(x, posterior = TRUE)[c(2, 1, 3)],
    NSEC=extract_NSEC(x, posterior = TRUE)[c(2, 1, 3)],
    NEC=x$NEC[1:3]    
  )   
})|> bind_rows(.id = "Model") |> 
  data.frame(Model="averaged",
    EC10=extract_ECx(x)[c(2, 1, 3)],
    NSEC=extract_NSEC(x)[c(2, 1, 3)],
    NEC=x$NEC[1:3]   
  )  
    
