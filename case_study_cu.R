library(jagsNEC)
library(tidyverse)

dat <- read.csv("ignore/C. armigera combined.csv") |> 
  dplyr::filter(Metal == "Cu") |> 
  dplyr::mutate(log.x = log(Conc + 1))
head(dat)

jagsfits <- fit.jagsMANEC(data = dat,
                          x.var = "log.x",
                          y.var = "Pcon")

check.chains(jagsfits,pdf.file="cu_chains.pdf")
jagsfits <- modify_jagsMANEC(jagsfits, drop.models = c("NECHormesis",
                                                       "ECxWeibull2", 
                                                       "ECx4param")) 

plot(jagsfits)
plot(jagsfits, all_models = TRUE)
jagsfits$mod.stats |> 
  dplyr::mutate(wi=round(wi,3),
                DIC=round(DIC,3)) |> 
  dplyr::select(DIC, wi) |> 
  write.csv("cs_cu_weights.csv")


#NEC function used to derive the no effect concentration where possible. 
library(drc)
Cu.nec <- drm(dat$Pcon, dat$Conc, data = dat, fct = NEC.3())
AIC(Cu.nec)
plot(Cu.nec, type='all')
summary(Cu.nec)
confint(Cu.nec, level = 0.95)


i=1
modnames <- names(jagsfits$mod.fits)
all_fit_dat <- lapply(jagsfits$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model") |> 
  rbind(data.frame(Model="averaged",jagsfits$pred.vals[c("x", "y", "up", "lw")]))

head(all_fit_dat)

all_endpoint_dat <- rbind(lapply(jagsfits$mod.fits, FUN = function(x){
    x$NEC })|> bind_rows(.id = "Model") |> data.frame(), 
    c(Model="averaged", jagsfits$NEC)) |> 
  mutate(estimate = if_else(grepl("NEC", Model), "NEC", 
                            if_else(grepl("ECx", Model), "NSEC",
                                    "N(S)EC")))


save(all_endpoint_dat, all_fit_dat, file = "cu_pred_vals.RData")    




