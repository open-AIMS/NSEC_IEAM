require(tidyverse)
require(jagsNEC)
require(ggplot2)


# urchin
load("ignore/Urchin.RData")
Urchin.sigmoidal <- Urchin.out2$mod.fits$NECsigmoidal
Urchin.out2.NEC <- modify_jagsMANEC(Urchin.out2, model.set = "NEC")
Urchin.out2 <- modify_jagsMANEC(Urchin.out2, 
                                drop.models = c("NECsigmoidal", "ECxsigmoidal"))
Urchin.out2.NEC <- modify_jagsMANEC(Urchin.out2.NEC, 
                                    drop.models = c("NECsigmoidal"))
EC10.Urchin <- extract_ECx(Urchin.out2, posterior = TRUE)
NEC.Urchin <- Urchin.out2.NEC$sims.list$NEC
NSEC.Urchin <- extract_NSEC(Urchin.out2, posterior = TRUE)

# coral
load("ignore/coral_out_becky.RData")
Coral.sigmoidal <- Coral.out4$mod.fits$NECsigmoidal
Coral.out4.NEC <- modify_jagsMANEC(Coral.out4, model.set = "NEC") 
Coral.out4.all <- modify_jagsMANEC(Coral.out4, 
                                   drop.models = c("NECHormesis", "ECxsigmoidal", "ECx4param", "NECsigmoidal" ))
Coral.out4.NEC <- modify_jagsMANEC(Coral.out4.NEC, 
                                   drop.models = c("NECHormesis", "NECsigmoidal"))
EC10.Coral <- extract_ECx(Coral.out4.all, posterior = TRUE) 
NEC.Coral <- Coral.out4.NEC$sims.list$NEC
NSEC.Coral <- extract_NSEC(Coral.out4.all, posterior = TRUE)

# ND
load("ignore/NDARP.RData")
nd.sigmoidal <- out.nd.all$mod.fits$NECsigmoidal
out.nd.nec <- modify_jagsMANEC(out.nd.all, model.set="NEC")
out.nd.all <- modify_jagsMANEC(out.nd.all, 
                               drop.models=c("NECsigmoidal", "NEC4param", "ECxWeibull2", "ECx4param"))
out.nd.nec <- modify_jagsMANEC(out.nd.nec, 
                               drop.models=c("NECsigmoidal", "NEC4param"))
EC10.nd <- extract_ECx(out.nd.all, ECx.val = 10, posterior = TRUE)
NEC.nd <- out.nd.nec$sims.list$NEC
NSEC.nd <- extract_NSEC(out.nd.all, posterior = TRUE)

# AA
load("ignore/AAARP.Rdata")
aa.sigmoidal <- out.aa.all$mod.fits$NECsigmoidal
out.nec.AA <- modify_jagsMANEC(out.aa.all, model.set="NEC")
out.aa.all <- modify_jagsMANEC(out.aa.all, 
                               drop.models=c("NECsigmoidal", "NEC4param", "ECxWeibull1", "ECx4param"))
out.nec.AA  <- modify_jagsMANEC(out.nec.AA, 
                                drop.models=c("NECsigmoidal", "NEC4param"))

EC10.aa <- extract_ECx(out.aa.all, ECx.val = 10, posterior = TRUE)
NEC.aa <- out.nec.AA$sims.list$NEC
NSEC.aa <- extract_NSEC(out.aa.all, posterior = TRUE)


endpoints <- list(Urchin=list(NEC=NEC.Urchin, 
                              EC10=EC10.Urchin, 
                              NSEC=NSEC.Urchin),
                  Coral=list(NEC=NEC.Coral, 
                             EC10=EC10.Coral, 
                             NSEC=NSEC.Coral),
                  aa=list(NEC=NEC.aa, 
                          EC10=EC10.aa, 
                          NSEC=NSEC.aa),
                  nd=list(NEC=NEC.nd, 
                          EC10=EC10.nd, 
                          NSEC=NSEC.nd))
endpoints <- lapply(endpoints, FUN = function(x){lapply(x, FUN = function(y){
  y[sample(1:length(y), 100, replace = TRUE)]
})})

mod.weights.dat <-list(
  "S. variolaris"=rownames_to_column(Urchin.out2$mod.stats, var="Model"),
  "A. millepora"=rownames_to_column(Coral.out4.all$mod.stats, var="Model"),
  "A. amphitrite"=rownames_to_column(out.aa.all$mod.stats, var="Model"),
  "N. dorsatus"=rownames_to_column(out.nd.all$mod.stats, var="Model")) %>% 
  bind_rows(.id = "Species") %>% 
  dplyr::mutate(wi=round(wi, 3)) %>% 
  dplyr::arrange(Species, desc(wi)) %>%
  dplyr::select(Species, Model, wi) %>% 
  dplyr::filter(Model != "NECHormesis") %>% 
  pivot_wider(values_from = wi, names_from = Species)  %>% 
  arrange(Model)

xform1 <- function(x){exp(x)-10}
xform2 <- function(x){exp(x)}

mod.weights.long <-list(
  "S. variolaris"=rownames_to_column(Urchin.out2$mod.stats, var="Model"),
  "A. millepora"=rownames_to_column(Coral.out4.all$mod.stats, var="Model"),
  "A. amphitrite"=rownames_to_column(out.aa.all$mod.stats, var="Model"),
  "N. dorsatus"=rownames_to_column(out.nd.all$mod.stats, var="Model")) %>% 
  bind_rows(.id = "Species") %>% 
  dplyr::mutate(wi=round(wi, 3)) %>% 
  dplyr::arrange(Species, desc(wi)) %>%
  dplyr::select(Species, Model, wi) %>% 
  dplyr::filter(wi>0) %>% 
  arrange(Species)



save(endpoints, mod.weights.dat, mod.weights.long, xform1, xform2, 
     Urchin.out2, Urchin.out2.NEC,
     Coral.out4.all, Coral.out4.NEC,
     out.aa.all, out.nec.AA,
     out.nd.all, out.nd.nec, 
     Urchin.sigmoidal, Coral.sigmoidal, aa.sigmoidal, nd.sigmoidal,
     file = "caste_study_results.RData")

