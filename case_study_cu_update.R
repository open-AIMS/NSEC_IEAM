library(jagsNEC)
library(tidyverse)

dat <- readxl::read_excel("ignore/Cryo combined with id.xlsx") %>% 
  group_by(id) %>% 
  mutate(norm=dbl/max(dbl)) %>% 
  mutate(sqrt.x=sqrt(conc)) %>% 
  data.frame()

head(dat)


dat %>% 
  ggplot()+
  geom_point(aes(x=sqrt.x, y=norm))

#non-normalised gamma, normalised gamma, and a normalised beta and check the 
#difference in ecx. Perhaps just run the normalised Beta first, and then we can 
#do that comparison only with the model with the highest weight.

#########Models ##########
#Run and load models: 

load("cu_norm_models.RData")

non_gamma_fit<-fit.jagsMANEC(data = dat,
                             x.var="sqrt.x",
                             y.var="dbl",
                             y.type="gamma")

norm_gamma <- fit.jagsMANEC(data = dat,
                            x.var = "sqrt.x",
                            y.var = "norm",
                            y.type="gamma")

norm_beta <- fit.jagsMANEC(data = dat,
                           x.var = "sqrt.x",
                           y.var = "norm",
                           y.type="beta")

#save(non_gamma_fit, norm_gamma, norm_beta,EC_out_avg,file = "cu_norm_models.RData")  


###non-normalised gamma 

non_gamma_fit_stats<- non_gamma_fit$mod.stats |> 
  dplyr::mutate(model=names(non_gamma_fit$mod.fits),               
                wi=round(wi,3),
                DIC=round(DIC,3),
                trial="non_gamma_fit") |> 
  dplyr::select(trial,model,DIC, wi)


plot(non_gamma_fit, all_models = TRUE)

non_gamma_endpoint <- rbind(lapply(non_gamma_fit$mod.fits, FUN = function(x){
  x$NEC })|> bind_rows(.id = "Model") |> data.frame(), 
  c(Model="averaged", non_gamma_fit$NEC)) |> 
  mutate(estimate = if_else(grepl("NEC", Model), "NEC", 
                            if_else(grepl("ECx", Model), "NSEC",
                                    "N(S)EC")))|>
  left_join(non_gamma_fit_stats %>% rename(Model=model) %>% select(Model, wi))

###normalised gamma

norm_gamma_stats<- norm_gamma$mod.stats |> 
  dplyr::mutate(model=names(norm_gamma$mod.fits),               
                wi=round(wi,3),
                DIC=round(DIC,3),
                trial="norm_gamma") |> 
  dplyr::select(trial,model,DIC, wi)

norm_gamma_endpoint <- rbind(lapply(norm_gamma$mod.fits, FUN = function(x){
  x$NEC })|> bind_rows(.id = "Model") |> data.frame(), 
  c(Model="averaged", norm_gamma$NEC)) |> 
  mutate(estimate = if_else(grepl("NEC", Model), "NEC", 
                            if_else(grepl("ECx", Model), "NSEC",
                                    "N(S)EC")))|>
  left_join(norm_gamma_stats %>% rename(Model=model) %>% select(Model, wi))

plot(norm_gamma, all_models = TRUE)

###normalised gamma

norm_beta_stats<- norm_beta$mod.stats |> 
  dplyr::mutate(model=names(norm_beta$mod.fits),               
                wi=round(wi,3),
                DIC=round(DIC,3),
                trial="norm_beta") |> 
  dplyr::select(trial,model,DIC, wi)

norm_beta_endpoint <- rbind(lapply(norm_beta$mod.fits, FUN = function(x){
  x$NEC })|> bind_rows(.id = "Model") |> data.frame(), 
  c(Model="averaged", norm_beta$NEC)) |> 
  mutate(estimate = if_else(grepl("NEC", Model), "NEC", 
                            if_else(grepl("ECx", Model), "NSEC",
                                    "N(S)EC")))|>
  left_join(norm_beta_stats %>% rename(Model=model) %>% select(Model, wi))

norm_beta_weight_stats <- lapply(norm_beta$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model")

plot(norm_beta, all_models = TRUE)

#comparing model weights 

weights<- rbind(non_gamma_fit_stats, norm_gamma_stats, norm_beta_stats) %>% 
  select(-DIC) %>% 
  pivot_wider(names_from = model, values_from = wi)

#########ECx values####

EC_vals<- c(1,10, 50)
EC_names<- c("EC1", "EC10", "EC50")

mod_list<- list(non_gamma_fit, norm_gamma, norm_beta)
mod_names<- c("non_gamma_fit", "norm_gamma", "norm_beta")


EC_out_avg<- lapply(mod_list, FUN=function(x){
  lapply(EC_vals, FUN=function(i){
  extract_ECx(x,ECx.val =i,posterior = F)
})
})

names(EC_out_avg)=mod_names

out1<- EC_out_avg[[1]] %>% unlist() %>%  data.frame()%>%rownames_to_column() %>%  mutate(trial = mod_names[1]) %>% rename("estimate"=".")
out2<- EC_out_avg[[2]] %>% unlist() %>%  data.frame()%>%rownames_to_column() %>%  mutate(trial = mod_names[2]) %>% rename("estimate"=".")
out3<- EC_out_avg[[3]] %>% unlist() %>%  data.frame()%>%rownames_to_column() %>%  mutate(trial = mod_names[3]) %>% rename("estimate"=".")

estimates<- rbind(out1, out2, out3) %>% 
  mutate(type=if_else(str_detect(rowname, "lw"), "lower", 
                      if_else(str_detect(rowname, "up"), "upper", "value"))) %>% 
  
  mutate(ECx=if_else(str_detect(rowname, "EC_10"), "EC10", 
                      if_else(str_detect(rowname, "EC_50"), "EC50", "EC1"))) %>% 
  pivot_wider(names_from=type, values_from=estimate, id_cols=c(-rowname)) %>% 
  mutate_if(is.double, .x^2)

estimates %>% 
  ggplot(aes(x=ECx, y=value, fill=trial))+
  geom_col(position="dodge")+
  geom_errorbar(aes(y=value, ymin=lower, ymax=upper, group=trial), position="dodge")+
  theme_bw()+
  scale_fill_brewer(type="qual")


#full posterior estimates
EC_out<- lapply(mod_list, FUN=function(h){
lapply(EC_vals, FUN=function(i){
  lapply(h$mod.fits, FUN = function(x) {
    extract_ECx(x,ECx.val =i,posterior = TRUE)
  })
  })
})


names(EC_out)=mod_names
names(EC_out[[1]])=EC_names
names(EC_out[[2]])=EC_names
names(EC_out[[3]])=EC_names



test<-  EC_out %>%
   map(~.x)
    mutate(ECx = names(EC_names),
           trial = names(mod_names))


map2_dfr(names())
map_dfr(EC_out, ~.x)

ECx <- EC_out %>% 
  tibble(pH=names(.), Value=.) %>% unnest() %>% 
  mutate(Label="EC")

EC_out<- lapply(EC_vals, FUN=function(i){
  lapply(norm_beta$mod.fits, FUN = function(x) {
    extract_ECx(x,ECx.val =i,posterior = TRUE)
  })
})

EC10.Zn <-  lapply(non_gamma_fit, FUN = function(x) {
  extract_ECx(x,ECx.val =10,posterior = TRUE)
})

EC50.Zn <-  lapply(Zn_fits_update, FUN = function(x) {
  extract_ECx(x,ECx.val =50, posterior = TRUE)
})


Zn_EC <- EC10.Zn %>% 
  tibble(pH=names(.), Value=.) %>% unnest() %>% 
  mutate(Label="EC")

Zn_stats <- rbind(Zn_EC, Zn_NEC)








weight_stats <- lapply(jagsfits$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model") |> 

check.chains(jagsfits,pdf.file="cu_chains.pdf")
#jagsfits <- modify_jagsMANEC(jagsfits, drop.models = c("NECHormesis",
 #                                                      "ECxWeibull2", 
  #                                                     "ECx4param")) 

plot(jagsfits)
plot(jagsfits, all_models = TRUE)
jagsfits$mod.stats |> 
  dplyr::mutate(wi=round(wi,3),
                DIC=round(DIC,3)) |> 
  dplyr::select(DIC, wi) |> 
  write.csv("cs_cu_weights.csv")


#NEC function used to derive the no effect concentration where possible. 
library(drc)
Cu.nec <- drm(dat$norm~ dat$conc, data = dat, fct = NEC.3())
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




