library(jagsNEC)
library(tidyverse)

dat <- readxl::read_excel("case_study1_data.xlsx") %>% 
  group_by(id) %>% 
  mutate(norm=dbl/max(dbl)) %>% 
  mutate(sqrt.x=sqrt(conc)) %>% 
  data.frame()

head(dat)


#check data: 

dat %>% 
  ggplot()+
  geom_point(aes(x=sqrt.x, y=norm))


######### Fit Models ##########
#comparing normalised vs non-normalised response variables
#(i.e. doublings per day or doublings as a percent of control doublings for each experiment)

#Run and load models: 
#cu_pred_vals.RData

non_gamma_fit<-fit.jagsMANEC(data = dat,
                             x.var="sqrt.x",
                             y.var="dbl",
                             y.type="gamma", 
                             )

norm_gamma <- fit.jagsMANEC(data = dat,
                            x.var = "sqrt.x",
                            y.var = "norm",
                            y.type="gamma")

norm_beta <- fit.jagsMANEC(data = dat,
                           x.var = "sqrt.x",
                           y.var = "norm",
                           y.type="beta")


save(non_gamma_fit, norm_gamma, norm_beta, file = "cu_mod_fits.RData")
load("Cu_mod_fits.RData")

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
  left_join(non_gamma_fit_stats %>% rename(Model=model) %>% dplyr::select(Model, wi))

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
  left_join(norm_gamma_stats %>% rename(Model=model) %>% dplyr::select(Model, wi))

plot(norm_gamma, all_models = TRUE)

###normalised beta

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
  left_join(norm_beta_stats %>% rename(Model=model) %>% dplyr::select(Model, wi))

norm_beta_weight_stats <- lapply(norm_beta$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model")

plot(norm_beta, all_models = TRUE)

######Model weights ######

weights<- rbind(non_gamma_fit_stats, norm_gamma_stats, norm_beta_stats) %>% 
  dplyr::select(-DIC) %>% 
  pivot_wider(names_from = model, values_from = wi)

#########Extracting ECx values####

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
  mutate_if(is.double, ~.^2)#convert back from sqrt transformation

#comparison of estimates and 
estimates %>% 
  mutate(trial= fct_recode(trial, 
    "Not normalised gamma" = "non_gamma_fit", 
    "Normalised beta"="norm_beta", 
    "Normalised gamma" = "norm_gamma")) %>% 
  ggplot(aes(x=ECx, y=value, fill=trial))+
  geom_col(position="dodge")+
  geom_errorbar(aes(y=value, ymin=lower, ymax=upper, group=trial), position="dodge")+
  theme_bw()+
  scale_fill_brewer(type="qual")+
  labs(y="Estimate (µg/L)")+
  theme(legend.position = "bottom")

#note that the original paper reported EC10 of 21.6 and EC50 of 63 

######progressing with the normalised beta
#removing the hormesis and NECsigmoidal models 
norm_beta_2 <- modify_jagsMANEC(norm_beta, 
                                drop.models = c("NECsigmoidal", "NECHormesis"))
dev.off()
plot(norm_beta_2)
plot(norm_beta_2, all_models = TRUE)


weight_stats <- lapply(norm_beta_2$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model") |> 
check.chains(norm_beta_2,pdf.file="cu_chains.pdf")

norm_beta_2$mod.stats |> 
  dplyr::mutate(wi=round(wi,3),
                DIC=round(DIC,3)) |> 
  dplyr::select(DIC, wi)|> 
 rownames_to_column("Model") |> 
  data.frame()|> 
  write.csv("Table3.csv")

i=1
modnames <- names(norm_beta_2$mod.fits)
all_fit_dat <- lapply(norm_beta_2$mod.fits, FUN = function(x){
  data.frame(x$pred.vals[c("x", "y", "up", "lw")])
}) |> bind_rows(.id = "Model") |> 
  rbind(data.frame(Model="averaged",norm_beta_2$pred.vals[c("x", "y", "up", "lw")]))

head(all_fit_dat)

all_endpoint_dat <- rbind(lapply(norm_beta_2$mod.fits, FUN = function(x){
    x$NEC })|> bind_rows(.id = "Model") |> data.frame(), 
    c(Model="averaged", norm_beta_2$NEC)) |> 
  mutate(estimate = if_else(grepl("NEC", Model), "NEC", 
                            if_else(grepl("ECx", Model), "NSEC",
                                    "N(S)EC")))

#extracting ECx values 

ECxs <-c(1,10)

#model average ECx 
Cu_avg_ECxs <- lapply(ECxs, FUN=function(x){
      extract_ECx(norm_beta_2,ECx.val =x,posterior = F)
  })
#individual model ECxs
Cu_avg_estimates<- rbind(Cu_avg_ECxs) %>% unlist() %>% data.frame() %>% rownames_to_column("lab") %>% 
  rename("value"=".") %>% 
  mutate(ECx=if_else(str_detect(lab, "EC_10"), "EC10", "EC1")) %>% 
  mutate(val=if_else(str_detect(lab, "lw"), "lw", 
                     if_else(str_detect(lab, "up"), "up",   "value"))) %>% 
  mutate(Model="Averaged")


Cu_ind_ECxs<- 
  lapply(norm_beta_2$mod.fits, FUN = function(i){
  lapply(ECxs, FUN = function(x){
  extract_ECx(i,ECx.val =x,posterior = F) 
    })
})
 
Cu_ind_estimates<-   Cu_ind_ECxs %>% unlist() %>%  data.frame()%>%rownames_to_column("lab") %>% 
  rename("value"=".") %>% 
  separate(col = lab, into = c("Model", "lab"), sep = "\\.") %>% 
  mutate(ECx=if_else(str_detect(lab, "EC_10"), "EC10", "EC1")) %>% 
  mutate(val=if_else(str_detect(lab, "lw"), "lw", 
                     if_else(str_detect(lab, "up"), "up",   "value"))) 


Cu_EC_estimates<- 
rbind(Cu_ind_estimates, Cu_avg_estimates) %>%  
  pivot_wider(names_from=val, values_from=value, id_cols=c(Model, ECx)) %>% 
  mutate_if(is.double, ~.^2)#convert back from sqrt transformation

  

save(norm_beta_2, all_endpoint_dat, all_fit_dat,Cu_EC_estimates, file = "cu_pred_vals.RData") 


#calculating averaged NEC 

norm_beta_nec <- modify_jagsMANEC(norm_beta, 
                                model.set = ("NEC"))
norm_beta_nec2<- modify_jagsMANEC(norm_beta_nec, 
                                 drop.models = c("NECHormesis", "NECsigmoidal"))

norm_beta_nec2$NEC


#Dunnets test

library(multcomp)
library(tidyverse)
library(bayesnec)
library(glmmTMB)
library(jagsNEC)
#dat_fc<- dat %>% mutate(conc=as.factor(conc)) %>% data.frame()

dunnett<- norm_beta_2$mod.fits$NEC3param$mod.dat %>% data.frame() %>% 
  mutate(x=as.factor(x))

fit.glmer<-glmmTMB(y~x, family=Beta(link = "logit", link_phi = "log"), data=dunnett)

summary(fit.glmer)
drop1(fit.glmer)
summary(glht(fit.glmer, linfct=mcp(x="Dunnett")))     


extract_ECx(norm_beta_2,ECx.val =1,posterior = F)

#EC1 = 2 (1-3)
#NSEC = 7 (3-11)

#save beta fit and endpoint and fit to a data file to be called in the rmd 

load("cu_pred_vals.RData")



