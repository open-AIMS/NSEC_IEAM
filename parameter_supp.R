require(tidyverse)

fig_dir <- "C:/Users/rfisher/Australian Institute of Marine Science/Program 2.3 - TERA - General/Projects/3158 Statistical Ecotoxicology/Estimating no effect/MS/revision1/"


all_mod_fits <- list(
  `Case study 1 _ C. armigera` = norm_beta_2$mod.fits,
  `Case study 2 _ A. amphitrite` = out.aa.all$mod.fits,
  `Case study 2 _ A. millepora` = Coral.out4.all$mod.fits,
  `Case study 2 _ N. dorsatus` = out.nd.all$mod.fits,
  `Case study 2 _ S. variolaris` = Urchin.out2$mod.fits  
)


all_params_out <- lapply(all_mod_fits, FUN = function(l){
  lapply(l, FUN = function(x){
      x$summary |> 
      data.frame()|> 
      rownames_to_column(var = "parameter") 
    }) |> 
    bind_rows(.id="model")  
}) |> 
  bind_rows(.id="species") |> 
  dplyr::filter(!grepl("SS", parameter)) |> 
  dplyr::filter(!grepl("deviance", parameter))

head(all_params_out)

write.csv(all_params_out, 
          file = paste(fig_dir,"Estimating_no_effect_Revision1_supp2.csv", sep=""),
          row.names = FALSE)
