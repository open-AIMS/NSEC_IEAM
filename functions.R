# FUNCTIONS ---------------
ECxsigmoidal <- function(x, top, beta, d) {
  y.pred <- top * exp(-beta * (x)^d)
  
  return(y.pred)
}

NEC3param <- function(x, NEC, top, beta) {
  pre.index <- which(x <= NEC)
  post.index <- which(x > NEC)
  x.seq.pre <- x[pre.index]
  x.seq.post <- x[post.index]
  
  y.pred <- rep(NA, length(x))
  
  y.pred.pre <- top
  y.pred.post <- top * exp(-beta * (x.seq.post - NEC))
  
  y.pred[pre.index] <- y.pred.pre
  y.pred[post.index] <- y.pred.post
  
  return(y.pred)
}

# binomial example ------
ECxsigmoidal.binom <- function(scen.dat){
  x.vec <- seq(x.range[1], x.range[2], length=as.numeric(scen.dat["n.treatments"]))
  cols.vec <- as.character(paste("X", 1:as.numeric(scen.dat["n"]), sep=""))
  tt <- rbinom(n=as.numeric(scen.dat["n"])*length(x.vec), 
               size=as.numeric(scen.dat["trials"]),
               prob=ECxsigmoidal(x=x.vec, 
                                 top=as.numeric(scen.dat["top"]), 
                                 beta=as.numeric(scen.dat["beta"]), 
                                 d=as.numeric(scen.dat["d"]))
  )   
  out <- data.frame(matrix(tt, ncol=as.numeric(scen.dat["n"]), nrow = length(x.vec))) %>%
    mutate(x=x.vec,
           trials=as.numeric(scen.dat["trials"]),
           top=scen.dat["top"],
           beta=scen.dat["beta"],
           d=scen.dat["d"]) %>%
    pivot_longer(all_of(cols.vec), names_to = "rep.number", values_to = "suc") %>%
    data.frame()
  
  return(out)
}

NEC3param.binom <- function(scen.dat){
  x.vec <- seq(x.range[1], x.range[2], length=as.numeric(scen.dat["n.treatments"]))
  cols.vec <- as.character(paste("X", 1:as.numeric(scen.dat["n"]), sep=""))
  tt <- rbinom(n=as.numeric(scen.dat["n"])*length(x.vec), 
               size=as.numeric(scen.dat["trials"]),
               prob=NEC3param(x=x.vec, 
                              top=as.numeric(scen.dat["top"]), 
                              beta=as.numeric(scen.dat["beta"]), 
                              NEC=as.numeric(scen.dat["NEC"]))
  )   
  out <- data.frame(matrix(tt, ncol=as.numeric(scen.dat["n"]), nrow = length(x.vec))) %>%
    mutate(x=x.vec,
           trials=as.numeric(scen.dat["trials"]),
           top=scen.dat["top"],
           beta=scen.dat["beta"],
           NEC=scen.dat["NEC"]) %>%
    pivot_longer(all_of(cols.vec), names_to = "rep.number", values_to = "suc") %>%
    data.frame()
  
  return(out)
}

extract_output <- function(fit){
  if(is.null(fit$mod.fits$NEC3param)!=TRUE){
    out.NEC3param <- list(
      NEC = unlist(fit$mod.fits$NEC3param$sims.list$NEC),
      top.NEC3param = fit$mod.fits$NEC3param$sims.list$top,   
      beta.NEC3param = fit$mod.fits$NEC3param$sims.list$beta, 
      NSEC01.NEC3param =  extract_NSEC(fit$mod.fits$NEC3param, posterior = TRUE),
      NSEC05.NEC3param =  extract_NSEC(fit$mod.fits$NEC3param, posterior = TRUE, sig.val = 0.05),      
      NSEC10.NEC3param =  extract_NSEC(fit$mod.fits$NEC3param, posterior = TRUE, sig.val = 0.1),         
      ECx10.NEC3param =  extract_ECx(fit$mod.fits$NEC3param, posterior = TRUE),
      ECx05.NEC3param =  extract_ECx(fit$mod.fits$NEC3param, posterior = TRUE, ECx.val = 5),      
      ECx01.NEC3param =  extract_ECx(fit$mod.fits$NEC3param, posterior = TRUE, ECx.val = 1),  
      pred.vals.NEC3param = predict(fit$mod.fits$NEC3param))
  }else{out.NEC3param <- list()}
  if(is.null(fit$mod.fits$ECxsigmoidal)!=TRUE){
    out.ECxsigmoidal <- list(   
      d = unlist(fit$mod.fits$ECxsigmoidal$sims.list$d),
      top.ECxsigmoidal = fit$mod.fits$ECxsigmoidal$sims.list$top,
      beta.ECxsigmoidal = fit$mod.fits$ECxsigmoidal$sims.list$beta,      
      NSEC01.ECxsigmoidal =  extract_NSEC(fit$mod.fits$ECxsigmoidal, posterior = TRUE),
      NSEC05.ECxsigmoidal =  extract_NSEC(fit$mod.fits$ECxsigmoidal, posterior = TRUE, sig.val = 0.05),     
      NSEC10.ECxsigmoidal =  extract_NSEC(fit$mod.fits$ECxsigmoidal, posterior = TRUE, sig.val = 0.1), 
      ECx10.ECxsigmoidal =  extract_ECx(fit$mod.fits$ECxsigmoidal, posterior = TRUE),
      ECx05.ECxsigmoidal =  extract_ECx(fit$mod.fits$ECxsigmoidal, posterior = TRUE, ECx.val = 5),      
      ECx01.ECxsigmoidal =  extract_ECx(fit$mod.fits$ECxsigmoidal, posterior = TRUE, ECx.val = 1), 
      pred.vals.ECxsigmoidal = predict(fit$mod.fits$ECxsigmoidal)) 
  }else{out.ECxsigmoidal <- list()}
  out <- list(mod.stats = fit$mod.stats, NEC3param = out.NEC3param, ECxsigmoidal= out.ECxsigmoidal)
  return(out)
}

create_sim_dat <- function(scenario.dat){
  # binomial ECxsigmoidal
  scenario.ECx <- scenario.dat %>%
    filter(model=="Sigmoidal") %>%
    unique()
  ECx.sim.dat.binom <- apply(scenario.ECx, MARGIN=1, FUN=ECxsigmoidal.binom)
  names(ECx.sim.dat.binom) <- scenario.ECx$key
  
  # binomial NEC3param
  scenario.NEC <- scenario.dat %>%
    filter(model=="NEC") %>%
    unique()
  NEC.sim.dat.binom <- apply(scenario.NEC, MARGIN=1, FUN=NEC3param.binom)
  names(NEC.sim.dat.binom) <- scenario.NEC$key
  
  # binomial all sim.dat
  sim.dat.binom <- c(ECx.sim.dat.binom, NEC.sim.dat.binom)
  return(sim.dat.binom)
}


numericcharacters <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(x)))) & is.character(x)
}


extract_sim_dat <- function(fitted.binom.allsims, n){
  NEC.vars <- c("NEC_X2.5.", "NEC_X50.", "NEC_X97.5.",  "top.NEC3param_X2.5.", "top.NEC3param_X50.",     
             "top.NEC3param_X97.5.", "beta.NEC3param_X2.5.", "beta.NEC3param_X50.", "beta.NEC3param_X97.5.", "NSEC01.NEC3param_X2.5.", 
             "NSEC01.NEC3param_X50." , "NSEC01.NEC3param_X97.5.", "NSEC05.NEC3param_X2.5.", "NSEC05.NEC3param_X50.", 
             "NSEC05.NEC3param_X97.5.", "NSEC10.NEC3param_X2.5.", "NSEC10.NEC3param_X50.", "NSEC10.NEC3param_X97.5.",
             "ECx10.NEC3param_X2.5.", "ECx10.NEC3param_X50.", "ECx10.NEC3param_X97.5.", "ECx05.NEC3param_X2.5.", "ECx05.NEC3param_X50.",
             "ECx05.NEC3param_X97.5.", "ECx01.NEC3param_X2.5.", "ECx01.NEC3param_X50.", "ECx01.NEC3param_X97.5.")
  ECx.vars <- c("d_X2.5.", "d_X50.", "d_X97.5.",  "top.ECxsigmoidal_X2.5.", "top.ECxsigmoidal_X50.",     
              "top.ECxsigmoidal_X97.5.", "beta.ECxsigmoidal_X2.5.", "beta.ECxsigmoidal_X50.", "beta.ECxsigmoidal_X97.5.", "NSEC01.ECxsigmoidal_X2.5.", 
              "NSEC01.ECxsigmoidal_X50." , "NSEC01.ECxsigmoidal_X97.5.", "NSEC05.ECxsigmoidal_X2.5.", "NSEC05.ECxsigmoidal_X50.", 
              "NSEC05.ECxsigmoidal_X97.5.", "NSEC10.ECxsigmoidal_X2.5.", "NSEC10.ECxsigmoidal_X50.", "NSEC10.ECxsigmoidal_X97.5.",
              "ECx10.ECxsigmoidal_X2.5.", "ECx10.ECxsigmoidal_X50.", "ECx10.ECxsigmoidal_X97.5.", "ECx05.ECxsigmoidal_X2.5.", "ECx05.ECxsigmoidal_X50.",
              "ECx05.ECxsigmoidal_X97.5.", "ECx01.ECxsigmoidal_X2.5.", "ECx01.ECxsigmoidal_X50.", "ECx01.ECxsigmoidal_X97.5.")
  
# out.all <- list()  
# for(l in 1:length(fitted.binom.allsims[[n]])){
#    x <- fitted.binom.allsims[[n]][[l]]
  tt2 <- lapply(fitted.binom.allsims[[n]], FUN=function(x){      
      # NEC model parameter dat
    if(x[1]!="error"){
      if(length(x$NEC3param)>0){       
         NEC.out <- do.call("rbind", lapply(x$NEC3param[-length(x$NEC3param)], 
                                                    FUN = quantile, probs=c(0.025, 0.5, 0.975))) %>%
             data.frame() %>%
             rownames_to_column(var = "parameter") %>%
             pivot_longer(cols = c("X2.5.", "X50.", "X97.5.")) %>%
             data.frame() %>%
             unite(rowname, parameter, name)
       }else{
         NEC.out <- data.frame("rowname"=NEC.vars,"value"=rep(NA, 27))}  
      # ECx model parameter dat  
      if(length(x$ECxsigmoidal)>0){
          ECx.out <- do.call("rbind", lapply(x$ECxsigmoidal[-length(x$ECxsigmoidal)], 
                                             FUN = quantile, probs=c(0.025, 0.5, 0.975))) %>%
             data.frame() %>%
             rownames_to_column(var = "parameter") %>%
             pivot_longer(cols = c("X2.5.", "X50.", "X97.5.")) %>%
             data.frame() %>%
             unite(rowname, parameter, name)
      }else{
          ECx.out <- data.frame("rowname"=ECx.vars,"value"=rep(NA, 27))} 
      out <- rbind(NEC.out, ECx.out)
   }else{
      out <- data.frame("rowname"=c(NEC.vars, ECx.vars), "value"=rep(NA, 27*2))
   }
      
    out$sim <- n
  #out.all <- c(out.all, list(out))
  
  return(out)
  })
  return.dat <- bind_rows(tt2, .id="key")
   return(return.dat)
 }


GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}