# Copyright (c) 2021 Mathyn Vervaart, Mark Strong, Karl Claxton, Nicky Welton, Torbjørn Wisløff, Eline Aas
# Licensed under the MIT License

#####################################################################################
# Load packages 
#####################################################################################
library(flexsurv)
library(survminer)
library(MASS)
library(mgcv)
library(doParallel)
library(rstan)
library(gridExtra)
library(data.table)
library(drc)

library(compiler)
enableJIT(3)

RNGkind("L'Ecuyer-CMRG") # "L'Ecuyer-CMRG" random number generator 


#####################################################################################
# Function for generating survival times
#####################################################################################
gen_times_fun <- function(npat, 
                          shape_w, 
                          scale_w, 
                          shape_g, 
                          rate_g, 
                          k) {
  
  set.seed(1)
  tt1 <- rweibull(npat/2 * k, shape_w, scale_w)
  tt2 <- rgamma(npat/2 * k, shape_g, rate_g)
  tt <- sort(c(tt1, tt2), decreasing = FALSE)
  tt <- colMeans(matrix(tt, k))
  return(tt)
}


#####################################################################################
# Function for truncating and right-censoring a dataset
#####################################################################################
ipd_cens_fun <- function(times, 
                         tmin, 
                         tmax, 
                         treat) {
  
  tt <- subset(times, times > tmin) # select times larger than left-truncation point
  event <- ifelse(tt > tmax, 0, 1) # event indicator
  tt <- pmin(tt, tmax) # right-censor times at tmax
  trunc <- rep(tmin, length(tt)) # left-truncate times at tmin
  treat <- treat
  ipd <-  as.data.frame(cbind(tt, event, trunc, treat))
  return(ipd)
}


#####################################################################################
# Function for selecting ggplot colors
#####################################################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#####################################################################################
# Function for plotting K-M data
#####################################################################################
plot_km_fun <- function (dataset, title, y_lab, ...) {
  
  fit_obj <- do.call(survfit, list(formula = Surv(tt, event) ~ treat, data = dataset))
  
  ggsurvplot(
    fit_obj,
    combine = T,
    conf.int = F,
    title = title,
    pval = F,
    font.submain = c(font_size, "plain"),
    palette = c(col[1], col[3]),#"lancet",
    xlab = "Time",
    ylab = y_lab,
    xlim = c(0, time_current),
    ylim = c(0.7,1),
    risk.table = TRUE,
    break.time.by = 3,
    
    legend = "top", 
    legend.labs = c("New treatment", "Standard care"),
    legend.title = "",
    
    ggtheme = theme_light(),
    tables.theme = theme_cleantable() + theme(panel.border = element_blank()),
    tables.y.text = F,
    risk.table.height = 0.2,
    
    font.x = c(font_size, "plain"),
    font.y = c(font_size, "plain"),
    font.legend = c(font_size-2, "plain"),
    font.tickslab = c(font_size-2, "plain", "black"),
    fontsize = 4
  )
  
}  


#####################################################################################
# Function to extract legend from plot
#####################################################################################
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


#####################################################################################
# Function to fit survival models
#####################################################################################
survfit_fun <- function (dists, 
                         ipd, 
                         t_horizon, 
                         names_scenarios) {  
  
  # create lists to store model summaries
  surv_summ <- vector(mode = "list", length = length(ipd))
  names(surv_summ) <- names_scenarios
  surv_summ <- lapply(surv_summ, function(x)
    x = vector(mode = "list", length = 4))
  for (i in 1:length(ipd)) {
    names(surv_summ[[i]]) <- c("fit", "aic", "Aw", "mean_survival")
  }
  
  for (i in 1:length(ipd)) {
    
    surv_summ[[i]]$fit <- vector(mode = "list", length = length(dists)) # list of models fits
    surv_summ[[i]]$aic <- matrix(NaN, length(dists), 1) # vectors of aic scores
    surv_summ[[i]]$Aw <- matrix(NaN, length(dists), 1) # vectors of aic weights
    surv_summ[[i]]$mean_survival <- matrix(NaN, length(dists), 1) # mean survival estimates
    surv_summ[[i]]$events <- matrix(NaN, length(dists), 1) # event numbers
    surv_summ[[i]]$trisk <- matrix(NaN, length(dists), 1) # time at risk
    
    names(surv_summ[[i]]$fit) <-
      names(surv_summ[[i]]$aic) <-
      names(surv_summ[[i]]$Aw) <-
      names(surv_summ[[i]]$mean_survival) <-
      names(surv_summ[[i]]$events) <-
      names(surv_summ[[i]]$trisk) <- c(dists)
  }  
  
  # fit the models
  for (j in 1:length(ipd)) {
    
    for (i in 1:length(dists)) {
      
      if(dists[i] == "exp" || dists[i] == "gengamma") {
        
        model <- flexsurvreg(Surv(tt, event) ~ 1, data = ipd[[j]], 
                             dist = dists[i], method = "BFGS")
        
      } else {
        
        model <- flexsurvreg(Surv(tt, event) ~ 1, data = ipd[[j]], 
                             dist = dists[i], method = "Nelder-Mead")
        
      }
      
      surv_summ[[j]]$fit[[i]] <- cbind(model$coefficients, model$cov)
      
      surv_summ[[j]]$aic[[i]] <- model$AIC
      
      surv_summ[[j]]$mean_survival[[i]] <- 
        eval(parse(text = paste0("rmst_", dists[i], "(", t_horizon, ",", 
                                 paste0(as.vector(model$res[,1]), collapse=","), ")")  )) 
      
      surv_summ[[j]]$events[[i]] <- model$events
      
      surv_summ[[j]]$trisk[[i]] <- model$trisk
      
    }
  }
  
  # compute AIC-based model weights
  aic_min <- lapply(surv_summ, function (x) {min(x$aic)} ) # identify the lowest AIC for each dataset
  
  # transform the AIC differences back to the scale of probabilities
  Ak <- lapply(1:length(ipd), 
               function (x) {exp(-0.5*(surv_summ[[x]]$aic-aic_min[[x]]))}) 
  # compute weights
  for (i in 1:length(ipd)) {
    surv_summ[[i]]$Aw <- apply(Ak[[i]], 2, function (y) {
      Ak[[i]] / sum(y)
    })
  }
  
  for (i in 1:length(ipd)) {
    names(surv_summ[[i]]$Aw) <- c(dists)
  }  
  
  # function output
  return(surv_summ)
}


#####################################################################################
# Function for sampling from a prior multivariate normal distribution
#####################################################################################
mvrnorm_fun <- function (surv_summ, 
                         dists, 
                         t_horizon, 
                         n, 
                         seed) {
  
  # set the seed
  set.seed(seed)
  
  # sample K survival distributions
  dist_select <- sample(x = dists, n, replace = T, prob = surv_summ$Aw)
  
  # matrix of transformed prior mean vectors
  mu <- lapply(dist_select, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,1]"))) 
  })
  
  # matrix of transformed prior variance matrices
  cov_matrix <- lapply(dist_select, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,2:","ncol(surv_summ$fit$", x,")]")))
  })
  
  # sample from priors on transformed scale
  prior_samples_trans <- lapply(1:n, function (x) {
    set.seed(seed + x)
    mvrnorm(1, mu = mu[[x]], Sigma = cov_matrix[[x]])
  })
  
  # transform parameter samples back to original scale
  prior_samples <- lapply(1:n, function (x) {
    sapply(1:length(prior_samples_trans[[x]]), function (y)  {
      eval(parse(text = paste0("flexsurv.dists$", dist_select[[x]],
                               "$inv.transforms[[y]](prior_samples_trans[[x]][y])")))
    })
  })
  
  prior_surv <- lapply(1:n, function (x) {
    eval(parse(text = paste0("rmst_", dist_select[[x]], "(", t_horizon, 
                             ",", paste(prior_samples[[x]], collapse = ",", sep = ","), ")")))
  })
  
  # return output
  output <- list(distribution = dist_select, theta = prior_samples, mean_surv = prior_surv)
  return(output)
}


#####################################################################################
# function for sampling from prior and generating datasets
#####################################################################################
gen_datasets_fun <- function(surv_summ, 
                             dists, 
                             time_current, 
                             t_new, 
                             t_trunc, 
                             loops_outer, 
                             seed) {
  
  # set the seed
  set.seed(seed)
  
  # sample survival distributions
  prior_mod_samples <- sample(x = dists, loops_outer, replace = T, prob = surv_summ$Aw)
  
  # list of transformed prior parameter vectors
  mu <- lapply(prior_mod_samples, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,1]"))) 
  })
  
  # list of transformed prior variance matrices
  cov_matrix <- lapply(prior_mod_samples, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,2:","ncol(surv_summ$fit$", x,")]")))
  })
  
  # sample from priors on transformed scale
  prior_par_samples_trans <- lapply(1:loops_outer, function (x) {
    set.seed(seed + x)
    mvrnorm(1, mu = mu[[x]], Sigma = cov_matrix[[x]])
  })
  
  # transform parameter samples back to original scale
  prior_par_samples <- lapply(1:loops_outer, function (x) {
    sapply(1:length(prior_par_samples_trans[[x]]), function (y)  {
      eval(parse(text = paste0("flexsurv.dists$", 
                               prior_mod_samples[[x]],"$inv.transforms[[y]](prior_par_samples_trans[[x]][y])")))
    })
  })
  
  # number of patients left at risk at the end of current follow up
  at_risk <- length(t_trunc) 
  
  # identify truncated cumulative density corresponding to current follow up time for each prior draw
  cond_cdf <- sapply(1:loops_outer, function (x) {
    pdist_parse <- parse(text = paste0("p", prior_mod_samples[[x]],
                                       "(", "y", ",", paste(prior_par_samples[[x]], collapse=", ", sep="") ,")"    ))
    pdist <- function (y) {eval(pdist_parse)}
    pdist(t_trunc)
  })
  cond_cdf <- t(cond_cdf)
  
  # generate random number sequence for conditional CDF for each prior draw
  rand_num <- apply(cond_cdf, 2, function (x) {
    runif(loops_outer, x,  1)
  })
  
  # plug the random numbers into the quantile function given each prior draw
  tt <- sapply(1:loops_outer, function (x) {
    qdist_parse <- parse(text = paste0("q", prior_mod_samples[[x]],
                                       "(", "y", ",", paste(prior_par_samples[[x]], collapse=", ", sep="") ,")"    ))
    qdist <- function (y) {eval(qdist_parse)}
    qdist(rand_num[x,])
  })
  tt <- t(tt)
  
  # event indicator = 0 if generated survival time > future follow up time, 1 otherwise
  event <- apply(tt, 2, function(x) {
    as.numeric(ifelse(x > t_new, 0, 1))
  })
  
  # right-censor the survival times
  tt <- pmin(tt, t_new)
  
  # create list of dataframes with survival datasets
  sim_datasets <- lapply(1:loops_outer, function(x) as.data.frame(matrix(, at_risk, 2)))
  for (i in 1:loops_outer) {
    sim_datasets[[i]] <- cbind(tt[i,], event[i,])
  }
  
  # return output
  sim_list <- list(prior_par_samples = prior_par_samples, prior_mod_samples = prior_mod_samples, sim_datasets = sim_datasets)
  return(sim_list)
}  


#####################################################################################
# function to compute EVSI using GAM regression
#####################################################################################
gam_evsi_fun <- function (outer_loop_tr1,
                          outer_loop_tr2,
                          ipd_tr1,
                          ipd_tr2,
                          time_add,
                          scenarios,
                          names_scenarios,
                          t_horizon) {
  
  # list to store evsi results
  evsi <-  matrix(NaN, length(time_add), scenarios)
  se <-  matrix(NaN, length(time_add), scenarios)
  
  for (s in 1:scenarios) {
    
    
    for (t in 1:length(time_add)) {
      
      print(paste("Computing EVSI for the", names_scenarios[s], "dataset", 
                  "for an additional follow-up of", time_add[t], "months."))
      
      # datasets generated in outer loop
      datasets_tr1 <- outer_loop_tr1[[s]][[t]]$sim_datasets 
      datasets_tr2 <- outer_loop_tr2[[s]][[t]]$sim_datasets 
      
      # prior model samples
      prior_mod_tr1 <- outer_loop_tr1[[s]][[t]]$prior_mod_samples
      prior_mod_tr2 <- outer_loop_tr2[[s]][[t]]$prior_mod_samples
      
      # prior parameter samples
      prior_par_samples_tr1 <- outer_loop_tr1[[s]][[t]]$prior_par_samples
      prior_par_samples_tr2 <- outer_loop_tr2[[s]][[t]]$prior_par_samples
      
      # left truncation times
      temp_LL_tr1 <- as.integer(subset(ipd_tr1[[s]][, 1], ipd_tr1[[s]][, 2] == 0)) 
      temp_LL_tr2 <- as.integer(subset(ipd_tr2[[s]][, 1], ipd_tr2[[s]][, 2] == 0)) 
      
      # prior survival
      nb_prior_tr1 <- sapply(1:length(prior_mod_tr1), function (x) {
        eval(parse(text = paste0("rmst_", prior_mod_tr1[[x]], "(", t_horizon, ",", 
                                 paste(prior_par_samples_tr1[[x]], collapse = ",", sep = ","), ")")))
      })
      
      nb_prior_tr2 <- sapply(1:length(prior_mod_tr2), function (x) {
        eval(parse(text = paste0("rmst_", prior_mod_tr2[[x]], "(", t_horizon, ",", 
                                 paste(prior_par_samples_tr2[[x]], collapse = ",", sep = ","), ")")))
      })
      
      
      ### obtain summary statistic (number of observed events and time at risk)
      summ_stat_tr1 <- sapply(datasets_tr1, function (x) {
        par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(temp_LL_tr1)))
      })
      
      summ_stat_tr2 <- sapply(datasets_tr2, function (x) {
        par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(temp_LL_tr2)))
      })
      
      
      ### GAM model for predicting posterior net benefit

      # fit GAM models
      gam_mod_tr1 <- gam(nb_prior_tr1 ~ te(summ_stat_tr1[1,], summ_stat_tr1[2,]))
      
      gam_mod_tr2 <- gam(nb_prior_tr2 ~ te(summ_stat_tr2[1,], summ_stat_tr2[2,]))
      
      # extract fitted values                  
      g_hat_tr1 <- gam_mod_tr1$fitted
      g_hat_tr2 <- gam_mod_tr2$fitted
      
      # compute the max of the means
      max_mean <- pmax(g_hat_tr1, g_hat_tr2)
      
      # compute evsi
      evsi[t,s] <- mean(max_mean) - max(mean(g_hat_tr1), mean(g_hat_tr2)) 
      
      
      ### compute standard error of the GAM estimator
      
      # extract the basis function values
      Xstar_tr1 <- model.matrix(gam_mod_tr1)
      Xstar_tr2 <- model.matrix(gam_mod_tr2)
      
      # extract coefficients
      beta_tr1 <- gam_mod_tr1$coef
      beta_tr2 <- gam_mod_tr2$coef
      
      # compute fitted values
      fits_tr1 <- Xstar_tr1 %*% beta_tr1
      fits_tr2 <- Xstar_tr2 %*% beta_tr2
      
      # covariance matrix
      v_tr1 <- vcov(gam_mod_tr1) 
      v_tr2 <- vcov(gam_mod_tr2) 
      
      # sample from the parameter distributions
      set.seed(123)
      parameter_draws_tr1 <- mvrnorm(10000, beta_tr1, v_tr1)
      
      set.seed(456)
      parameter_draws_tr2 <- mvrnorm(10000, beta_tr2, v_tr2)
      
      # from these draws, calculate draws from fitted values
      fitted_draws_tr1 <- Xstar_tr1 %*% t(parameter_draws_tr1)
      fitted_draws_tr2 <- Xstar_tr2 %*% t(parameter_draws_tr2)
      
      # compute EVSI for each sample
      evsi_samples <-
        sapply(1:ncol(fitted_draws_tr1), function (x) {
          v_max_mean <- pmax(fitted_draws_tr1[, x], fitted_draws_tr2[, x])
          evsi_samples <-
            mean(v_max_mean) - max(mean(fitted_draws_tr1[, x]), mean(fitted_draws_tr2[, x]))
          return(evsi_samples)
        })
      
      # compute sample standard deviation
      se[t,s]  <- sd(evsi_samples)
      
      print(paste("EVSI is", round(evsi[t, s], 2), "(", round(se[t, s], 3), ")",
                  "for an additional follow-up of", time_add[t], "months for the", 
                  names_scenarios[s], "dataset."
      )
      )
      
    }
  }
  
  # return function output
  output <- list()
  output$evsi <- evsi
  output$se <- se
  return(output)
}


#####################################################################################
# function for generating posterior samples using mcmc
#####################################################################################
gen_mcmc_post <-
  function(sim_datasets,
           smod,
           ipd,
           prior,
           time_add,
           dists,
           dist_temp,
           rmst,
           t_horizon,
           mcmc_length,
           warmup,
           chains,
           seed) {
    
    # initialize variables
    theta0 <- prior[, 1] # prior means
    sigma0 <- prior[, 2:ncol(prior)] # prior variance-covariance matrix
    LL <- as.integer(subset(ipd[, 1], ipd[, 2] == 0)) # left truncation time
    K <- as.integer(nrow(prior))  # number of model parameters
    loops_outer <- length(sim_datasets) # number of outer loops
    
    # register parallel backend
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    
    # export global variables and functions to each cluster worker
    ex <- list(setdiff(ls(),
                       c("sim_datasets",
                         "outer_loop_tr1",
                         "outer_loop_tr2",
                         "m_prior_tr1",
                         "m_prior_tr2",
                         "post_list_tr1",
                         "post_list_tr2"
                       )
    )) #ls(globalenv())
    clusterExport(cl, ex)
    
    # save the generated datasets to a rds file, which will be loaded on each worker within the parallel loop
    saveRDS(sim_datasets, file = here("sim_datasets.rds"))
    
    # sample from Stan models in parallel
    fits <-
      foreach(
        i = 1:loops_outer,
        .combine = cbind,
        .packages = c("here", "rstan", "bridgesampling")
        
      )  %dopar%    {
        
        # load the generated datasets on each worker
        sim_datasets <- readRDS(here("sim_datasets.rds"))
        
        # prepare data for Stan
        n <- as.integer(length(sim_datasets[[i]][,1])) # total number at risk
        t <- as.vector(sim_datasets[[i]][,1]) # times
        d <- as.vector(sim_datasets[[i]][,2]) # censoring indicator
        trisk0 <- as.integer(sum(ipd[,1]))  # prior time at risk
        events0 <- as.integer(sum(ipd[,2])) # prior number of events
        
        # data for gamma model
        nobs <- as.integer(length(subset(sim_datasets[[i]][,1], sim_datasets[[i]][,2]==1))) # number of observed events
        ncens <- as.integer(length(subset(sim_datasets[[i]][,1], sim_datasets[[i]][,2]==0))) # number of censored events
        yobs <- as.array(subset(sim_datasets[[i]][,1], sim_datasets[[i]][,2]==1)) # observed events
        ycens <- as.array(subset(sim_datasets[[i]][,1], sim_datasets[[i]][,2]==0)) # censored events
        
        # list of data for Stan
        surv_data <- list(
          n = n,
          t = t,
          d = d,
          nobs = nobs,
          ncens = ncens,
          yobs = yobs,
          ycens = ycens,
          LL = LL,
          K = K,
          theta0 = theta0,
          sigma0 = sigma0,
          trisk0 = trisk0,
          events0 = events0
        )
        
        # sample from the Stan model
        stan_fit <-
          sampling(
            smod,
            data = surv_data,
            cores = 1,
            chains = chains,
            iter = mcmc_length,
            seed = seed,
            warmup = warmup
          )
        
        if (length(dists) == 1) {
          logml <- -1
          
        } else {
          
          # bridge sampling estimate of the log marginal likelihood
          set.seed(seed + 1)
          log_marg_like <-
            bridge_sampler(stan_fit, method = "warp3", silent = TRUE)
          logml <- log_marg_like$logml
        }
        
        # extract posterior samples
        fit <- sapply(1:nrow(prior), function (j) {
          eval(parse(
            text = paste(
              "extract(stan_fit, par = c('",
              rownames(prior)[j],
              "'), permuted = TRUE)$",
              rownames(prior)[j],
              sep = ""
            )
          ))
        })
        
        # compute posterior expected survival
        fit <- mean(apply(fit, 1, function (k) {
          if (dist_temp == "gengamma") {
            rmst(t_horizon, k[1], k[2], k[3])
          }
          if (dist_temp == "exp") {
            rmst(t_horizon, k[1])
          } else {
            rmst(t_horizon, k[1], k[2])
          }
          
        }))
        
        return(list(fit, logml))
        
      } # close "dopar" loop
    
    # stop cluster
    stopCluster(cl)
    
    # extract and return return function output
    odd_indexes <- seq(1,loops_outer*2,2)
    even_indexes <- seq(2,loops_outer*2,2)
    post_samples <- list()
    post_samples$post_fit <- lapply(odd_indexes, function (x) fits[[x]])
    post_samples$logml <- lapply(even_indexes, function (x) fits[[x]])
    return(post_samples)
  }


#####################################################################################
# function to compute posterior model weights
#####################################################################################
update_weights_fun <-
  function (post_list,
            surv_summ,
            scenarios,
            time_add,
            dists,
            loops_outer) {
    
    lapply(1:scenarios, function (s) {
      lapply(1:length(time_add), function (t) {
        if (length(dists) == 1) {
          return(matrix(1, 1, loops_outer))
          
        } else {
          # combine the estimates of the log marginal likelihoods into a matrix
          m_logml <- sapply(post_list, function (x) {
            unlist(x[[s]][[t]]$logml)
          })
          
          # transform the estimates of the log marginal likelihoods to the probability scale
          logml_k <- apply(m_logml, 2, function (x) {
            exp(x)
          })
          
          # compute posterior weights
          post_weights_k <- sapply(1:ncol(logml_k), function (x) {
            surv_summ[[s]]$Aw[[x]] * logml_k[, x]
          })
          
          post_weights <- t(apply(post_weights_k, 1, function (x) {
            x / sum(x)
          }))
          
          if (nrow(post_weights) > 1)  {
            colnames(post_weights) <- dists
          }
          
          return(post_weights)
        }
      })
    })
  }


#####################################################################################
# function to multiply posterior expected survival with posterior model weights
#####################################################################################
w_post_ex_surv_fun <- function (scenarios, 
                                time_add, 
                                dists, 
                                post_surv, 
                                post_weights) {
  
  lapply(1:scenarios, function (s) {
    
    lapply(1:length(time_add), function (t) {
      
      sapply(1:length(dists), function (d) {
        
        unlist(post_surv[[d]][[s]][[t]]) * post_weights[[s]][[t]][, d]
        
      })
    })
  })
}


#####################################################################################
# function to compute posterior model-averaged expected survival
#####################################################################################
post_ma_ex_surv_fun <- function (scenarios, 
                                 time_add, 
                                 w_post_surv) {
  
  lapply(1:scenarios, function (s) {
    
    lapply(1:length(time_add), function (t) {
      
      rowSums(w_post_surv[[s]][[t]])
      
    })
  })
}


#####################################################################################
# function for computing marginal EVSI and cost functions
#####################################################################################
marg_fun <- function (evsi, time_add, n_pat, trial_cost, incr_nb) {
  
  df <- as.data.frame(cbind("y" = evsi, "x" = time_add))
  
  # interpolate EVSI using asymptotic regression
  ar_mod_evsi <- drm(y ~ x, data = df, fct = AR.3(fixed = c(NA, NA, NA)), pmodels = list(~1, ~1,~1))
  new_times <- data.frame("x" = seq.int(min(time_add), max(time_add), 1))
  evsi_smooth <- predict(ar_mod_evsi, newdata = new_times)
  
  # plot marginal EVSI + marginal cost
  x_times <- seq.int(min(time_add)+1, max(time_add), 1)
  
  marg_evsi <- as.numeric(diff(evsi_smooth))
  n_times <- seq.int(1, max(time_add)-min(time_add), 1)
  n_times_desc <- sort(n_times, decreasing = T)
  marg_evsi_pop <- cbind(x_times, (n_times_desc*marg_evsi)*n_pat) 
  marg_evsi_pop <- data.table(times = x_times, NB = (n_times_desc*marg_evsi)*n_pat, group = 1) 
  
  # marginal trial costs
  marg_cost_trial <- data.table(times = x_times, NB = rep(trial_cost, length(marg_evsi)), group = 2) 
  
  # marginal costs due to withholding access
  marg_cost_wait  <- new_times$x*(incr_nb*n_pat)
  marg_cost_wait <- diff(marg_cost_wait)
  marg_cost_wait <- data.table(times = x_times, NB = (marg_cost_wait + marg_cost_trial$NB), group = 3)
  
  df <- rbind(marg_evsi_pop, marg_cost_trial, marg_cost_wait)
  df$group <- as.factor(df$group)
  
  return(df)
}


#####################################################################################
# Function for finding the intersect between two vectors
#####################################################################################
intersect_fun <- function (x, y, x2, y2, step) {
  v1 <- approx(x, y, xout = seq.int(min(x), max(x), step))
  v2 <- approx(x2, y2, xout = seq.int(min(x), max(x), step))
  
  intersect <- v1$x[which.min(abs(v1$y - v2$y))]
  return(intersect)
}


#####################################################################################
# ENBS plotting function
#####################################################################################
plot_enbs_fun <- function (df, ylab, leg_pos, x1, x2, y1, y2, anno1, anno2) {
  
  ggplot(data=subset(df, df$times >= 6 & df$times <=84), aes(x=times, y=NB, group = group)) +
    theme_light() + 
    geom_line(aes(color=group), size = 0.5) +
    xlab("Additional follow-up (months)") +
    ylab(ylab) +
    theme(legend.position = leg_pos) +
    ggtitle("") + scale_size(range=c(0.1, 2), guide=FALSE) + 
    scale_x_continuous(breaks = seq(12,84,24), limits = c(6,84)) +
    scale_y_continuous(breaks = seq(0,200,50), limits = c(0,210)) +
    
    geom_point(aes(x = x1, y = y1), colour="black", size = 1) +
    geom_point(aes(x = x2, y = y2), colour="black", size = 1) +
    annotate("text", label = anno1, x = x1 + 4, y = y1+12, size = 4) +
    annotate("text", label = anno2, x = x2 + 4, y = y2+12, size = 4) +
    
    theme(axis.text.x = element_text(color="black")) +
    theme(axis.text.y = element_text(color="black")) +
    theme(text = element_text(size = 14))
}