# Copyright (c) 2021 Mathyn Vervaart, Mark Strong, Karl Claxton, Nicky Welton, Torbjørn Wisløff, Eline Aas
# Licensed under the MIT License

#rm(list=ls(all=TRUE))	# clear workspace

#####################################################################################
# Load packages and functions
#####################################################################################

# auto-compile functions
library(compiler)
enableJIT(3) 

# load functions
library(here)
source(here("analysis", "functions.R"))

library(rstan)


#####################################################################################
# Analysis settings
#####################################################################################

# number of patients enrolled in each trial arm
npat_tr1 <- 200
npat_tr2 <- 200

# current follow-up time
time_current <- 12

# additional follow-up times
time_add <- c(12, 24, 36, 48)

# scenario names (number of scenarios must match the number of scenario names)
names_scenarios <- c("Increasing hazard", "Decreasing hazard") 
scenarios <- length(names_scenarios)

# Weibull shapes and scales for simulating the hypothetical datasets
shape_w <- c(1.1, 0.6)
scale_w_tr1 <- c(70, 80)
scale_w_tr2 <- c(50, 57)

# gamma shapes and rates for simulating the hypothetical datasets
shape_g <- c(1.8, 0.8) 
rate_g_tr1 <- c(0.04, 0.01)
rate_g_tr2 <- c(0.04, 0.01)

# to generate smooth datasets based on evenly spaced quantiles of the distribution
k <- 10e3

# names of survival distributions to be used in the analysis
# names need to match the distribution names used by the "flexsurv" package
# currently only the exponential, Weibull, gamma, lognormal and log-logistic models are supported for the MCMC approach
# change "dists <- weibull" to perform the Weibull EVSI analysis in the paper
dists <- c("weibull", "gamma", "lnorm", "llogis") 

# time horizon in months for computing expected survival
t_horizon <- 240

# number of datasets to be generated in the outer loop of the EVSI procedure
loops_outer <- 6e3

# number of prior draws for computing EVPI (not used in EVSI calculations)
n_draws <- 100e3

# MCMC settings
mcmc_length <- 2000 # number of posterior samples per chain
warmup <- 1000 # number of posterior samples to be discarded as warm-up per chain
chains <- 2 # number of chains


#####################################################################################
# Generate IPD and right-censor data at current follow-up time
#####################################################################################

times_tr1 <- mapply(gen_times_fun, npat_tr1, shape_w, scale_w_tr1, shape_g, rate_g_tr1, k)
ipd_tr1 <- apply(times_tr1, 2, ipd_cens_fun,  tmin = 0, tmax = time_current, treat = 1)

times_tr2 <- mapply(gen_times_fun, npat_tr2, shape_w, scale_w_tr2, shape_g, rate_g_tr2, k)
ipd_tr2 <- apply(times_tr2, 2, ipd_cens_fun,  tmin = 0, tmax = time_current, treat = 2)

# merge datasets
ipd <- vector(mode = "list", length = scenarios)
for (i in 1:scenarios) {
  ipd[[i]] <- rbind(ipd_tr1[[i]], ipd_tr2[[i]])
}

names(ipd) <- names(ipd_tr1) <- names(ipd_tr2) <- names_scenarios


#####################################################################################
# Visualize K-M data
#####################################################################################

font_size <- 14
col <- gg_color_hue(3) 

# plot increasing hazard dataset
plot1 <- plot_km_fun(ipd[[1]], title = "", y_lab = "Survival probability", font_size = font_size, risktab = TRUE)
mylegend <- g_legend(plot1$plot) # extract legend
plot1$plot <- plot1$plot + theme(legend.position = "none") # remove legend in plot

# plot decreasing hazard dataset
plot2 <- plot_km_fun(ipd[[2]], title = "", y_lab = "", font_size = font_size, risktab = TRUE)
plot2$plot <- plot2$plot + theme(legend.position = "none") # remove legend in plot

# arrange plots in a grid
gA <- ggplotGrob(plot1$plot)
gB <- ggplotGrob(plot1$table)

gC <- ggplotGrob(plot2$plot)
gD <- ggplotGrob(plot2$table)


maxWidth <- grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
gD$widths[2:5] <- as.list(maxWidth)

m_layout <- rbind(c(1,1), 
                  do.call("rbind", replicate(8, c(2,3), simplify = FALSE)),
                  c(4,5),
                  c(4,5))

grid.arrange(arrangeGrob(mylegend), gA, gC, gB, gD, layout_matrix = m_layout)


#####################################################################################
# Fit survival models using MLE and obtain summaries
#####################################################################################

surv_summ_tr1 <- survfit_fun(dists, ipd_tr1, t_horizon, names_scenarios)
surv_summ_tr2 <- survfit_fun(dists, ipd_tr2, t_horizon, names_scenarios)


#####################################################################################
# Compute prior model averaged expected survival
#####################################################################################

prior_ma_ex_surv_tr1 <- lapply(surv_summ_tr1, function (x) {
  colSums(x$Aw*x$mean_surv)
})
  
prior_ma_ex_surv_tr2 <- lapply(surv_summ_tr2, function (x) {
  colSums(x$Aw*x$mean_surv)
})

# incremental model averaged expected survival
prior_ma_ex_surv_tr1[[1]][1] - prior_ma_ex_surv_tr2[[1]][1]
prior_ma_ex_surv_tr1[[2]][1] - prior_ma_ex_surv_tr2[[2]][1]


#####################################################################################
# Compute EVPI 
#####################################################################################

m_prior_tr1 <- lapply(surv_summ_tr1, function (x) {
  mvrnorm_fun(x, dists, t_horizon, n_draws, seed = 123)
})
m_prior_tr2 <- lapply(surv_summ_tr2, function (x) {
  mvrnorm_fun(x, dists, t_horizon, n_draws, seed = 456)
})

# prior survival
prior_surv_tr1 <- lapply(m_prior_tr1, function (x) unlist(x$mean_surv))
prior_surv_tr2 <- lapply(m_prior_tr2, function (x) unlist(x$mean_surv))

# compute the maximum of the prior model averaged expected survival
prior_max <- lapply(1:scenarios, function (x) {
  pmax(prior_surv_tr1[[x]], prior_surv_tr2[[x]])
})

# compute EVPI
evpi <- lapply(1:scenarios, function (x) {
  mean(prior_max[[x]]) - max(mean(prior_surv_tr1[[x]]), mean(prior_surv_tr2[[x]]))
})
rm(m_prior_tr1, m_prior_tr2, prior_surv_tr1, prior_surv_tr2, prior_max)

evpi


#####################################################################################
# Outer loop - generate datasets
#####################################################################################

# create lists to store the generated datasets
list_temp <- vector(mode = "list", length = scenarios)
list_temp <- lapply(list_temp, function(x)
  x = vector(mode = "list", length = length(time_add)))
names(list_temp) <- names_scenarios
for ( i in 1:(scenarios)) {names(list_temp[[i]])  <- c(time_current + time_add) }

outer_loop_tr1 <- outer_loop_tr2 <- list_temp

# obtain truncation times
t_trunc_tr1 <- lapply(ipd, function(x)
  subset(subset(x, x$treat == 1)$tt, subset(x, x$treat == 1)$event == 0))
t_trunc_tr2 <- lapply(ipd, function(x)
  subset(subset(x, x$treat == 2)$tt, subset(x, x$treat == 2)$event == 0))

# generate the datasets
ptm_datasets <- Sys.time() # start the clock
for (s in 1:scenarios) {
  
  for (t in 1:length(time_add)) {
    
    outer_loop_tr1[[s]][[t]] <-
      gen_datasets_fun(
        surv_summ_tr1[[s]],
        dists,
        time_current,
        time_current + time_add[[t]],
        t_trunc_tr1[[s]],
        loops_outer,
        seed = 123
      )
    
    outer_loop_tr2[[s]][[t]] <-
      gen_datasets_fun(
        surv_summ_tr2[[s]],
        dists,
        time_current,
        time_current + time_add[[t]],
        t_trunc_tr2[[s]],
        loops_outer,
        seed = 456
      )
  }
}
ptm_datasets <- Sys.time() - ptm_datasets # Stop the clock
ptm_datasets


#####################################################################################
# GAM approach to EVSI
#####################################################################################

ptm_gam <- Sys.time() # Start the clock

evsi_gam <-
  gam_evsi_fun(
    outer_loop_tr1,
    outer_loop_tr2,
    ipd_tr1,
    ipd_tr2,
    time_add,
    scenarios,
    names_scenarios,
    t_horizon
  )

ptm_gam <- Sys.time() - ptm_gam # Stop the clock
ptm_gam

evsi_gam


#####################################################################################
# MCMC approach to EVSI 
#####################################################################################

# create lists to store the posterior samples
post_list_tr1 <- post_list_tr2 <- rep(list(list_temp), length(dists))
names(post_list_tr1) <- names(post_list_tr2) <- dists
temp_tr1 <- temp_tr2 <- post_list_tr1

# start the mcmc procedure
ptm_mcmc <- Sys.time() # Start the clock

for (d in 1:length(dists)) {
  
  # select restricted mean survival function
  dist_temp <- dists[d]
  rmst_temp <- eval(parse(text = paste("rmst_", dist_temp, sep = "")))
  rmst_temp <- Vectorize(rmst_temp)
  
  # fit the model in Stan
  options(mc.cores = 1)
  rstan_options("auto_write" = TRUE)
  smod <- stan_model(here(paste(dists[d], "_model.stan", sep = "")))
  
  # sample from the Stan model
  for (s in 1:scenarios) {
    
    for (t in 1:length(time_add)) {
      
      print(paste0("Sampling for model '", dists[d], "', ", names_scenarios[s], 
                   " dataset, additional follow-up time of ", time_add[t], " months.", collapse = ""))
      
      current_time <- Sys.time()
      print(paste0("Current run time is " , round(difftime(current_time, ptm_mcmc, units='mins'), 1), " minutes.", collapse = ""))
      
      post_list_tr1[[d]][[s]][[t]] <- gen_mcmc_post(
        sim_datasets = outer_loop_tr1[[s]][[t]]$sim_datasets,
        smod = smod,
        ipd = subset(ipd[[s]], ipd[[s]]$treat == 1),
        prior = eval(parse(text = paste("surv_summ_tr1[[s]]$fit$", dists[d], sep = ""))),
        time_add = time_add[t],
        dists = dists,
        dist_temp = dist_temp,
        rmst = rmst_temp,
        t_horizon = t_horizon,
        mcmc_length = mcmc_length,
        warmup = warmup,
        chains = chains,
        seed = 123
      )
      
      post_list_tr2[[d]][[s]][[t]] <- gen_mcmc_post(
        sim_datasets = outer_loop_tr2[[s]][[t]]$sim_datasets,
        smod = smod,
        ipd = subset(ipd[[s]], ipd[[s]]$treat == 2),
        prior = eval(parse(text = paste("surv_summ_tr2[[s]]$fit$", dists[d], sep = ""))),
        time_add = time_add[t],
        dists = dists,
        dist_temp = dist_temp,
        rmst = rmst_temp,
        t_horizon = t_horizon,
        mcmc_length = mcmc_length,
        warmup = warmup,
        chains = chains,
        seed = 456
      )
    } # close "t" loop
    
  } # close "s" loop
  
} # close "d" loop

ptm_mcmc <- Sys.time() - ptm_mcmc # Stop the clock
ptm_mcmc

# extract posterior expected survival for each model
post_surv_tr1 <- post_surv_tr2 <- rep(list(list_temp), length(dists))
names(post_surv_tr1) <- names(post_surv_tr2) <- dists

for (d in 1:length(dists)) {
  
  for (s in 1:scenarios) {
    
    for (t in 1:length(time_add)) {
      
      post_surv_tr1[[d]][[s]][[t]] <- post_list_tr1[[d]][[s]][[t]]$post_fit
      
      post_surv_tr2[[d]][[s]][[t]] <- post_list_tr2[[d]][[s]][[t]]$post_fit
      
    } # close "t" loop
    
  } # close "s" loop
  
} # close "d" loop

# Update the model weights 
post_weights_tr1 <- update_weights_fun(post_list_tr1, surv_summ_tr1, scenarios, time_add, dists, loops_outer)
post_weights_tr2 <- update_weights_fun(post_list_tr2, surv_summ_tr2, scenarios, time_add, dists, loops_outer)

# multiply posterior expected survival with posterior weights for each model
w_post_ex_surv_tr1 <- w_post_ex_surv_fun(scenarios, time_add, dists, post_surv_tr1, post_weights_tr1)
w_post_ex_surv_tr2 <- w_post_ex_surv_fun(scenarios, time_add, dists, post_surv_tr2, post_weights_tr2)

# compute posterior model averaged expected survival
post_ma_ex_surv_tr1 <- post_ma_ex_surv_fun(scenarios, time_add, w_post_ex_surv_tr1)
post_ma_ex_surv_tr2 <- post_ma_ex_surv_fun(scenarios, time_add, w_post_ex_surv_tr2)

# name lists
names(post_ma_ex_surv_tr1) <-
  names(post_ma_ex_surv_tr2) <-
  names(post_weights_tr1) <-
  names(post_weights_tr2) <- names_scenarios

for (i in 1:(scenarios)) {
  names(post_ma_ex_surv_tr1[[i]])  <-
    names(post_ma_ex_surv_tr2[[i]])  <-
    names(post_weights_tr1[[i]]) <-
    names(post_weights_tr1[[i]]) <- c(time_current + time_add)
}
rm(w_post_ex_surv_tr1, w_post_ex_surv_tr2)

# compute the maximum of posterior model averaged expected survival
post_max <- lapply(1:scenarios, function (s) {
  
  lapply(1:length(time_add), function (t) {
    
    pmax(post_ma_ex_surv_tr1[[s]][[t]], post_ma_ex_surv_tr2[[s]][[t]])
    
  })
})

# compute EVSI
evsi_mcmc <- lapply(1:scenarios, function (s) {
  
  lapply(1:length(time_add), function (t) {
    
    mean(post_max[[s]][[t]]) - max(mean(post_ma_ex_surv_tr1[[s]][[t]]), mean(post_ma_ex_surv_tr2[[s]][[t]]))
    
  })
})

# compute standard errors for the MCMC estimator
max_mean <- lapply(1:scenarios, function (s) {
  
  lapply(1:length(time_add), function (t) {
    
    if (max(mean(post_ma_ex_surv_tr1[[s]][[t]]), mean(post_ma_ex_surv_tr1[[s]][[t]])) == mean(post_ma_ex_surv_tr1[[s]][[t]])) {
      
      post_ma_ex_surv_tr1[[s]][[t]]
      
    } else {
      
      post_ma_ex_surv_tr2[[s]][[t]]
      
    }
  })
})

se_mcmc <- lapply(1:scenarios, function (s) {
  
  lapply(1:length(time_add), function (t) {
    
    sqrt(1/loops_outer*(var(post_max[[s]][[t]]) + var(max_mean[[s]][[t]]) - 2*cov(post_max[[s]][[t]], max_mean[[s]][[t]])))
    
  })
})

# store EVSI and se in a list
evsi <- matrix(unlist(evsi_mcmc), nrow = length(time_add), ncol = scenarios)
se <- matrix(unlist(se_mcmc), nrow = length(time_add), ncol = scenarios)
evsi_mcmc <- list(evsi = evsi, se = se)
rm(evsi,se)

evsi_mcmc