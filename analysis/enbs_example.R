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


#####################################################################################
# Analysis settings
#####################################################################################

# number of patients enrolled in each trial arm
npat_tr1 <- 200
npat_tr2 <- 200

# current follow-up time
time_current <- 12

# additional follow-up times
time_add <- c(5, 12,  24,  36,  48,  72,  96, 120) 

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
dists <- c("weibull", "gamma", "lnorm", "llogis") 

# time horizon in months for computing expected survival
t_horizon <- 240

# number of datasets to be generated in the outer loop of the EVSI procedure
loops_outer <- 6e3

# number of prior draws for computing EVPI (not used in EVSI calculations)
n_draws <- 100e3


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

col <- gg_color_hue(3) # set color scheme
font_size <- 14

# plot increasing hazard dataset and extract legend
plot1 <- plot_km_fun(ipd[[1]], title = "", y_lab = "Survival probability", font_size = font_size, risktab = TRUE)
mylegend <- g_legend(plot1$plot) # extract legend
plot1$plot <- plot1$plot + theme(legend.position = "none") # remove legend in plot

# plot decreasing hazard dataset
plot2 <- plot_km_fun(ipd[[2]], title = "", y_lab = "", font_size = font_size, risktab = TRUE)
plot2$plot <- plot2$plot + theme(legend.position = "none") # remove legend in plot

# generate ggplot2 grobs
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
# Simple ENBS example
#####################################################################################

n_pat <- 5 # number of patients treated each month

trial_cost <- 5 # trial costs per month

###########################
# Increasing hazard dataset
###########################

# prior incremental net benefit
incr_nb1 <- prior_ma_ex_surv_tr1$`Increasing hazard` - prior_ma_ex_surv_tr2$`Increasing hazard`

# create a dataframe with marginal EVSI and cost functions
marg_df1 <- marg_fun(evsi_gam$evsi[,1], time_add, n_pat, trial_cost, incr_nb1)

# find intersections between marginal EVSI and cost functions
x_1 <- intersect_fun(x = subset(marg_df1$times, marg_df1$group ==1), 
                     y = subset(marg_df1$NB, marg_df1$group ==1),
                     x2 = subset(marg_df1$times, marg_df1$group ==2),
                     y2 = subset(marg_df1$NB, marg_df1$group ==2),
                     step = 0.01)

x_2 <- intersect_fun(x = subset(marg_df1$times, marg_df1$group ==1), 
                     y = subset(marg_df1$NB, marg_df1$group ==1),
                     x2 = subset(marg_df1$times, marg_df1$group ==3),
                     y2 = subset(marg_df1$NB, marg_df1$group ==3),
                     step = 0.01)

y_1 <- subset(marg_df1$NB, marg_df1$group ==2)[1]
y_2 <- subset(marg_df1$NB, marg_df1$group ==3)[1]

anno1 <- bquote(~ italic("t"^"*") == .(round(x_1,0)))
anno2 <- bquote(~ italic("t"^"*") == .(round(x_2,0)))

plot1 <- plot_enbs_fun(marg_df1, "Net benefit", "top", x_1, x_2, y_1, y_2, anno1, anno2) +
  scale_colour_manual(
    "",
    values = col,
    breaks = c("1", "2", "3"),
    labels = c(expression(paste("MB"[EVSI], sep = '')),
               expression(paste("MC"[AWR], sep = '')),
               expression(paste("MC"[OIR], sep = '')))
  )

# extract legend
mylegend <- g_legend(plot1) 

# remove legend from plot
plot1 <- plot1 + theme(legend.position = "none") 


###########################
# Decreasing hazard dataset
###########################

# prior incremental net benefit
incr_nb2 <- prior_ma_ex_surv_tr1$`Decreasing hazard` - prior_ma_ex_surv_tr2$`Decreasing hazard`

# create a dataframe with marginal EVSI and cost functions
marg_df2 <- marg_fun(evsi_gam$evsi[,2], time_add, n_pat, trial_cost, incr_nb2)

# find intersections between marginal EVSI and cost functions
x_1 <- intersect_fun(x = subset(marg_df2$times, marg_df2$group ==1), 
                     y = subset(marg_df2$NB, marg_df2$group ==1),
                     x2 = subset(marg_df2$times, marg_df2$group ==2),
                     y2 = subset(marg_df2$NB, marg_df2$group ==2),
                     step = 0.01)

x_2 <- intersect_fun(x = subset(marg_df2$times, marg_df2$group ==1), 
                     y = subset(marg_df2$NB, marg_df2$group ==1),
                     x2 = subset(marg_df2$times, marg_df2$group ==3),
                     y2 = subset(marg_df2$NB, marg_df2$group ==3),
                     step = 0.01)

y_1 <- subset(marg_df2$NB, marg_df2$group ==2)[1]
y_2 <- subset(marg_df2$NB, marg_df2$group ==3)[1]

anno1 <- bquote(~ italic("t"^"*") == .(round(x_1,0)))
anno2 <- bquote(~ italic("t"^"*") == .(round(x_2,0)))

plot2 <- plot_enbs_fun(marg_df2, "", "none", x_1, x_2, y_1, y_2, anno1, anno2)

###########################
# arrange plots in a grid
###########################

gA <- ggplotGrob(plot1)
gB <- ggplotGrob(plot2)

m_layout <- rbind(c(1,1), do.call("rbind", replicate(8, c(2,3), simplify = FALSE)))

grid.arrange(arrangeGrob(mylegend), gA, gB, layout_matrix = m_layout)