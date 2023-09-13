#######################################
### Data for Model Fitting
#######################################

nimData <- list(Z_period = Z_period,
    Z_age = Z_age,
    num_period = num_period,
    adj_period = adj_period,
    weights_period = weights_period,
    age_lookup_f = age_lookup_f,
    age_lookup_m = age_lookup_m,
    # period_effect_survival = rep(NA,nT_period_overall),
    y_hunt_foi = 1,
    y_hunt_foi_teststatus = d_fit_hunt$teststatus,
    n_cases_foi_hunt = d_fit_hunt$n_cases,
    a_fit_hunt = d_fit_hunt$ageweeks,
    sex_fit_hunt = d_fit_hunt$sex,
    age2date_fit_hunt = d_fit_hunt$birthweek,
    y_sus_foi_teststatus = d_fit_sus_foi$teststatus,
    left_sus_foi = d_fit_sus_foi$left,
    right_sus_foi = d_fit_sus_foi$right,
    sex_sus_foi = d_fit_sus_foi$sex,
    age2date_sus_foi = d_fit_sus_foi$age2date,
    y_icap_foi_teststatus = d_fit_icap_foi$teststatus,
    left_icap_foi = d_fit_icap_foi$left,
    sex_icap_foi = d_fit_icap_foi$sex,
    age2date_icap_foi = d_fit_icap_foi$age2date,
    y_rec_neg_cens_postno_foi_teststatus = d_fit_rec_neg_cens_postno_foi$teststatus,
    left_rec_neg_cens_postno_foi = d_fit_rec_neg_cens_postno_foi$left,
    right_rec_neg_cens_postno_foi = d_fit_rec_neg_cens_postno_foi$right,
    sex_rec_neg_cens_postno_foi = d_fit_rec_neg_cens_postno_foi$sex,
    age2date_rec_neg_cens_postno_foi = d_fit_rec_neg_cens_postno_foi$age2date,
    y_recap_foi_teststatus = d_fit_recap_foi$teststatus,
    left_recap_foi = d_fit_recap_foi$left,
    right_recap_foi = d_fit_recap_foi$right,
    sex_recap_foi = d_fit_recap_foi$sex,
    age2date_recap_foi = d_fit_recap_foi$age2date,
    y_idead_foi_teststatus = d_fit_idead_foi$teststatus,
    left_idead_foi = d_fit_idead_foi$left,
    right_idead_foi = d_fit_idead_foi$right,
    sex_idead_foi = d_fit_idead_foi$sex,
    age2date_idead_foi = d_fit_idead_foi$age2date,
    #################################################################
    #### survival data
    #################################################################
    y_sus_surv = d_fit_sus$surv_censor,
    left_sus = d_fit_sus$left,
    right_sus = d_fit_sus$right,
    sex_sus = d_fit_sus$sex,
    age2date_sus = d_fit_sus$age2date,
    y_icap_surv = d_fit_icap$surv_censor,
    left_icap = d_fit_icap$left,
    right_icap = d_fit_icap$right,
    sex_icap = d_fit_icap$sex,
    age2date_icap = d_fit_icap$age2date,
    y_idead = d_fit_idead$surv_censor,
    idead_left_age_e = d_fit_idead$left_age_e,
    idead_right_age_r = d_fit_idead$right_age_r,
    idead_right_age_s = d_fit_idead$right_age_s,
    idead_left_period_e = d_fit_idead$left_period_e,
    idead_right_period_s = d_fit_idead$right_period_s,
    idead_sex = d_fit_idead$sex,
    idead_age2date =  d_fit_idead$age2date,
    y_rec_pos_cens = 1,
    rec_pos_cens_left_age_e = d_fit_rec_pos_cens$left_age_e,
    rec_pos_cens_right_age_r = d_fit_rec_pos_cens$right_age_r,
    rec_pos_cens_ageweek_recap = d_fit_rec_pos_cens$ageweek_recap,
    rec_pos_cens_left_period_e = d_fit_rec_pos_cens$left_period_e,
    rec_pos_cens_periodweek_recap = d_fit_rec_pos_cens$periodweek_recap,
    rec_pos_cens_sex = d_fit_rec_pos_cens$sex,
    rec_pos_cens_age2date = d_fit_rec_pos_cens$age2date,
    y_rec_pos_mort = rep(1,n_fit_rec_pos_mort),
    rec_pos_mort_left_age_e = d_fit_rec_pos_mort$left_age_e,
    rec_pos_mort_right_age_r = d_fit_rec_pos_mort$right_age_r,
    rec_pos_mort_right_age_s = d_fit_rec_pos_mort$right_age_s,
    rec_pos_mort_ageweek_recap = d_fit_rec_pos_mort$ageweek_recap,
    rec_pos_mort_left_period_e = d_fit_rec_pos_mort$left_period_e,
    rec_pos_mort_periodweek_recap = d_fit_rec_pos_mort$periodweek_recap,
    rec_pos_mort_sex = d_fit_rec_pos_mort$sex,
    rec_pos_mort_age2date = d_fit_rec_pos_mort$age2date,
    y_sus_draw = rep(1,n_fit_sus_draw),
    sus_draw_left_age_e = d_fit_sus_draw$left_age_e,
    sus_draw_right_age_r = d_fit_sus_draw$right_age_r,
    sus_draw_left_period_e = d_fit_sus_draw$left_period_e ,
    sus_draw_right_period_r = d_fit_sus_draw$right_period_r,
    sus_draw_sex = d_fit_sus_draw$sex,
    sus_draw_age2date = d_fit_sus_draw$age2date,
    y_sus_mort_postno = rep(1,n_fit_sus_mort_postno),
    sus_mort_postno_left_age_e = d_fit_sus_mort_postno$left_age_e,
    sus_mort_postno_right_age_r = d_fit_sus_mort_postno$right_age_r,
    sus_mort_postno_right_age_s = d_fit_sus_mort_postno$right_age_s,
    sus_mort_postno_left_period_e = d_fit_sus_mort_postno$left_period_e ,
    sus_mort_postno_right_period_s = d_fit_sus_mort_postno$right_period_s,
    sus_mort_postno_sex = d_fit_sus_mort_postno$sex,
    sus_mort_postno_age2date = d_fit_sus_mort_postno$age2date,
    y_rec_neg_cens_postno = rep(1,n_fit_rec_neg_cens_postno),
    rec_neg_cens_postno_left_age_e = d_fit_rec_neg_cens_postno$left_age_e,
    rec_neg_cens_postno_right_age_r = d_fit_rec_neg_cens_postno$right_age_r,
    rec_neg_cens_postno_ageweek_recap = d_fit_rec_neg_cens_postno$ageweek_recap,
    rec_neg_cens_postno_left_period_e = d_fit_rec_neg_cens_postno$left_period_e,
    rec_neg_cens_postno_right_period_r = d_fit_rec_neg_cens_postno$right_period_r,
    rec_neg_cens_postno_periodweek_recap = d_fit_rec_neg_cens_postno$periodweek_recap,
    rec_neg_cens_postno_sex = d_fit_rec_neg_cens_postno$sex,
    rec_neg_cens_postno_age2date = d_fit_rec_neg_cens_postno$age2date,
    mort_hh = d_fit_hh$mort_h,
    sex_cause = d_fit_hh$sex,
    Z_cause_gun = Z_cause_gun,
    Z_cause_ng = Z_cause_ng,
    Cage_less = Cage_less,
    Cage_ant = Cage_ant,
    # O = Ototal,
    lobs = log(Ototal),
    f_logpop_sus = f_logpop_sus,
    f_logpop_inf = f_logpop_inf,
    m_logpop_sus = m_logpop_sus,
    m_logpop_inf = m_logpop_inf,
    obs_ct_fd_alpha = obs_ct_fd_alpha,
    obs_ct_fd_beta = obs_ct_fd_beta,
    Nfawn = fawndoe_df$overall_fawn,
    Ndoe = fawndoe_df$overall_doe#,
    # f_period_foi = c(rep(0,8),rep(NA,n_year - 8)),
    # m_period_foi = c(rep(0,8),rep(NA,n_year - 8)),
    # f_period_foi_temp = c(rep(0,8),rep(NA,n_year - 8)),
    # m_period_foi_temp = c(rep(0,8),rep(NA,n_year - 8))
    )

#######################################
### Constants for MCMC
#######################################

nimConsts <- list(n_year = n_year,
    n_year_precollar = n_year_precollar,
    n_study_area = n_study_area,
    n_sex = n_sex,
    n_agef = n_agef,
    n_agem = n_agem,
    n_ageclassf = n_ageclassf,
    n_ageclassm = n_ageclassm,
    n_sex = n_sex,
    # cal = d_fit_season$ng_end - d_fit_season$yr_start,
    sizeCage_f = sizeCage_f,
    sizeCage_m = sizeCage_m,
    report_hyp_all = report_hyp_all,
    report_hyp_y = report_hyp_y,
    nT_period_overall = nT_period_overall,
    nT_period_precollar = nT_period_precollar,
    nT_period_collar = nT_period_collar,
    nT_period_overall_hunt = nT_period_overall_hunt,
    nT_age_surv = nT_age_surv,
    nT_age_surv_aah_f = nT_age_surv_aah_f,
    nT_age_surv_aah_m = nT_age_surv_aah_m,
    nT_age_short_f = nT_age_short_f,
    nT_age_short_m = nT_age_short_m,
    n_year_fec_early = n_year_fec_early,
    n_age_lookup_f = length(age_lookup_f),
    n_age_lookup_m = length(age_lookup_m),
    nknots_age = nknots_age,
    nknots_period = nknots_period,
    n_adj_period = n_adj_period,
    period_lookup_foi = period_lookup_foi,
    period_lookup_foi_hunt = period_lookup_foi_hunt,
    ng_start = d_fit_season$ng_start,
    gun_start = d_fit_season$gun_start,
    gun_end = d_fit_season$gun_end,
    ng_end = d_fit_season$ng_end,
    yr_start = d_fit_season$yr_start,
    yr_end = d_fit_season$yr_end,
    n_fit_hunt = n_fit_hunt,
    n_fit_sus_foi = n_fit_sus_foi,
    n_fit_icap_foi = n_fit_icap_foi,
    n_fit_rec_neg_cens_postno_foi = n_fit_rec_neg_cens_postno_foi,
    n_fit_recap_foi = n_fit_recap_foi,
    n_fit_idead_foi = n_fit_idead_foi,
    n_fit_idead= n_fit_idead,
    n_fit_sus = n_fit_sus,
    n_fit_icap = n_fit_icap,
    n_fit_rec_pos_mort = n_fit_rec_pos_mort,
    n_fit_sus_draw = n_fit_sus_draw,
    n_fit_sus_mort_postno = n_fit_sus_mort_postno,
    n_fit_rec_neg_cens_postno = n_fit_rec_neg_cens_postno,
    sect_hunt = d_fit_hunt$ew,
    sect_sus_foi = d_fit_sus_foi$study_area,
    sect_icap_foi = d_fit_icap_foi$study_area,
    sect_rec_neg_cens_postno_foi = d_fit_rec_neg_cens_postno_foi$study_area,
    sect_recap_foi = d_fit_recap_foi$study_area,
    sect_idead_foi = d_fit_idead_foi$study_area,
    sect_idead = d_fit_idead$study_area,
    sect_rec_pos_cens = d_fit_rec_pos_cens$study_area,
    sect_rec_pos_mort = d_fit_rec_pos_mort$study_area,
    sect_sus_draw = d_fit_sus_draw$study_area,
    sect_sus_mort_postno = d_fit_sus_mort_postno$study_area,
    sect_rec_neg_cens_postno = d_fit_rec_neg_cens_postno$study_area,
    records_cause = records_cause,
    interval_cause = d_fit_hh$right_period_s - 1,
    period_annual_survival = rep(0,n_year_precollar + 1),
    # f_period_foi = rep(0, n_year),
    # m_period_foi = rep(0, n_year)
    # indx_mat_pe_surv = indx_mat_pe_surv,
    # intvl_step_yr = intvl_step_yr_weekly#,
    # num_foi_cal = num_foi_cal
    Z_collar_gun = Z_collar_gun,
    Z_collar_ng = Z_collar_ng,
    nconst = 1 / sqrt(2 * pi)
    )


#######################################
### Initial Values for MCMC
#######################################

initsFun <- function()list(beta_male = rnorm(1, -.5, .01),
    beta0_sus_temp = rnorm(1, -5.5, 0.0001),
    # beta0_survival_sus = rnorm(1, -8, 0.0001),
    sus_mix = 1,
    beta0_inf_temp = rnorm(1, -4, 0.0001),
    # beta0_survival_inf = rnorm(1, -7, 0.0001),
    inf_mix = 1,
    ln_b_age_survival = rnorm(nknots_age) * 10^-4,
    # b_period_survival = rnorm(nknots_period) * 10^-4,
    # tau_period_survival = runif(1, .1, 1),
    beta_harvest_gun = rnorm(1, 0, sd = 1),
    beta_harvest_ng = rnorm(1, 0, sd = 1),
    mix_survival = 1,
    ln_sk_period = runif(1, .1, 1),
    sda_period = runif(1, .1, 1),
    alpha_period = rnorm(nknots_period, 0, 1),
    tau_age_survival = runif(1, .1, .4),
    tau_age_foi_male = runif(1, 1.5, 1.7),
    tau_age_foi_female = runif(1, 2.7, 4.2),
    tau_period_foi_male = runif(1, 4.2, 6.8),
    tau_period_foi_female = runif(1, 2.27, 3.44),
    m_period_foi = seq(-.25, .25, length = n_year),
    f_period_foi = seq(-.25, .25, length = n_year),
    m_age_foi = c(rnorm(1, -6, sd = .1),
                  rnorm(1, -5.5, sd = .1),
                  rnorm(1, -5, sd = .1),
                  rnorm(1, -5.5, sd = .1),
                  rnorm(1, -6, sd = .1),
                  rnorm(1, -7.2, sd = .1)) - 1.5,
    f_age_foi = c(rnorm(1, -6, sd = .1),
                  rnorm(1, -5.5, sd = .1),
                  rnorm(1, -6, sd = .1),
                  rnorm(1, -6.5, sd = .1),
                  rnorm(1, -6.8, sd = .1),
                  rnorm(1, -7.2, sd = .1),
                  rnorm(1, -8, sd = .1)) - 1.5,
    # tau_period_precollar = rgamma(1, 1, 1),
    # period_annual_survival = rnorm(n_year_precollar + 1, .1),
    beta0_cause = rnorm(1, -2.8, .1),
    beta_cause_male = rnorm(1, 0, .1),
    beta_cause_gun = rnorm(1, 1.5, .1),
    beta_cause_ng = rnorm(1, 3, .1),
    # tau_obs = matrix(runif(4, 1, 3), 2, 2),
    # tau_obs = runif(1, .01, .5),
    # tau_pop = runif(2, .05, .25),
    # report_overall = report_overall_init,
    # report = report_init,
    # fec_epsilon = fec_eps_init,
    # mu_fec = rnorm(1, mu_fec_init, .01),
    # fec_prec_eps = runif(1, 5, 10),
    space_temp = rnorm(1, -.55, .01),
    space_mix = 1
    )
nimInits <- initsFun()

########################################################
### Build and compile model in R
########################################################

# start_Rmodel <- Sys.time()
# Rmodel <- nimbleModel(code = modelcode,
#                       constants = nimConsts,
#                       data = nimData,
#                       inits = initsFun(),
#                       calculate = FALSE,
#                       check = FALSE
#                       )
# # end_Rmodel <- Sys.time() - start_Rmodel
# Rmodel$initializeInfo()

# Cnim <- compileNimble(Rmodel)

# for(i in 1:10){beepr::beep(1)}

#######################################
### Parameters to trace in MCMC
#######################################

parameters <- c(
              "beta_male",
              "tau_age_foi_male",
              "tau_age_foi_female",
            #   "tau1_age_foi_male",
            #   "tau1_age_foi_female",
              "m_age_foi",
              "f_age_foi",
              "m_age_foi_mu",
              "f_age_foi_mu",
              "tau_period_foi_male",
              "tau_period_foi_female",
              "f_period_foi",
              "m_period_foi",
              "space",
              "space_temp",
              "space_mix",
            #   "beta0_sus_temp",
            #   "sus_mix",
              "beta0_survival_sus",
              "tau_age_survival",
              "age_effect_survival",
              "ln_b_age_survival",
              # "b_period_survival",
              # "tau_period_survival",
              "mix_survival",
              "ln_sk_period",
              "sdk_period",
              "tauk_period",
              "stauk_period",
              "sda_period",
              "taua_period",
              "alpha_period",
              "alphau_period",
              "ratioinf_period",
              "beta_harvest_gun",
              "beta_harvest_ng",
              # "period_effect_surv",
              # "tau_period_precollar",
              "period_effect_survival",
            #   "beta0_inf_temp",
            #   "inf_mix",
              "beta0_survival_inf",
              "beta0_cause",
              "beta_cause_gun",
              "beta_cause_ng",
              "beta_cause_male",
              "p_nogun_f",
              "p_gun_f",
              "p_nogun_m",
              "p_gun_m",
              # "report",
              # "fec",
              # "mu_fec",
              # "fec_prec_eps",
              # "fec_epsilon",
              "psi",
              "psi_hat",
              "sn_inf",
              "sn_sus",
              "sh_inf",
              "sh_sus"#,
              # "mu_obs",
              # "tau_obs",
              # "tau_pop"
               )

# confMCMC <- configureMCMC(Rmodel,
#                          monitors = parameters,
#                          thin = 1,
#                          # enableWAIC = TRUE,
#                          useConjugacy = FALSE)
# nimMCMC <- buildMCMC(confMCMC)
# CnimMCMC <- compileNimble(nimMCMC,
#                          project = Rmodel)
# for(i in 1:10){beepr::beep(1)}

# set.seed(7654321)
# starttime <- Sys.time()
# mcmcout <- runMCMC(CnimMCMC,
#                   niter = 1000,
#                   nburnin = 0,
#                   nchains = 1,
#                   inits = initsFun,
#                   samplesAsCodaMCMC = TRUE,
#                   summary = TRUE
#                   )
# runtime <- difftime(Sys.time(),
#                     starttime,
#                     units = "min")
# for (i in 1:10) {beepr::beep(1)}


# end_Rmodel
# endtime_rmodel_compile
# endtime_mcmc
# runtime

# sink("runtime_allsteps.txt")
# cat("Rmodel:\n")
# end_Rmodel
# cat("\nCompile Rmodel:\n")
# endtime_rmodel_compile
# cat("\nCompile MCMC:\n")
# endtime_mcmc
# cat("\nRun MCMC 100 iter: ",runtime)
# sink()

#############################################################
###
### Running in parallel, restartable
###
#############################################################

ni  <- 100
nb <- .5
bin <- ni * nb
nt <- 1
nc <- 3
starttime <- Sys.time()
cl <- makeCluster(nc, timeout = 5184000)

clusterExport(cl, c("modelcode",
                    "initsFun",
                    "nimData",
                    "nimConsts",
                    "parameters",
                    "ni",
                    "bin",
                    "nt"
                    ))
for (j in seq_along(cl)) {
  set.seed(j + 1000)
  init <- initsFun()
  clusterExport(cl[j], "init")
}
for (i in 1:10) {beepr::beep(1)}

starttime <- Sys.time()
mcmcout1 <-  mcmc.list(clusterEvalQ(cl, {
  library(nimble)
  library(coda)

  source("13_distributions.R")
  source("14_calculations.R")

  ##############################################################
  ###
  ### Execute MCMC
  ###
  ##############################################################

  Rmodel <- nimbleModel(code = modelcode,
                        name = "modelcode",
                        constants = nimConsts,
                        data = nimData,
                        inits = init)
  Cnim <- compileNimble(Rmodel)
  confMCMC <- configureMCMC(Rmodel,
                            monitors = parameters,
                            thin = nt,
                            useConjugacy = FALSE)
  # confMCMC$removeSamplers(c("beta0_survival_sus", "beta0_survival_inf")) 
  # confMCMC$addSampler(target = "beta0_survival_sus", type = "slice")
  # confMCMC$addSampler(target = "beta0_survival_inf",  type = "slice")
  nimMCMC <- buildMCMC(confMCMC)
  CnimMCMC <- compileNimble(nimMCMC,
                            project = Rmodel)

  CnimMCMC$run(ni, reset = FALSE)

  return(as.mcmc(as.matrix(CnimMCMC$mvSamples)))
}))
(runtime1 <- difftime(Sys.time(), starttime, units = "min"))



for(chn in 1:nc) { # nc must be > 1
  ind_keep <- c()
  for(p in 1:length(parameters)) {
    ind_keep <- c(ind_keep,
        which(str_detect(dimnames(mcmcout1[[chn]])[[2]],
                         parameters[p]))) %>%
        unique()
  }
  mcmcout1[[chn]] <- mcmcout1[[chn]][,ind_keep]
}

## Check convergence ##
out1 <- mcmcout1
ni_saved <- nrow(out1[[1]])
for(chn in 1:nc) { # nc must be > 1
  
  if(nb < 1) {
    nb_real <- (round(ni_saved * nb)+1)
  } else {
    nb_real <- (round(nb/nt)+1)
  }
  out1[[chn]] <- out1[[chn]][nb_real:ni_saved,]
}
out_mcmc1 <- coda::as.mcmc.list(lapply(out1, coda::as.mcmc))
mod <- mcmcOutput::mcmcOutput(out_mcmc1)
sumTab <- summary(mod,
                  MCEpc = FALSE,
                  Rhat = TRUE,
                  n.eff = TRUE, 
                  f = TRUE,
                  overlap0 = TRUE,
                  verbose = FALSE)

sumTab <- sumTab %>%
      as_tibble() %>%
      mutate(Parameter = row.names(sumTab)) %>%
      select(Parameter, mean:f)

sumTab

mcmcout2 <- clusterEvalQ(cl, {
  CnimMCMC$run(ni, reset = FALSE, resetMV = TRUE) # Resume sampling.
  return(as.mcmc(as.matrix(CnimMCMC$mvSamples)))
  gc(verbose = F)
})

for(chn in 1:nc) { # nc must be > 1
  ind.keep <- c()
  for(p in 1:length(parameters)) ind.keep <-
      c(ind.keep, which(str_detect(dimnames(mcmcout2[[chn]])[[2]], parameters[p]))) %>% unique()
  mcmcout2[[chn]] <- mcmcout2[[chn]][,ind.keep]
}
mod.nam = "mod"
n.runs <- 2
R.utils::saveObject(mcmcout2, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.



# par.ignore.Rht=c()
# if(length(par.ignore.Rht) == 0) {
#       mxRht <- sumTab %>% pull(Rhat) %>% max(na.rm = T)
#       mn.neff <- sumTab %>% pull(n.eff) %>% min(na.rm = T)
# } else {
#     ind.ignore <- c()
#     for(p in 1:length(par.ignore.Rht)) ind.ignore <-
#         c(ind.ignore, which(str_detect(sumTab$Parameter, par.ignore.Rht[p]))) %>%
#         unique()
#     if(length(ind.ignore) > 0) {
#         mxRht <- sumTab %>% slice(-ind.ignore) %>% pull(Rhat) %>% max(na.rm = T)
#         mn.neff <- sumTab %>% slice(-ind.ignore) %>% pull(n.eff) %>% min(na.rm = T)
#     } else {
#         mxRht <- sumTab %>% pull(Rhat) %>% max(na.rm = T)
#         mn.neff <- sumTab %>% pull(n.eff) %>% min(na.rm = T)
#     }
# }

















# stopCluster(cl)

# save(mcmcmcmcout1, file = "mcmcmcmcout1.Rdata")
# save(runtime1, file = "runtime1.Rdata")
# save(endtime_rmodel_compile, file = "endtime_rmodel_compile.Rdata")
# save(endtime_mcmc, file = "endtime_mcmc.Rdata")

###
### save model run
###

# save(runtime,file="results/runtime.Rdata")
# save(mcmcout,file="results/mcmcout.Rdata")