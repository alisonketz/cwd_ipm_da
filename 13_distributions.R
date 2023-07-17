#######################################################################
###
###   Likelihoods for each kind of data listed below
###
#######################################################################

# d_fit_hunt
# d_fit_endlive
# d_fit_sus_cens_posttest
# d_fit_sus_cens_postno
# d_fit_idead
# d_fit_sus_mort_posttest
# d_fit_sus_mort_postno
# d_fit_rec_neg_mort
# d_fit_rec_neg_cens
# d_fit_rec_pos_mort
# d_fit_rec_pos_cens
# d_fit_icap_cens
# d_fit_icap_mort

#######################################################################
###
###   User defined distribution for likelihood for
###  harvested deer, revised for multiple deer
###
###   d_fit_hunt
###
#######################################################################

dFOIhunt <- nimble::nimbleFunction( # nolint
    run = function( # nolint
                   ### argument type declarations
                   x = integer(0),
                   n_cases = double(1),
                   n_samples = integer(0), # number of samples in dataset
                   a = double(1), # age (weeks) at harvest
                   sex = double(1),
                   age2date = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double(0)) {

        # start the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {

            # intitialize scalars

            #################################################
            ### loop over ages and accumulate sums over 1:a-1
            ### have to loop separately for lam_inf
            #################################################

            if (sex[i] == 0) { # age loops for females
                lik_foi <- 0
                for (j in 1:(a[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        f_age_foi[age_lookup_f[j]] +
                        f_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_foi <- lik_foi + lam_foij 
                }
                p <- 1 - exp(-lik_foi)
                lik_temp <- dbinom(x,1,p,log=TRUE)
            } else { # age loops for males
                for (j in 1:(a[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        m_age_foi[age_lookup_f[j]] +
                        m_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_foi <- lik_foi + lam_foij 
                }
                p <- 1 - exp(-lik_foi)
                lik_temp <- dbinom(x,1,p,log=TRUE)
            } # end if(sex)
            
            #######################################
            ### accumulate the joint likelihood
            #######################################

            sumllik <- sumllik + lik_temp * n_cases[i]
        } # end the loop for individual i
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        }
    }
)

nimble::registerDistributions(list(
    dFOIhunt = list(
        BUGSdist = "dFOIhunt(n_cases,n_samples,a,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_cases = double(1)",
            "n_samples = integer(0)",
            "a = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "period_lookup_foi = double(1)",
            "space = double(1)",
            "sect = double(1)",
            "returnType = double(0)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

# for a user-defined distribution
# assign("dFOIhunt",
#     dFOIhunt,
#     envir = .GlobalEnv
# )

# starttime <- Sys.time()
# test <- dFOIhunt(
#     x = 1,
#     n_cases = d_fit_hunt$n_cases,
#     n_samples = nrow(d_fit_hunt),
#     a = d_fit_hunt$ageweeks,
#     sex = d_fit_hunt$sex,
#     age2date = d_fit_hunt$birthweek - 1,
#     f_age_foi = f_age_foi,
#     m_age_foi = m_age_foi,
#     age_lookup_f = age_lookup_f,
#     age_lookup_m = age_lookup_m,
#     period_lookup_foi = period_lookup_foi,
#     f_period_foi = f_period_foi,
#     m_period_foi = m_period_foi,
#     space = c(0, -.55),
#     sect = d_fit_hunt$ew,
#     log = TRUE
# )
# (end <- Sys.time() - starttime)
# test


########################################################################
###
###   User defined distribution for FOI from collared deer 
###   with antemortem tests (+) at capture
### 
########################################################################

dFOIicap <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(0),
        left = double(0),
        sex = double(0),
        age2date = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period = double(1),
        m_period = double(1),
        space = double(0),
        log = double()
        ) {

    logL<- 0 #intialize log-likelihood
    gam <- nimNumeric(left - 1)

    for (k in 1:(left - 1)) {
        if (sex == 0) {
            gam[k] <- space +
                      f_age[age_lookup_f[k]] +
                      f_period[period_lookup[k + age2date - 1]]
        } else {
            gam[k] <- space +
                      m_age[age_lookup_m[k]] +
                      m_period[period_lookup[k + age2date - 1]]
        }
    }
    #total probability of getting infected
    p <- 1 - exp(-sum(exp(gam[1:(left-1)])))
    logL <- dbinom(x,1,p,log=TRUE)

    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIicap = list(
        BUGSdist = 'dFOIicap(left,sex,age2date,f_age,m_age,age_lookup_f,age_lookup_m,f_period,m_period,period_lookup,space)',
        types = c('p = double(0)',
                  'left=double(0)',
                  'sex=double(0)',
                  'age2date=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup_f=double(1)',
                  'age_lookup_m=double(1)',
                  'f_period=double(1)',
                  'm_period=double(1)',
                  'period_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))

assign('dFOIicap', dFOIicap, envir = .GlobalEnv)


#############################################################
###
###   User defined distribution for FOI
###
############################################################


dFOIcollar <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(),
        left = double(0),
        right = double(0),
        sex = double(0),
        age2date = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period=double(1),
        m_period=double(1),
        space = double(0),
        log = double()
        ) {

    logL<-0 #intialize log-likelihood
    gam <-nimNumeric(right)
    for (k in left:(right-1)) {
        if(sex == 0){
            gam[k] <- space +
                      f_age[age_lookup_f[k]] +
                      f_period[period_lookup[k + age2date]]
        } else {
             gam[k] <- space +
                       m_age[age_lookup_m[k]] +
                       m_period[period_lookup[ k + age2date]]
        }
    }
    #total probability of getting infected
    p <- 1 - exp(-sum(exp(gam[left:(right - 1)])))
    logL <- dbinom(x, 1, p, log = TRUE)

    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIcollar = list(
        BUGSdist = 'dFOIcollar(left,right,age2date,sex,f_age,m_age,age_lookup_f,age_lookup_m,f_period,m_period,period_lookup,space)',
        types = c('p = double(0)',
                  'gam = double(1)',
                  'left=double(0)',
                  'right=double(0)',
                  'age2date=double(0)',
                  'sex=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup_f=double(1)',
                  'age_lookup_m=double(1)',
                  'f_period=double(1)',
                  'm_period=double(1)',
                  'period_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign('dFOIcollar', dFOIcollar, envir = .GlobalEnv)

#################################
###
###  User Defined Distribution
###  Age-Period Survival Model
###
#################################

dSurvival <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        left = double(0),
        right = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect = double(1),
        period_effect = double(1),
        nT_age = double(0),
        beta0 = double(0),
        beta_sex = double(0),
        log = double()
        ) {
    
    logL<-0 #intialize log-likelihood
   
    UCH <-nimNumeric(nT_age)

    for (k in left:(right-1)) {
        UCH[k] <- exp(beta0 + 
                      age_effect[k] +
                      period_effect[k + age2date] +
                      beta_sex * sex)
    }
    # total prob of surviving
    p <- exp(-sum(UCH[left:(right-1)]))
    logL <- dbinom(x, 1, p, log = TRUE)

    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival = list(
        BUGSdist = 'dSurvival(left,right,sex,age2date,age_effect,period_effect,nT_age,beta0,beta_sex)',
        types = c('left = double(0)',
                  'right = double(0)',
                  'sex = double(0)',
                  'age2date = double(0)',
                  'age_effect = double(1)',
                  'period_effect = double(1)',
                  'nT_age = double(0)',
                  'beta0 = double(0)',
                  'beta_sex = double(0)'
                  ),
        discrete = TRUE
    )
))
 
## for a user-defined distribution
assign('dSurvival', dSurvival, envir = .GlobalEnv)













# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   Uninfected radio-marked deer right censor:
# ###   Test neg at cap and censoring
# ###
# ###   d_fit_sus_cens_posttest
# ###   d_fit_sus_cens_postno
# ###   d_fit_endlive
# ###
# #######################################################################

# dSusCensFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0),
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    sect = double(1),
#                    space = double(1),
#                    log = double()) {

#         # starttime the loop through individuals
#         sumllik <- 0
#         for (i in 1:n_samples) {
#             # intitialize vectors
#             lam_foi <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # age loops for females
#                 for (j in e[i]:(r[i] - 1)) {
#                     # sum up foi
#                     lam_foi <- lam_foi +
#                         exp(space[sect[i]] +
#                             f_age_foi[age_lookup_f[j]] +
#                             f_period_foi[period_lookup_foi[age2date[i] + j]])
                   
#                 }
#             } else { # age loops for males
#                 for (j in e[i]:(r[i] - 1)) {
#                     # sum up foi
#                     lam_foi <- lam_foi +
#                         exp(space[sect[i]] +
#                             m_age_foi[age_lookup_m[j]] +
#                             m_period_foi[period_lookup_foi[age2date[i] + j]])
#                 }
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik - lam_foi
#         } # end loop over individuals

#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )


# nimble::registerDistributions(list(
#     dSusCensFOI = list(
#         BUGSdist = "dSusCensFOI(n_samples,e,r,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,sect,space)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "sex = double(1)",
#             "age2date = double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "period_lookup_foi = double(1)",
#             "f_period_foi = double(1)",
#             "m_period_foi = double(1)",
#             "sect = double(1)",
#             "space = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))


# ### for a user-defined distribution
# # assign("dSusCensTest", dSusCensTest, envir = .GlobalEnv)

# 		# n_samples = nrow(d_fit_sus_cens_posttest)
#         # e = d_fit_sus_cens_posttest$left_age_e
#         # r = d_fit_sus_cens_posttest$right_age_r
#         # sex = d_fit_sus_cens_posttest$sex
#         # age2date = sus_cens_posttest_age2date
#         # beta_male = beta_male
#         # beta0_sus = beta0_survival_sus
#         # age_effect_surv = age_effect_survival_test
#         # period_effect_surv = period_effect_survival_test
#         # f_age_foi = f_age_foi
#         # m_age_foi = m_age_foi
#         # age_lookup_f = age_lookup_f
#         # age_lookup_m = age_lookup_m
#         # period_lookup_foi = period_lookup_foi
#         # f_period_foi = f_period_foi
#         # m_period_foi = m_period_foi
#         # sect = d_fit_sus_cens_posttest$study_area
#         # space = c(0,-.55)

# # starttime <- Sys.time()
# # test <- dSusCensFOI(
# #         x = 1,
# # 		n_samples = nrow(d_fit_sus_cens_posttest),
# #         e = d_fit_sus_cens_posttest$left_age_e,
# #         r = d_fit_sus_cens_posttest$right_age_r,
# #         sex = d_fit_sus_cens_posttest$sex,
# #         age2date = sus_cens_posttest_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         sect = d_fit_sus_cens_posttest$study_area,
# #         space = c(0,-.55),
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test

# # starttime <- Sys.time()
# # test <- dSusCensFOI(
# #         x = 1,
# # 		  n_samples = nrow(d_fit_sus_cens_postno),
# #         e = d_fit_sus_cens_postno$left_age_e,
# #         r = d_fit_sus_cens_postno$right_age_r,
# #         sex = d_fit_sus_cens_postno$sex,
# #         age2date = sus_cens_postno_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         sect = d_fit_sus_cens_postno$study_area,
# #         space = c(0,-.55),
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# # starttime <- Sys.time()
# # test <- dSusCensFOI(
# #         x = 1,
# # 		  n_samples = nrow(d_fit_endlive),
# #         e = d_fit_endlive$left_age_e,
# #         r = d_fit_endlive$right_age_r,
# #         sex = d_fit_endlive$sex,
# #         age2date = endlive_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         sect = d_fit_endlive$study_area,
# #         space = c(0,-.55),
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test

# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected radio-marked deer mortalities:
# ###   test neg at cap and tested mort
# ###   or no test at mortality, 
# ###   assuming these individuals didn't become pos
# ###
# ###   d_fit_sus_mort_posttest
# ###   d_fit_sus_mort_postno
# ###
# ###
# #######################################################################

# dSusMortFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = double(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    space = double(1),
#                    sect = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0
#             lam_sus <- 0
#             lam_susD <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) {
#                 for (j in e[i]:(s[i] - 1)) {
#                     #  up foi hazard from 1  to j
#                     lam_foi <- lam_foi +
#                                exp(space[sect[i]] +
#                                f_age_foi[age_lookup_f[j]] +
#                                f_period_foi[period_lookup_foi[age2date[i] + j]])
#                 }
#             } else {
#                 for (j in e[i]:(s[i] - 1)) {
#                     #  up foi hazard from 1  to j
#                     lam_foi <- lam_foi +
#                                exp(space[sect[i]] +
#                                m_age_foi[age_lookup_f[j]] +
#                                m_period_foi[period_lookup_foi[age2date[i] + j]])
#                 }
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik - lam_foi
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dSusMortFOI = list(
#         BUGSdist = "dSusMortFOI(n_samples,e,r,s,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
#         types = c(
#             "value = double(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "sex = double(1)",
#             "age2date = double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "period_lookup_foi=double(1)",
#             "f_period_foi=double(1)",
#             "m_period_foi=double(1)",
#             "space = double(1)",
#             "sect = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# # ###for a user-defined distribution
# # assign("dSusMortFOI", dSusMortFOI, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dSusMortFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_sus_mort_posttest),
# #         e = d_fit_sus_mort_posttest$left_age_e,
# #         r = d_fit_sus_mort_posttest$right_age_r,
# #         s = d_fit_sus_mort_posttest$right_age_s,
# #         sex = d_fit_sus_mort_posttest$sex,
# #         age2date = sus_mort_posttest_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         space = c(0,-.55),
# #         sect = d_fit_sus_mort_posttest$study_area,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# # starttime <- Sys.time()
# # test <- dSusMortFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_sus_mort_postno),
# #         e = d_fit_sus_mort_postno$left_age_e,
# #         r = d_fit_sus_mort_postno$right_age_r,
# #         s = d_fit_sus_mort_postno$right_age_s,
# #         sex = d_fit_sus_mort_postno$sex,
# #         age2date = sus_mort_postno_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         space = c(0,-.55),
# #         sect = d_fit_sus_mort_postno$study_area,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test

# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   infected deer mortalities for radio marked deer that
# ###   enter the study as test positive at capture
# ###
# ###   d_fit_icap_cens
# ###   d_fit_icap_mort
# ###
# #######################################################################

# dIcapFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    space = double(1),
#                    sect = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # age loops for females
#                for (j in 1:(e[i] - 1)) {
#                     lam_foi <- exp(space[sect[i]] +
#                         f_age_foi[age_lookup_f[j]] +
#                         f_period_foi[period_lookup_foi[age2date[i] + j]])
#                 }
#             } else { # age loops for males
#                 for (j in 1:(e[i] - 1)) {
#                     lam_foi <- exp(space[sect[i]] +
#                         m_age_foi[age_lookup_m[j]] +
#                         m_period_foi[period_lookup_foi[age2date[i] + j]])
#                 }
               
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             p <- 1 - exp(-lam_foi)
#             sumllik <- sumllik + sum(dbinom(1,1,p,log=TRUE))
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )


# nimble::registerDistributions(list(
#     dIcapFOI = list(
#         BUGSdist = "dIcapFOI(n_samples,e,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
#         types = c(
#             "value=integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "sex = double(1)",
#             "age2date = double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "period_lookup_foi=double(1)",
#             "f_period_foi=double(1)",
#             "m_period_foi=double(1)",
#             "space = double(1)",
#             "sect = double(1)",
#             "log = double()"
#         ),
#         discrete = TRUE
#     )
# ))

# # for a user-defined distribution
# # assign("dIcapFOI", dIcapFOI, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dIcapFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_icap_cens),
# #         e = d_fit_icap_cens$left_age_e,
# #         sex = d_fit_icap_cens$sex,
# #         age2date = icap_cens_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         space = c(0,-.55),
# #         sect = d_fit_icap_cens$study_area,
# #         log = TRUE
# #         )
# # (endtime <- Sys.time() - starttime)
# # test

# # starttime <- Sys.time()
# # test <-  dIcapFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_icap_mort),
# #         e = d_fit_icap_mort$left_age_e,
# #         sex = d_fit_icap_mort$sex,
# #         age2date = icap_mort_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         space = c(0,-.55),
# #         sect = d_fit_icap_mort$study_area,
# #         log = TRUE
# #         )
# # (endtime <- Sys.time() - starttime)
# # test

# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected deer that were test neg at capture,
# ###   then test negative at recap, that are right censored,
# ###   and have been tested post censoring
# ###
# ###   d_fit_rec_neg_cens_posttest
# ###
# #######################################################################

# dRecNegCensTestFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    dn1 = double(1), # dn1, recapture, last test negative
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    sect = double(1),
#                    space = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0
#              lik_temp <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # females
#                 # force of infection infection hazard
#                 for (j in dn1[i]:(r[i] - 1)) {
#                     lam_foi <- lam_foi + exp(space[sect[i]] +
#                         f_age_foi[age_lookup_f[j]] +
#                         f_period_foi[period_lookup_foi[j + age2date[i]]])
#                 }
#             } else { # males
#                 # force of infection infection hazard
#                 for (j in dn1[i]:(r[i] - 1)) {
#                     lam_foi <- lam_foi + exp(space[sect[i]] +
#                         m_age_foi[age_lookup_m[j]] +
#                         m_period_foi[period_lookup_foi[j + age2date[i]]])
#                 }
#             }
#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik - lam_foi
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecNegCensTestFOI = list(
#         BUGSdist = "dRecNegCensTestFOI(n_samples,e,r,dn1,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "dn1 = double(1)",
#             "sex = double(1)",
#             "age2date = double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "f_period_foi=double(1)",
#             "m_period_foi=double(1)",
#             "period_lookup_foi=double(1)",
#             "space = double(1)",
#             "sect = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# ## for a user-defined distribution
# assign("dRecNegCensTestFOI", dRecNegCensTestFOI, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <- dRecNegCensTestFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_rec_neg_cens_posttest),
# #         e = d_fit_rec_neg_cens_posttest$left_age_e,
# #         r = d_fit_rec_neg_cens_posttest$right_age_r,
# #         dn1 = d_fit_rec_neg_cens_posttest$ageweek_recap,
# #         sex = d_fit_rec_neg_cens_posttest$sex,
# #         age2date = rec_neg_cens_posttest_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         sect = d_fit_rec_neg_cens_posttest$study_area,
# #         space = c(0,-.55),
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test



# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   uninfected deer that were test neg at capture,
# ###   then test negative at recap,
# ###   that die
# ###
# ###   d_fit_rec_neg_mort
# ###
# #######################################################################

# dRecNegMortFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    dn1 = double(1), #age at recapture and tested neg
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    sect = double(1),
#                    space = double(1),
#                    log = double(0)) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- 0

#             #############################################
#             # hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # age loops for females
#                 # force of infection infection hazard
#                 for (j in dn1[i]:(s[i] - 1)) {
#                     lam_foi <- lam_foi + exp(space[sect[i]] +
#                         f_age_foi[age_lookup_f[j]] +
#                         f_period_foi[period_lookup_foi[(j + age2date[i])]])
#                 }
#             } else { # males
#                # force of infection infection hazard
#                 for (j in dn1[i]:(s[i] - 1)) {
#                     lam_foi <- lam_foi + exp(space[sect[i]] +
#                         m_age_foi[age_lookup_m[j]] +
#                         m_period_foi[period_lookup_foi[(j + age2date[i])]])
#                 }
#             }

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             sumllik <- sumllik - lam_foi
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecNegMortFOI = list(
#         BUGSdist = "dRecNegMortFOI(n_samples,e,r,s,dn1,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "dn1 = double(1)",
#             "sex =  double(1)",
#             "age2date =  double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "period_lookup_foi = double(1)",
#             "f_period_foi = double(1)",
#             "m_period_foi = double(1)",
#             "space = double(1)",
#             "sect = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# # ###for a user-defined distribution
# # assign("dRecNegMortFOI", dRecNegMortFOI, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dRecNegMortFOI(
#         x = 1,
#         n_samples = nrow(d_fit_rec_neg_mort),
#         e = d_fit_rec_neg_mort$left_age_e,
#         r = d_fit_rec_neg_mort$right_age_r,
#         s = d_fit_rec_neg_mort$right_age_s,
#         dn1 = d_fit_rec_neg_mort$ageweek_recap,
#         sex = d_fit_rec_neg_mort$sex,
#         age2date = rec_neg_mort_age2date,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_rec_neg_mort$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   deer that were test neg at capture,
# ###   then test positive at recap,
# ###   than die
# ###
# ###   the following are combined into d_fit_rec_pos
# ###   d_fit_rec_pos_mort
# ###   d_fit_rec_pos_cens
# ###
# #######################################################################

# dRecPosFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    dn = double(1), # week of age recaptured and test positive
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    space = double(1),
#                    sect = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(dn[i] - 1, init = FALSE)

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # females
#                 # force of infection infection hazard
#                 lam_foi[e[i]:(dn[i] - 1)] <- exp(space[sect[i]] +
#                     f_age_foi[age_lookup_f[e[i]:(dn[i] - 1)]] +
#                     f_period_foi[period_lookup_foi[(e[i] + age2date[i]):
#                                           ((dn[i] - 1) + age2date[i])]])
#             } else { # males
#                 # force of infection infection hazard
#                 lam_foi[e[i]:(dn[i] - 1)] <- exp(space[sect[i]] +
#                     m_age_foi[age_lookup_m[e[i]:(dn[i] - 1)]] +
#                     m_period_foi[period_lookup_foi[(e[i] + age2date[i]):
#                                           ((dn[i] - 1) + age2date[i])]])
#             }

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#             p <- 1 - exp(-lam_foi[e[i]:(dn[i] - 1)])
#             sumllik <- sumllik + sum(dbinom(1,1,p,log=TRUE))
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

# nimble::registerDistributions(list(
#     dRecPosFOI = list(
#         BUGSdist = "dRecPosFOI(n_samples,e,dn,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "dn1 = double(1)",
#             "sex = double(1)",
#             "age2date = double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "period_lookup_foi=double(1)",
#             "f_period_foi=double(1)",
#             "m_period_foi=double(1)",
#             "space = double(1)",
#             "sect = double(1)",
#             "log = double(0)"
#         ),
#         discrete = TRUE
#     )
# ))

# ### Global Declaration so Nimble can access
# # assign("dRecPosMortFOI", dRecPosMortFOI, envir = .GlobalEnv)

# # starttime <- Sys.time()
# # test <-  dRecPosFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_rec_pos),
# #         e = d_fit_rec_pos$left_age_e,
# #         dn = d_fit_rec_pos$ageweek_recap,
# #         sex = d_fit_rec_pos$sex,
# #         age2date = rec_pos_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         sect = d_fit_rec_pos$study_area,
# #         space = c(0,-.55),
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test

# #######################################################################
# ###
# ###   User defined distribution for likelihood for
# ###   infected deer mortalities for radio marked deer that
# ###   enter the study as test negative at capture
# ###
# ###   d_fit_idead
# ###
# #######################################################################

# dNegCapPosMortFOI <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = integer(0),
#                    n_samples = integer(0), # number of samples in dataset
#                    e = double(1), # e, age of entry
#                    r = double(1), # r, age of last known alive
#                    s = double(1), # s, age of known mortality
#                    sex = double(1),
#                    age2date = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    space = double(1),
#                    sect = double(1),
#                    log = double()) {
#         sumllik <- 0 # intialize log-likelihood
#         for (i in 1:n_samples) {
#             lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
#             lik_temp <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # females
#                 # force of infection infection hazard
#                 lam_foi[e[i]:(s[i] - 1)] <- exp(space[sect[i]] +
#                     f_age_foi[age_lookup_f[e[i]:(s[i] - 1)]] +
#                     f_period_foi[period_lookup_foi[(e[i] + age2date[i]):
#                                                (s[i]- 1 + age2date[i])]])
#             } else { # males
#                 # force of infection infection hazard
#                 lam_foi[e[i]:(s[i] - 1)] <- exp(space[sect[i]] +
#                     m_age_foi[age_lookup_m[e[i]:(s[i] - 1)]] +
#                     m_period_foi[period_lookup_foi[(e[i] + age2date[i]):
#                                            (s[i] - 1 + age2date[i])]])
#             }

#             #######################################
#             ### calculating the joint likelihood
#             #######################################
#         p <- 1 - exp(-lam_foi[e[i]:(s[i] - 1)])
#         sumllik <- sumllik + sum(dbinom(1,1,p,log=TRUE))
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )


# nimble::registerDistributions(list(
#     dNegCapPosMortFOI = list(
#         BUGSdist = "dNegCapPosMortFOI(n_samples,e,r,s,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
#         types = c(
#             "value = integer(0)",
#             "n_samples = integer(0)",
#             "e = double(1)",
#             "r = double(1)",
#             "s = double(1)",
#             "sex = double(1)",
#             "age2date = double(1)",
#             "f_age_foi = double(1)",
#             "m_age_foi = double(1)",
#             "age_lookup_f = double(1)",
#             "age_lookup_m = double(1)",
#             "period_lookup_foi = double(1)",
#             "f_period_foi = double(1)",
#             "m_period_foi = double(1)",
#             "space = double(1)",
#             "sect = double(1)",
#             "log = double()"
#         ),
#         discrete = TRUE
#     )
# ))

# ### for a user-defined distribution
# # assign("dNegCapPosMortFOI", dNegCapPosMortFOI, envir = .GlobalEnv)


# # starttime <- Sys.time()
# # test <-  dNegCapPosMortFOI(
# #         x = 1,
# #         n_samples = nrow(d_fit_idead),
# #         e = d_fit_idead$left_age_e,
# #         r = d_fit_idead$right_age_r,
# #         s = d_fit_idead$right_age_s,
# #         sex = d_fit_idead$sex,
# #         age2date = idead_age2date,
# #         f_age_foi = f_age_foi,
# #         m_age_foi = m_age_foi,
# #         age_lookup_f = age_lookup_f,
# #         age_lookup_m = age_lookup_m,
# #         period_lookup_foi = period_lookup_foi,
# #         f_period_foi = f_period_foi,
# #         m_period_foi = m_period_foi,
# #         space = c(0,-.55),
# #         sect = d_fit_idead$study_area,
# #         log = TRUE
# #         )
# # (end<- Sys.time()-starttime)
# # test


# #################################
# ###
# ###  User Defined Distribution
# ###  Age-Period Survival Model
# ###
# #################################

# dSurvival <- nimble::nimbleFunction(
#     run = function(
#         ### argument type declarations
#         x = double(),
#         left = double(0),
#         right = double(0),
#         age_effect = double(1),
#         period_effect = double(1),
#         age2date = double(0),
#         nT_age = double(0),
#         beta0 = double(0),
#         beta_sex = double(0),
#         sex = double(0),
#         log = double()
#         ) {
    
#     logL<-0 #intialize log-likelihood
   
#     UCH <-nimNumeric(nT_age)

#     for (k in left:(right-1)) {
#         UCH[k] <- exp(beta0 + 
#                       age_effect[k] +
#                       period_effect[k + age2date] +
#                       beta_sex * sex)
#     }
#     # total prob of surviving
#     p <- exp(-sum(UCH[left:(right-1)]))
    
#     logL <- log(p)

#     returnType(double())
#     if(log) return(logL) else return(exp(logL))    ## return log-likelihood
#   })

# nimble::registerDistributions(list(
#     dSurvival = list(
#         BUGSdist = 'dSurvival(left,right,age_effect,period_effect,age2date,nT_age,beta0,beta_sex,sex)',
#         types = c('left = double(0)',
#                   'right = double(0)',
#                   'age_effect = double(1)',
#                   'period_effect = double(1)',
#                   'age2date = double(0)',
#                   'nT_age = double(0)',
#                   'beta0 = double(0)',
#                   'beta_sex = double(0)',
#                   'sex = double(0)'
#                   ),
#         discrete = TRUE
#     )
# ))
 
# ## for a user-defined distribution
# assign('dSurvival', dSurvival, envir = .GlobalEnv)
