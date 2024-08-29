# CHAPTER 3: SCENARIO ANALYSIS #################################################
################################################################################
### 1. Load packages ####
library(pacman)
p_load(deSolve, tidyverse, gridExtra, readxl, ggplot2)

### 2. Input definitions ####
# 2.1. define the number of days to run the model ####
t_start<- 0
n <- 10
t_end <- n*360
step <- 30

times <- seq(t_start, t_end, step)

# 2.2. Parameters ####
falciparum <- read_excel("Parameters and initial condition - falciparum.xlsx", sheet = "default") # loading excel data
parameters_value.vector <- c(falciparum$value)  # taking the value of parameters from data frame into a vector
parameters_name.vector <- c(falciparum$name)    # taking the name of parameters from data frame into a vector

names(parameters_value.vector) <- parameters_name.vector  # giving the parameters value with name

parms <- parameters_value.vector

# 2.3. ITN data ####
itndata <- read_excel("itndata.xlsx")

# 2.4. Initial conditions ####
initSf <- S0.fit  # initial susceptible population
initEf <- 1000  # initial number of exposed compartment
initAf <- 0     # initial number of asymptomatic malaria
initCf <- 1000  # initial number of malaria cases with clinical symptoms
initCtf <- 0    # initial number of symptomatic malaria who are destined to be treated
initTf <- 0     # initial number of treated cases
initRf <- 3000  # initial number of recovered individuals
initSev <- 0    # initial number of severe malria cases
initH <- 0      # initial number of hospitalised patients
initCinc <- 0   # initial value for incidence counter
initCtr <- 0    # initial value for treatment counter
initCtrt_total <- 0 # initial value for total treatment counter (outpatients + hospitalised)
initCtest <- 0  # initial value for test counter

start <- c(Sf = S0.fit,          # initial susceptible cases
           Ef = 1000,            # initial exposed cases
           Af = 0,               # initial asymptomatic cases
           Cf = 1000,            # initial uncomplicated cases who are not going to be treated
           Ctf = 0,              # initial uncomplicated cases who are destined to be treated
           Tf = 0,               # initial treated patients
           Sev = 0,              # initial severe patients
           H = 0,                # initial hospitalised patients
           Rf = 3000,            # initial recovered inidividuals
           Rpf = 0,              # initial
           Cinc = 0,             # initial cumulative incidence 
           Ctrt = 0,             # initial number of uncomplicated malaria patients who receive treatment (outpatient)
           Chosp = 0,            # initial number of complicated malaria pastients who receive treatment (hospitalised)
           Ctrt_total = 0,       # initial number of those receiving treatment by all means
           Ctest = 0,            # initial number of tests
           ITN = 0,              # assigning the initial number of ITN
           Death = 0,            # initial death
           Sfm = 0,              # 
           Tfm = 0,              # 
           Rfm = 0,              # 
           Hm = 0,               # 
           Sfh = 0,              # initial number of those in state of susceptibility for those on MDA but not protected
           Rfh = 0,              # initial state of recovery for those on MDA but not protected
           Efh = 0,              # 
           Afh = 0,              # 
           Cfh = 0,              # 
           Ctfh = 0,             # 
           Sevh = 0,             # initial number of severe cases after
           CMDA = 0,             # MDA administered at time 0
           CSnap = 0)            # those who are snap back to the original compartment

### 3. Define dynamic Pf Human-static Vector Plasmodium falciparum (pf) model ####
## 3.1. the model with behaviour change ####
behaviour_change.model <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    
    # defining total population
    Pf = Sf+Ef+Af+Cf+Ctf+Tf+Sev+H+Rf+Rpf+Sfm+Tfm+Rfm+Hm+Sfh+Rfh+Efh+Afh+Cfh+Ctfh+Sevh          # the population does not account for death
    
    ## A. INTERVENTION PART ####
    ## A.1. Bednets: Baseline ####
    # ITN effectiveness
    factor = 1
    if (t>360*5) {factor = net.coverage}
    
    itn_cov = itndata$nets/360*ppn.fit/Pf*factor
    itn_t = itndata$time
    itn_dist<-approx(itn_t, itn_cov, t, method="constant")$y
    eta <- -log(1-(1-0.6))/(3*360) #loss due to attrition rate
    
    # plot(1:10000, approx(itn_t, itn_cov, 1:10000,method="constant")$y , ty="l")
    # points(itn_t, itn_cov, col="red")
    
    itn = min(ITN,1)*p_use*p_effectiveness    # gives the cap limit to one as coverage cannot be more than 100%
    
    ## A.2. Seasonal Mass Drug Administration (MDA) ####
    # Several rounds: MDA takes place in January (for 30 days) for 5 360s
    
    rounds = 5            # 5 rounds of MDA campaign (1 round: 30 days in a year)
    mdur = 30             # 30 days
    pulse = c(rep(0, 360*5),
              rep(c(rep(1, mdur), rep(0, (360-mdur))),rounds),
              rep(0, (360*5+rounds*360)))
    
    # par(mfrow = c(1,1))
    # plot(pulse, ty="l")
    
    mrate <- snap <-0
    
    #browser()
    
    # the first 5 years data are discarded, 5 years starts from year-5 to year-10
    if(t> (360*5) & t <= (360*10)) {
      mrate = approx(1:length(pulse), (-log(1-mcov)/mdur)*pulse, t, method="constant", rule=2)$y
      snap = approx(1:length(pulse), (snapr*lag(pulse, default=0, n=mdur)), t, method="constant", rule=2)$y
    }
    
    ## A.3. Behaviour change ####
    # behaviour change intervention: it takes 2 phases
    # the 1st phase:
    int.start<- 360*5         # since the first 5 years data are discarded, so the behaviour change intervention begins at the beginning of the 360
    int.lag <- 360*1          # there is a lag time before an intervention takes an effect (e.g. due to preparation of the program)
    int.target_time <- 360*2  # it is assumed that the target will be reached within 2 years from the point when the effect starts
    
    # the target for bednets usage, 80% of people use bed nets and 80% of people seek for treatment
    # int.target_net
    # int.target_seek
    
    if(t > (int.start + int.lag)){  
      p_use = log(t-int.target_time)/log(int.target_time)*(int.target_net - 0.54) + 0.54  # it follows log function, with the cap limit is 0.8
      pseek = log(t-int.target_time)/log(int.target_time)*(int.target_seek - 0.60) + 0.60 # it follows log function, with the cap limit is 0.8
    }
    
    # plot((log(1:730))/log(730)*(0.8-0.54)+0.54, type = "l")   # this plot is to check the dynamics of behaviour change for nets usage
    # plot((log(1:730))/log(730)*(0.8-0.6)+0.6, type = "l")     # this plot is to check the dynamics of behaviour change for treatment seeking
    
    # the 2nd phase:
    else if (t > (int.start + int.lag + int.target_time)){  # the intervention continuous at the second year
      p_use = int.target_net   # the proportion of net use increases to 80%
      pseek = int.target_seek   # the proportion of patients who seek for treatment increases to 80%
    }
    
    # B. MODEL PART ####
    # treatment cascade
    pi_t = pseek*(ptest_slide*psens_slide + ptest_RDT*psens_RDT)*ptreat
    ptrt = pi_t
    
    # force of infection
    Infectious = Cf+Ctf+Sev+zeta_a*Af+zeta_t*(Tf+H) # Infectious contribution to the population
    seas <- 1+amp*cos(2*pi*(t/360 - phi))^peak           # seasonality, assuming there is no phase angle meaning the peak occurs at the beginning of each cycle cycle
    lambdaf = (1-itn)*seas*(a^2*b*c*m*Infectious/Pf)/(a*c*Infectious/Pf+mu_m)*(gamma_m/(gamma_m+mu_m))
    
    # compartment model
    dSf = mu_h*Pf - lambdaf*Sf + rho_f*Rf - mu_h*Sf - mrate*Sf + snap*Sfh
    dEf = lambdaf*Sf - (gamma_hf + mu_h)*Ef - mrate*Ef + snap*Efh
    dAf = pa_f*gamma_hf*Ef + omega_f*Cf - (delta_f + mu_h)*Af - mrate*Af + snap*Afh
    dCf = (1-pa_f)*(1-ptrt)*gamma_hf*Ef +  epsilon_f*Sev - (omega_f + nu_f + mu_h)*Cf - mrate*Cf + snap*Cfh
    dCtf = (1-pa_f)*ptrt*gamma_hf*Ef - (tau_f+ mu_h)*Ctf - mrate*Ctf + snap*Ctfh
    dSev = nu_f*Cf - (epsilon_f + tau_sev + mu_sev + mu_h)*Sev - mrate*Sev + snap*Sevh
    dTf = tau_f*Ctf - (r_f + mu_h)*Tf
    dH = tau_sev*Sev - (r_hosp+ mu_hosp + mu_h)*H
    dRf = delta_f*Af  - (rho_f + mu_h)*Rf - mrate*Rf + snap*Rfh + kappa*Rpf
    
    # for this MDA model, instead of Tf and H that go to Rf, they go to Rpf first then to Rf, with the rate of change from Rpf to Rf = kappa.
    dRpf = r_f*Tf + r_hosp*H - (kappa + mrate + mu_h)*Rpf
    
    dSfm = mrate*Sf - (mprotect + mu_h)*Sfm
    dTfm = mrate*(Ef+Af+Cf+Ctf) + tau_f*Ctfh - (r_f + mu_h)*Tfm
    dRfm = mrate*(Rpf+Rf) + r_f*Tfm +r_hosp*Hm - (mprotect + mu_h)*Rfm
    dHm = mrate*Sev + tau_sev*Sevh - (r_hosp + mu_hosp + mu_h)*Hm   #
    
    dSfh = mprotect*Sfm - (lambdaf + snap + mu_h)*Sfh # state of susceptibility for those on MDA but not protected
    dRfh = mprotect*Rfm + delta_f*Afh - (snap + mu_h)*Rfh # state of Recovered for those on MDA but not protected
    dEfh = lambdaf*Sfh - (gamma_hf + snap + mu_h)*Efh
    dAfh = pa_f*gamma_hf*Efh + omega_f*Cfh - (delta_f + snap + mu_h)*Afh
    dCfh = (1-pa_f)*(1-ptrt)*gamma_hf*Efh + epsilon_f*Sevh - (omega_f + nu_f + snap + mu_h)*Cfh
    dCtfh = (1-pa_f)*ptrt*gamma_hf*Efh - (tau_f + snap + mu_h)*Ctfh
    dSevh = nu_f*Cfh - (epsilon_f + tau_sev + snap + mu_sev + mu_h)*Sevh
    
    # Counters
    dD = mu_sev*(Sev+Sevh) + mu_hosp*(H+Hm)   # death counter, as now death also comes from those in the MDA compartments (Sevh and Hm)
    
    dCinc = lambdaf*(Sf+Sfh)             # cumulative incidence
    dCtrt = tau_f*(Ctf+Ctfh)             # cumulative treated patients with uncomplicated malaria
    dChosp = tau_sev*(Sev+Sevh)          # cumulative hospitalised patients
    dCtrt_total = tau_f*(Ctf+Ctfh) + tau_sev*(Sev+Sevh)   # cumulative treated patients by all means
    
    dCMDA = mrate*(Sf+Ef+Af+Cf+Rpf+Rf+Ctf+Sev)
    dCsnap = snap*(Sfh+Rfh+Efh+Afh+Cfh+Ctfh+Sevh)
    
    dCtest = tau_f*(Ctf + Ctfh)/ptreat/(ptest_RDT*RDT_pos + ptest_slide*slide_pos) + # cumulative number of tests
      tau_sev*(Sev + Sevh)/1/(ptest_RDT*RDT_pos + ptest_slide*slide_pos) # all severe cases are assumed to be treated
    
    # ITN
    dITN = itn_dist - (eta + mu_net)*ITN # coverage of the population with potential to be protected by nets CURRENTLY IN CIRCULATION
    
    # returning the rate of change
    output <- c(dSf, dEf, dAf, dCf, dCtf, dTf, dSev, dH, dRf, dRpf, dCinc, dCtrt, dChosp, dCtrt_total, dCtest, dITN, dD,
                dSfm, dTfm, dRfm, dHm, dSfh, dRfh, dEfh, dAfh, dCfh, dCtfh, dSevh, dCMDA, dCsnap)
    list(output)
  })
}

## 3.2. the model without behaviour change ####
no_behaviour_change.model <- function(t,x,parms){
  parms['int.target_net'] <- parms['p_use']
  parms['int.target_seek'] <- parms['pseek']
  behaviour_change.model(t,x,parms)
}

## 3.3. Run models ####
ncov1 <- parms['net.coverage'] <- 1
mcov1 <- parms['mcov'] <- 0

hm.0 <- ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)

# 3.4. check the basic model ####
# mutate
df.0 <- as_tibble(as.data.frame(hm.0)) %>% 
  mutate(P = Sf+Ef+Af+Cf+Ctf+Sev+Tf+H+Rf+Rpf+Sfm+Tfm+Rfm+Hm+Sfh+Rfh+Efh+Afh+Cfh+Ctfh+Sevh,
         P_total = P + Death,
         C_total = Cf + Ctf + Cfh + Ctfh+ Sev,
         T_total = Tf + Tfm + H+ Hm,
         Deathf = c(0, diff(Death)),
         Incf = c(0, diff(Cinc)),
         Trtf = c(0, diff(Ctrt))) %>% 
  pivot_longer(names_to = "variable", cols = !1)

# Population check
df.0  %>% 
  filter(variable %in% c("P_total")) %>% 
  ggplot()+
  geom_line(aes(x = time, y=value))+
  theme_minimal() +
  labs(title = "Populations", y =("population"))

# run the models again 
hm.behaviour_change <- ode(times=times, y=start, func=behaviour_change.model,parms=parms)
hm.no_behaviour_change <- ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)

### 4. health benefit analysis (HBA) ####
## 4.1. HBA in scenarios with behaviour change ####
# preparing for an array to contain the outputs
netcoverage.grid <- seq(0, 1, 0.1)
mdacoverage.grid <- seq(0, 0.9, 0.1)  # as in the function, 100% coverage would result in inf, so here the max coverage is set to be 0.9
out.name <- c("nhb", "nmb", "cost", "daly", "diff.cost", "daly.gained")

dim.net <- length(netcoverage.grid)
dim.mda <- length(mdacoverage.grid)
dim.out <- length(out.name)

nhb.grid <- array(NA, 
                  dim=c(dim.net, dim.mda, dim.out), 
                  dimnames = list(netcoverage.grid, mdacoverage.grid, out.name))

hba.bc <- function(x) {
  # preparing for an array to contain the outputs ####
  netcoverage.grid <- seq(0, 1, 0.1)
  mdacoverage.grid <- seq(0, 0.9, 0.1)  # as in the function, 100% coverage would result in inf, so here the max coverage is set to be 0.9
  out.name <- c("nhb", "nmb", "cost", "daly", "diff.cost", "daly.gained")
  
  dim.net <- length(netcoverage.grid)
  dim.mda <- length(mdacoverage.grid)
  dim.out <- length(out.name)
  
  nhb.grid <- array(NA, 
                    dim=c(dim.net, dim.mda, dim.out), 
                    dimnames = list(netcoverage.grid, mdacoverage.grid, out.name))
  
  # setting up for threshold for health benefit measurement
  indonesian.gdp.per.capita <- 3820 # https://data.worldbank.org/indicator/NY.GDP.PCAP.KD?locations=ID
  threshold <- 3*indonesian.gdp.per.capita
  
  # defining variable and parameters for baseline ####
  ncov <- parms['net.coverage'] <- 1
  parms['mcov'] <- 0
  
  hm <-ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)
  
  # sum uncomplicated cases
  uncomplicated <- sum(hm[x, ('Cf')] + hm[x,'Ctf'] + hm[x,'Cfh'] + hm[x,'Ctfh'])
  
  # sum severe cases
  complicated <- sum(hm[x, ('Sev')] + hm[x,'Sevh'])
  
  # sum deaths
  death <- sum(diff(hm[x, 'Death']))
  
  # sum hospitalised
  hospitalised <- sum(diff(hm[x, 'Chosp']))
  
  # sum outpatients
  outpatient <- sum(diff(hm[x, 'Ctrt']))
  
  # sum tests
  test <- sum(diff(hm[x, 'Ctest']))
  
  # sum nets
  net <- sum(itndata$nets[5:9])
  
  # sum mda
  mda <- sum(diff(hm[x, 'CMDA']))
  
  # Health Cost ####
  # treatment cost: hospitalised*(33.5) + outpatient*(5.5 + 0.293)
  treatment.cost <- hospitalised*(33.5) + outpatient*(5.5 + 0.293)
  
  # test cost: test*(1.2)
  test.cost <- test*(1.2)
  
  # bed-net cost: nets*(5+2)
  net.cost <- net*(5+2)
  
  # mda cost: mda*(5.5 + 0.293)
  mda.cost <- mda*(5.5 + 0.293)
  
  # health system cost: total_population*no_of_year*(4.5)
  health.system.cost <- sum(start)*5*(4.5)
  
  # total cost: accounting 3% discount rate for 5 years as this is the total investment in 2017
  total.cost.baseline <- (treatment.cost + test.cost + net.cost*ncov+ mda.cost + health.system.cost)/((1+0.3)^5)
  
  # DALY ####
  # DALY: uncomplicated*(0.211) + complicated*(0.436) + death*(50) [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3991722/#!po=36.1111]
  daly.baseline <- uncomplicated*(0.211) + complicated*(0.436) + death*(50)
  
  # sensitivity analysis ####
  for(i in 1:dim.net){
    for(j in 1:dim.mda){
      # defining variable and parameters ####
      ncov <- parms['net.coverage'] <- netcoverage.grid[i]
      parms['mcov'] <- mdacoverage.grid[j]
      
      hm <-ode(times=times, y=start, func=behaviour_change.model,parms=parms)
      
      # sum uncomplicated cases
      uncomplicated <- sum(hm[x, ('Cf')] + hm[x,'Ctf'] + hm[x,'Cfh'] + hm[x,'Ctfh'])
      
      # sum severe cases
      complicated <- sum(hm[x, ('Sev')] + hm[x,'Sevh'])
      
      # sum deaths
      death <- sum(diff(hm[x, 'Death']))
      
      # sum hospitalised
      hospitalised <- sum(diff(hm[x, 'Chosp']))
      
      # sum outpatients
      outpatient <- sum(diff(hm[x, 'Ctrt']))
      
      # sum tests
      test <- sum(diff(hm[x, 'Ctest']))
      
      # sum nets
      net <- sum(itndata$nets[5:9])
      
      # sum mda
      mda <- sum(diff(hm[x, 'CMDA']))
      
      # Health Cost ####
      # treatment cost: hospitalised*(33.5) + outpatient*(5.5 + 0.293)
      treatment.cost <- hospitalised*(33.5) + outpatient*(5.5 + 0.293)
      
      # test cost: test*(1.2)
      test.cost <- test*(1.2)
      
      # bed-net cost: nets*(5+2)
      net.cost <- net*(5+2)
      
      # mda cost: mda*(5.5 + 0.293)
      mda.cost <- mda*(5.5 + 0.293)
      
      # behaviour change cost: no of village * cost per village per year * year = 8747*500*5
      bc.cost <- parms['bc']
      
      # health system cost: total_population*no_of_year*(4.5)
      health.system.cost <- sum(start)*5*(4.5)
      
      # total cost: accounting 3% discount rate for 5 years as this is the total investment in 2017
      total.cost <- (treatment.cost + test.cost + net.cost*ncov + mda.cost + bc.cost + health.system.cost)/((1+0.3)^5)
      
      # cost difference
      diff.cost <- total.cost - total.cost.baseline
      
      # DALY ####
      # DALY: uncomplicated*(0.211) + complicated*(0.436) + death*(50) [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3991722/#!po=36.1111]
      daly <- uncomplicated*(0.211) + complicated*(0.436) + death*(50)
      
      # daly gained
      daly.gained <- daly.baseline - daly
      
      # HEALTH BENEFIT MEASUREMENT ####
      # Net health benefit (NHB): daly gained - total cost/threshold
      nhb <- daly.gained - (diff.cost)/threshold
      nmb <- nhb*threshold
      
      nhb.grid[i,j,1] <- nhb
      nhb.grid[i,j,2] <- nmb
      nhb.grid[i,j,3] <- total.cost
      nhb.grid[i,j,4] <- daly
      nhb.grid[i,j,5] <- diff.cost
      nhb.grid[i,j,6] <- daly.gained
      
      print(c("total cost baseline" = total.cost.baseline))
      print(c("nhb" = nhb))
      print(c('net' = parms['net.coverage'], 'mda' = parms['mcov'], 'cost' = total.cost, 'diff.cost' = diff.cost, 'daly.gained' = daly.gained))
    }
  }
  return(nhb.grid)
}

## 4.2. HBA in scenarios without behaviour change ####
hba.nbc <- function(x) {
  # preparing for an array to contain the outputs ####
  netcoverage.grid <- seq(0, 1, 0.1)
  mdacoverage.grid <- seq(0, 0.9, 0.1)  # as in the function, 100% coverage would result in inf, so here the max coverage is set to be 0.9
  out.name <- c("nhb", "nmb", "cost", "daly", "diff.cost", "daly.gained")
  
  dim.net <- length(netcoverage.grid)
  dim.mda <- length(mdacoverage.grid)
  dim.out <- length(out.name)
  
  nhb.grid <- array(NA, 
                    dim=c(dim.net, dim.mda, dim.out), 
                    dimnames = list(netcoverage.grid, mdacoverage.grid, out.name))
  
  # setting up for threshold for health benefit measurement
  indonesian.gdp.per.capita.2017 <- 3820
  threshold <- 3*indonesian.gdp.per.capita.2017 # WHO threshold for a cost-effective intervention (no more than 3 gdp)
  
  # defining variable and parameters for baseline ####
  ncov <- parms['net.coverage'] <- 1
  parms['mcov'] <- 0
  
  hm <-ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)
  
  # sum uncomplicated cases
  uncomplicated <- sum(hm[x, ('Cf')] + hm[x,'Ctf'] + hm[x,'Cfh'] + hm[x,'Ctfh'])
  
  # sum severe cases
  complicated <- sum(hm[x, ('Sev')] + hm[x,'Sevh'])
  
  # sum deaths
  death <- sum(diff(hm[x, 'Death']))
  
  # sum hospitalised
  hospitalised <- sum(diff(hm[x, 'Chosp']))
  
  # sum outpatients
  outpatient <- sum(diff(hm[x, 'Ctrt']))
  
  # sum tests
  test <- sum(diff(hm[x, 'Ctest']))
  
  # sum nets
  net <- sum(itndata$nets[5:9])
  
  # sum mda
  mda <- sum(diff(hm[x, 'CMDA']))
  
  # Health Cost ####
  # treatment cost: hospitalised*(33.5) + outpatient*(5.5 + 0.293)
  treatment.cost <- hospitalised*(33.5) + outpatient*(5.5 + 0.293)
  
  # test cost: test*(1.2)
  test.cost <- test*(1.2)
  
  # bed-net cost: nets*(5+2)
  net.cost <- net*(5+2)
  
  # mda cost: mda*(5.5 + 0.293)
  mda.cost <- mda*(5.5 + 0.293)
  
  # health system cost: total_population*no_of_year*(4.5)
  health.system.cost <- sum(start)*5*(4.5)
  
  # total cost: accounting 3% discount rate for 5 years as this is the total investment in 2017
  total.cost.baseline <- (treatment.cost + test.cost + net.cost*ncov + mda.cost + health.system.cost)/((1+0.3)^5)
  
  # DALY ####
  # DALY: uncomplicated*(0.211) + complicated*(0.436) + death*(50) [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3991722/#!po=36.1111]
  daly.baseline <- uncomplicated*(0.211) + complicated*(0.436) + death*(50)
  
  # sensitivity analysis ####
  for(i in 1:dim.net){
    for(j in 1:dim.mda){
      # defining variable and parameters ####
      ncov <- parms['net.coverage'] <- netcoverage.grid[i]
      parms['mcov'] <- mdacoverage.grid[j]
      
      hm <-ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)
      
      # sum uncomplicated cases
      uncomplicated <- sum(hm[x, ('Cf')] + hm[x,'Ctf'] + hm[x,'Cfh'] + hm[x,'Ctfh'])
      
      # sum severe cases
      complicated <- sum(hm[x, ('Sev')] + hm[x,'Sevh'])
      
      # sum deaths
      death <- sum(diff(hm[x, 'Death']))
      
      # sum hospitalised
      hospitalised <- sum(diff(hm[x, 'Chosp']))
      
      # sum outpatients
      outpatient <- sum(diff(hm[x, 'Ctrt']))
      
      # sum tests
      test <- sum(diff(hm[x, 'Ctest']))
      
      # sum nets
      net <- sum(itndata$nets[5:9])
      
      # sum mda
      mda <- sum(diff(hm[x, 'CMDA']))
      
      # Health Cost ####
      # treatment cost: hospitalised*(33.5) + outpatient*(5.5 + 0.293)
      treatment.cost <- hospitalised*(33.5) + outpatient*(5.5 + 0.293)
      
      # test cost: test*(1.2)
      test.cost <- test*(1.2)
      
      # bed-net cost: nets*(5+2)
      net.cost <- net*(5+2)
      
      # mda cost: mda*(5.5 + 0.293)
      mda.cost <- mda*(5.5 + 0.293)
      
      # behaviour change cost: no of village * cost per village per year * year = 8795*500*5
      bc.cost <- 0
      
      # health system cost: total_population*no_of_year*(4.5)
      health.system.cost <- sum(start)*5*(4.5)
      
      # total cost: accounting 3% discount rate for 5 years as this is the total investment in 2017
      total.cost <- (treatment.cost + test.cost + net.cost*ncov + mda.cost + bc.cost + health.system.cost)/((1+0.3)^5)
      
      # cost difference
      diff.cost <- total.cost - total.cost.baseline
      
      # DALY ####
      # DALY: uncomplicated*(0.211) + complicated*(0.436) + death*(50) [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3991722/#!po=36.1111]
      daly <- uncomplicated*(0.211) + complicated*(0.436) + death*(50)
      
      # daly gained
      daly.gained <- daly.baseline - daly
      
      # HEALTH BENEFIT MEASUREMENT ####
      # Net health benefit (NHB): daly gained - total cost/threshold
      nhb <- daly.gained - (diff.cost)/threshold
      nmb <- nhb*threshold
      
      nhb.grid[i,j,1] <- nhb
      nhb.grid[i,j,2] <- nmb
      nhb.grid[i,j,3] <- total.cost
      nhb.grid[i,j,4] <- daly
      nhb.grid[i,j,5] <- diff.cost
      nhb.grid[i,j,6] <- daly.gained
      
      print(c("total cost baseline" = total.cost.baseline))
      print(c("nhb" = nhb))
      print(c('net' = parms['net.coverage'], 'mda' = parms['mcov'], 'cost' = total.cost, 'diff.cost' = diff.cost, 'daly.gained' = daly.gained))
    }
  }
  return(nhb.grid)
}

## 4.3. Running the HBA ####
interval = 61:121   # running from the start of 2017 to the end of 2021

start.time <- Sys.time()

# running the function
nhb.grid.bc <- hba.bc(interval)
nhb.grid.nbc <- hba.nbc(interval)

# time needed for running the code
finish.time <- Sys.time()
finish.time - start.time

# check the output for net health benefit
nhb.grid.bc[,,'nhb']
nhb.grid.nbc[,,'nhb']

## 4.4. Analysing which one has the highest net health benefit (NHB): without behaviour change ####
# 4.4.1. net health benefit analysis for behaviour change intervention ####
nbc.nhb <- matrix(0, nrow = dim.net, ncol = dim.mda)
rownames(nbc.nhb) <- as.character(netcoverage.grid)
colnames(nbc.nhb) <-as.character(mdacoverage.grid)
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if (nhb.grid.nbc[i,j,'nhb'] <= 0) {nbc.nhb[i,j] <- 0}
    else if (0 < nhb.grid.nbc[i,j,'nhb'] && nhb.grid.nbc[i,j,'nhb'] <= 50000) {nbc.nhb[i,j] <- 1}
    else if (50000 < nhb.grid.nbc[i,j,'nhb'] && nhb.grid.nbc[i,j,'nhb'] <= 100000) {nbc.nhb[i,j] <- 2}
    else {nbc.nhb[i,j] <- 3}
  }
}

# 4.4.2. total cost analysis for behaviour change intervention ####
nbc.cost <- matrix(0, nrow = dim.net, ncol = dim.mda) 
rownames(nbc.cost) <- netcoverage.grid
colnames(nbc.cost) <-mdacoverage.grid
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if (nhb.grid.nbc[i,j,'cost'] <= (nhb.grid.nbc[11,1,'cost'])) {nbc.cost[i,j] <- 1}
    else {nbc.cost[i,j] <- 0}
  }
}
nbc.cost

# 4.4.3. nhb under budget constraints analysis for behaviour change intervention ####
nbc.nhb.cost.accepted <- matrix(0, nrow = dim.net, ncol = dim.mda)
rownames(nbc.nhb.cost.accepted) <- netcoverage.grid
colnames(nbc.nhb.cost.accepted) <-mdacoverage.grid
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if (nbc.cost[i,j] == 1 & nbc.nhb[i,j] > 0) {nbc.nhb.cost.accepted[i,j] <- nhb.grid.nbc[i,j,'nhb']}
    else {nbc.nhb.cost.accepted[i,j] <- 0}
  }
}
nbc.nhb.cost.accepted

## 4.5. Analysing which one has the highest net health benefit (NHB): with behaviour change ####
# 4.5.1. net health benefit analysis for behaviour change intervention ####
bc.nhb <- matrix(0, nrow = dim.net, ncol = dim.mda)
rownames(bc.nhb) <- netcoverage.grid
colnames(bc.nhb) <-mdacoverage.grid
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if (nhb.grid.bc[i,j,'nhb'] <= 0) {bc.nhb[i,j] <- 0}
    else if (0 < nhb.grid.bc[i,j,'nhb'] && nhb.grid.bc[i,j,'nhb'] <= 50000) {bc.nhb[i,j] <- 1}
    else if (50000 < nhb.grid.bc[i,j,'nhb'] && nhb.grid.bc[i,j,'nhb'] <= 100000) {bc.nhb[i,j] <- 2}
    else {bc.nhb[i,j] <- 3}
  }
}

# 4.5.2. total cost analysis for behaviour change intervention ####
bc.cost <- matrix(0, nrow = dim.net, ncol = dim.mda) 
rownames(bc.cost) <- netcoverage.grid
colnames(bc.cost) <-mdacoverage.grid
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if (nhb.grid.bc[i,j,'cost'] <= (nhb.grid.bc[11,1,'cost'])) {bc.cost[i,j] <- 1}
    else {bc.cost[i,j] <- 0}
  }
}
bc.cost

# 4.5.3. nhb under budget constraints analysis for behaviour change intervention ####
bc.nhb.cost.accepted <- matrix(0, nrow = dim.net, ncol = dim.mda)
rownames(bc.nhb.cost.accepted) <- netcoverage.grid
colnames(bc.nhb.cost.accepted) <-mdacoverage.grid
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if (bc.cost[i,j] == 1 & bc.nhb[i,j] > 0) {bc.nhb.cost.accepted[i,j] <- nhb.grid.bc[i,j,'nhb']}
    else {bc.nhb.cost.accepted[i,j] <- 0}
  }
}

bc.nhb.cost.accepted

### 5. Plots ####
# WITHOUT BEHAVIOUR CHANGE ####
par(mfrow = c(2,3))
image(netcoverage.grid, mdacoverage.grid, nbc.nhb, col = hcl.colors(4, palette = "YlOrRd"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Net Health Benefit: No Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -nbc.cost, col = hcl.colors(2, palette = "greens"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Cost: No Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -nbc.nhb.cost.accepted, col = hcl.colors(10, palette = "blues"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "NHB and Cost: No Behaviour Change")

# WITH BEHAVIOUR CHANGE ####
image(netcoverage.grid, mdacoverage.grid, bc.nhb, col = hcl.colors(4, palette = "YlOrRd"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Net Health Benefit: With Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -bc.cost, col = hcl.colors(2, palette = "greens"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Cost: With Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -bc.nhb.cost.accepted, col = hcl.colors(10, palette = "blues"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "NHB and Cost: With Behaviour Change")
