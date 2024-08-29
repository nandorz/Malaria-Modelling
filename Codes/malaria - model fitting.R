# CHAPTER 2: MODEL FITTING #####################################################
################################################################################
### 1. Model Fitting using maximum likelihood method to find the best fit for S0 and ppn ####
# S0 = initial susceptible
# ppn = person per net

time <- 0                              # set initial time

t_start<- 0
n <- 10
t_end <- n*360
step <- 1

times <- seq(t_start, t_end, step)

## 1.1. Parameters ####
falciparum <- read_excel("Parameters and initial condition - falciparum.xlsx", sheet = "default") # loading excel data
parameters_value.vector <- c(falciparum$value)  # taking the value of parameters from data frame into a vector
parameters_name.vector <- c(falciparum$name)    # taking the name of parameters from data frame into a vector

names(parameters_value.vector) <- parameters_name.vector  # giving the parameters value with name

parms <- parameters_value.vector

## 1.2. ITN data ####
itndata <- read_excel("itndata.xlsx")

## 1.3. Initial conditions ####
initSf <- 5000000  # initial susceptible population
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

initITN <- 0     # assigning the initial number of ITN

start <- c(Sf = initSf, 
           Ef = initEf, 
           Af = initAf, 
           Cf = initCf, 
           Ctf = initCtf,
           Tf = initTf, 
           Sev = initSev,
           H = initH,
           Rf = initRf, 
           Rpf = 0,
           Cinc = initCinc, 
           Ctrt = initCtr,
           Chosp = 0,
           Ctrt_total = initCtrt_total,
           Ctest = initCtest,
           ITN = initITN,
           Death = 0,
           Sfm = 0,
           Tfm = 0,
           Rfm = 0,
           Hm = 0,
           Sfh = 0,
           Rfh = 0,
           Efh = 0,
           Afh = 0,
           Cfh = 0,
           Ctfh = 0,
           Sevh = 0,
           CMDA = 0,
           CSnap = 0)

### 2. Define dynamic Pf Human-static Vector Plasmodium falciparum (pf) model ####
## 2.1. the model with behaviour change ####
behaviour_change.model <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    
    # defining total population
    Pf = Sf+Ef+Af+Cf+Ctf+Tf+Sev+H+Rf+Rpf+Sfm+Tfm+Rfm+Hm+Sfh+Rfh+Efh+Afh+Cfh+Ctfh+Sevh          # the population does not account for death
    
    ## A. INTERVENTION PART ####
    ## A.1. Bednets: Baseline ####
    # ITN effectiveness
    factor = 1
    if (t>360*5) {factor = net.coverage}
    
    itn_cov = itndata$nets/360*ppn/Pf*factor
    itn_t = itndata$time
    itn_dist<-approx(itn_t, itn_cov, t, method="constant")$y
    eta <- -log(1-(1-0.6))/(3*360) #loss due to attrition rate
    
    # plot(1:10000, approx(itn_t, itn_cov, 1:10000,method="constant")$y , ty="l")
    # points(itn_t, itn_cov, col="red")
    
    itn = min(ITN,1)*p_use*p_effectiveness    # gives the cap limit to one as coverage cannot be more than 100%
    
    ## A.2. Seasonal Mass Drug Administration (MDA) #### not necessary for model fitting
    # Several rounds: MDA takes place in January (for 30 days) for 5 360s
    
    # rounds = 5            # 5 rounds of MDA campaign (1 round: 30 days in a year)
    # mdur = 30             # 30 days
    # pulse = c(rep(0, 360*5),
    #          rep(c(rep(1, mdur), rep(0, (360-mdur))),rounds),
    #          rep(0, (360*5+rounds*360)))
    
    # par(mfrow = c(1,1))
    # plot(pulse, ty="l")
    
    mrate <- snap <-0
    
    #browser()
    
    # the first 5 years data are discarded, 5 years starts from year-5 to year-10
    # if(t> (360*5) & t <= (360*10)) {
    #  mrate = approx(1:length(pulse), (-log(1-mcov)/mdur)*pulse, t, method="constant", rule=2)$y
    #  snap = approx(1:length(pulse), (snapr*lag(pulse, default=0, n=mdur)), t, method="constant", rule=2)$y
    #}
    
    ## A.3. Behaviour change #### not necessary for model fitting
    # behaviour change intervention: it takes 2 phases
    # the 1st phase:
    # int.start<- 360*5         # since the first 5 years data are discarded, so the behaviour change intervention begins at the beginning of the 360
    # int.lag <- 360*1          # there is a lag time before an intervention takes an effect (e.g. due to preparation of the program)
    # int.target_time <- 360*2  # it is assumed that the target will be reached within 2 years from the point when the effect starts
    
    # the target for bednets usage, 80% of people use bed nets and 80% of people seek for treatment
    # int.target_net
    # int.target_seek
    
    # if(t > (int.start + int.lag)){  
    #  p_use = log(t-int.target_time)/log(int.target_time)*(int.target_net - 0.54) + 0.54  # it follows log function, with the cap limit is 0.8
    #  pseek = log(t-int.target_time)/log(int.target_time)*(int.target_seek - 0.60) + 0.60 # it follows log function, with the cap limit is 0.8
    #}
    
    # plot((log(1:730))/log(730)*(0.8-0.54)+0.54, type = "l")   # this plot is to check the dynamics of behaviour change for nets usage
    # plot((log(1:730))/log(730)*(0.8-0.6)+0.6, type = "l")     # this plot is to check the dynamics of behaviour change for treatment seeking
    
    # the 2nd phase:
    # else if (t > (int.start + int.lag + int.target_time)){  # the intervention continuous at the second year
    # p_use = int.target_net   # the proportion of net use increases to 80%
    # pseek = int.target_seek   # the proportion of patients who seek for treatment increases to 80%
    #}
    
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

## 2.2. the model without behaviour change ####
no_behaviour_change.model <- function(t,x,parms){
  parms['int.target_net'] <- parms['p_use']
  parms['int.target_seek'] <- parms['pseek']
  behaviour_change.model(t,x,parms)
}

## 2.3. Run models ####
ncov1 <- parms['net.coverage'] <- 1
mcov1 <- parms['mcov'] <- 0

hm.0 <- ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)

### 3. data generated ####
pf.cases <- c(b*pf.incidence$'2017', b*pf.incidence$'2018', b*pf.incidence$'2019', b*pf.incidence$'2020', b*pf.incidence$'2021')

N1 <- pf.cases

### 4. creating a function to calculate negative likelihood (NLL) ####
calcNLL <- function(params){
  # returns NLL assuming binomial observations with N trials and d successes  
  # with initial 'Sf' and 'ppn' specified as the first 2 elements of the params vector
  
  # constructing the model 
  start['Sf'] <- exp(params[1])
  parms['ppn'] <- params[2]    # this value is also going to be fitted
  
  model_output.df <- data.frame(lsoda(
    y = start,               # Initial conditions for population
    times = times,             # Timepoints for evaluation
    func = no_behaviour_change.model,                   # Function to evaluate
    parms = parms                 # Vector of parameters
  ))
  
  ## construct the new infection: lambda
  # extracting the number of row for the model output
  num.timepoints <- dim(model_output.df)[1]
  
  # taking the difference of the cumulative incidence returns to the incidence
  new.infections<-round(c(model_output.df$Cinc[1],
                          model_output.df$Cinc[2:num.timepoints] - model_output.df$Cinc[1:(num.timepoints-1)] ))
  lambda <- new.infections
  
  ## returning the output NLL value given the inputs
  # the output that are fitted are those in equilibrium. It is assumed that the dynamics of the model becomes stable after 5 years, so the first 5 years are discarded
  # then, the the following 5 years are fitted with the past 5 years of incidence report
  NLL <- -sum(dpois(round(N1),  lambda[(360*5+1):(360*10)], log = TRUE))
  return(NLL)           
}

# 5. Fitting the model using optim ####
# running the optim function from certain starting values of params[1] and params[2]
set.seed(123)
fit1<-optim(c(log(5e6), 2),
            calcNLL,
            control=list(reltol=1e-12), 
            hessian = T)
fit1

# backtransforming the parameters
S0.fit <- exp(fit1$par[1])
ppn.fit <- fit1$par[2]


# calculating 95%CI
sds.for.parameters <- sqrt(diag(solve(fit1$hessian)))

S0.fit_lower.95CI <- S0.fit - 1.96*exp(sds.for.parameters[1])
S0.fit_upper.95CI <- S0.fit + 1.96*exp(sds.for.parameters[1])

ppn_lower.95CI <- ppn.fit - 1.96*sds.for.parameters[2]
ppn_upper.95CI <- ppn.fit + 1.96*sds.for.parameters[2]

# initial Sf value
S0.fit
S0.fit_lower.95CI
S0.fit_upper.95CI

# person per net
ppn.fit
ppn_lower.95CI 
ppn_upper.95CI

# initial Sf value
S0.fit

# person per net
ppn.fit

## 5. contour plot ####

S0.grid <- seq(4000000, 9000000, 100000)
ppn.grid <- seq(1, 5, 0.2)
dim1 <- length(S0.grid)
dim2 <- length(ppn.grid)
NLL.grid <- matrix(NA, nrow = dim1, ncol = dim2)

for(i in 1:dim1){
  for(j in 1:dim2){
    NLL.grid[i,j] <- calcNLL(c(log(S0.grid[i]), ppn.grid[j]))
    print(i)
    print(j)
    print(S0.grid[i])
    print(ppn.grid[j])
    print(NLL.grid[i,j])
  }
}

#NLL.df <- as.data.frame(NLL.grid)
#colnames(NLL.grid) <- ppn.grid
#rownames(NLL.grid) <- S0.grid

S0.position <- c()
ppn.position <- c()
for(i in 1:dim1){
  for(j in 1:dim2){
    if(NLL.grid[i,j] == min(NLL.grid)) {
      best.fit.S0 <- S0.grid[i]
      best.fit.ppn <- ppn.grid[j]
      S0.position <- i
      ppn.position <- j
    }
  }
}

par(mfrow = c(1,2))
contour(NLL.grid, nlevels = 250,x = S0.grid, y = ppn.grid, 
        xlab = "Initial susceptible individuals", ylab = "Person per net", main = "Contour Plot: Grid Search for Maximum Likelihood")
points(x = best.fit.S0, y = best.fit.ppn, pch= 20, col = "red", cex = 2)
points(x = S0.fit, y = ppn.fit, pch= 20, col = "green", cex = 2)


image(S0.grid, ppn.grid, -log(NLL.grid), col = heat.colors(10000),
      xlab="Initial susceptible individuals",ylab="Person per net", main = "Heat Map: Grid Search for Maximum Likelihood")
points(x = best.fit.S0, y = best.fit.ppn, pch= 20, col = "red", cex = 2)
points(x = S0.fit, y = ppn.fit, pch= 20, col = "green", cex = 2)


best.fit.S0
best.fit.ppn

S0.fit
ppn.fit

# 6. Model check  ####
ncov1 <- parms['net.coverage'] <- 1
mcov1 <- parms['mcov'] <- 0
start['Sf'] <- S0.fit
parms['ppn'] <- ppn.fit

hm.fit <- ode(times=times, y=start, func=no_behaviour_change.model,parms=parms)

df.fit <- as_tibble(as.data.frame(hm.fit)) %>% 
  mutate(P = Sf+Ef+Af+Cf+Ctf+Sev+Tf+H+Rf+Rpf+Sfm+Tfm+Rfm+Hm+Sfh+Rfh+Efh+Afh+Cfh+Ctfh+Sevh,
         P_total = P + Death,
         C_total = Cf + Ctf + Sev,
         T_total = Tf + H,
         Deathf = c(0, diff(Death)),
         Incf = c(0, diff(Cinc)),
         Trtf = c(0, diff(Ctrt))) %>% 
  pivot_longer(names_to = "variable", cols = !1)

# Population check
df.fit  %>% 
  filter(variable %in% c("P_total")) %>% 
  ggplot()+
  geom_line(aes(x = time, y=value))+
  theme_minimal() +
  labs(title = "Populations", y =("population"))

# Incidence and nets ####
plot_incidence.fit <- df.fit %>%
  filter(variable %in% c("Incf")) %>% 
  filter(time > 360*5) %>% 
  ggplot()+
  geom_line(aes(x = time, y=value), colour = "cyan")+
  geom_point(aes(x = time, y = pf.cases), colour = "orange") + 
  theme_minimal() +
  labs(title = "Incidence - Baseline", y =("population")) 

plot_incidence.fit