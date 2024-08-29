# CHAPTER 4: PLOTTING ##########################################################
################################################################################

library(pacman)
p_load(deSolve, tidyverse, gridExtra, readxl, ggplot2)

# model preparation for plotting ####
# 1.1. setting up parameters for each scenario ####
scn <- c("scenario 3", "scenario 2", "scenario 1", "baseline", "reverse 1", "reverse 2")
date <- seq(as.Date("2011/12/1"), as.Date("2021/12/1"), "month")
ncov3 <- c()
mcov3 <- c()
ncov2 <- c()
mcov2 <- c()

# indexing which combination of net coverage and mda coverage has maximum 
# net health benefit under budget constraint
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if(nbc.nhb.cost.accepted[i,j] == max(nbc.nhb.cost.accepted)){      
      ncov3 <-parms['net.coverage'] <- netcoverage.grid[i]
      mcov3 <-parms['mcov'] <- mdacoverage.grid[j]
    }
  }
}

# indexing which combination of net coverage and mda coverage has maximum 
# net health benefit without budget constraint
for(i in 1:dim.net){
  for(j in 1:dim.mda){
    if(bc.nhb.cost.accepted[i,j] == max(bc.nhb.cost.accepted)){      
      ncov2 <-parms['net.coverage'] <- netcoverage.grid[i]
      mcov2 <-parms['mcov'] <- mdacoverage.grid[j]
    }
  }
}

ncov <- c(ncov3, ncov2, 1, 1, 0.5, 0)
mcov <- c(mcov3, mcov2, 0, 0, 0, 0)

# 1.2. seting up time ####
t_start<- 0
n <- 10
t_end <- n*360
step <- 30

times <- seq(t_start, t_end, step)

# 1.3. create empty data frames for combining data frame of each scenario ####
date.df <- data.frame('date' = date)
pop.check.all <- data.frame()
all.cases.all <- data.frame()
clinical.cases.all <- data.frame()
treated.cases.all <- data.frame()
death.all <- data.frame()
incidence.all <- data.frame()

# 1.4. running all scenarios
for(i in 1:length(scn)){
  # defining net and mda coverage for each parameter in each scenario (scn)
  parms['net.coverage'] <- ncov[i]
  parms['mcov'] <- mcov[i]
  
  # running the output
  if(scn[i] == "scenario 3" | scn[i] == "scenario 1"){
    out <- ode(times=times, y=start, func=behaviour_change.model,parms=parms)
  }
  else (out <- ode(times=times, y=start, func=no_behaviour_change.model,parms=parms))
  
  # binding date data frame into output
  out <- cbind(out, date.df)
  
  # mutate the output
  out.df <-as_tibble(as.data.frame(out)) %>% 
    mutate(P = Sf+Ef+Af+Cf+Ctf+Sev+Tf+H+Rf+Rpf+Sfm+Tfm+Rfm+Hm+Sfh+Rfh+Efh+Afh+Cfh+Ctfh+Sevh,
           P_total = P + Death,
           C_total = Cf + Ctf + Cfh + Ctfh + Sev,
           T_total = Tf + Tfm + H + Hm,
           Deathf = c(0, diff(Death)),
           MDA = c(0, diff(CMDA)),
           Snap = c(0, diff(CSnap)),
           Hosp = c(0, diff(Chosp)),
           Test = c(0, diff(Ctest)),
           Incf = c(0, diff(Cinc)),
           Trtf = c(0, diff(Ctrt))) %>%
    pivot_longer(names_to = "variable", cols = !c(1,32))%>% 
    mutate(model = "clinical",
           Scenario = scn[i])
  
  # population check preparation
  pop.check <- out.df %>% 
    filter(variable %in% c("P_total"), time > 0)
  
  pop.check.all <- rbind(pop.check.all, pop.check)

  # all cases plot
  all.cases <- out.df  %>% 
    filter(variable %in% c("Cf", "Ctf", "Sev", "Af")) %>% 
    group_by(variable) %>%
    filter(time > 360*4)
  
  all.cases.all <- rbind(all.cases.all, all.cases)
  
  # clinical cases
  clinical.cases <- out.df %>% 
    filter(variable %in% c("Cf", "Ctf", "Sev", "C_total")) %>% 
    group_by(variable) %>%
    filter(time > 360*4)
  
  clinical.cases.all <- rbind(clinical.cases.all, clinical.cases)
  
  # treated cases
  treated.cases <- out.df  %>% 
    filter(variable %in% c("Tf", "H", "T_total")) %>% 
    group_by(variable) %>%
    filter(time > 360*4)
  
  treated.cases.all <- rbind(treated.cases.all, treated.cases)
  
  # deaths
  death <- out.df  %>% 
    filter(variable %in% c("Deathf")) %>% 
    group_by(variable) %>%
    filter(time > 360*4)
  
  death.all <- rbind(death.all, death)
  
  # incidence
  incidence <- out.df %>%
    filter(variable %in% c("Incf")) %>% 
    filter(time > 360*4)
  
  incidence.all <- rbind(incidence.all, incidence)
}

## 2. Disease Dynamics plot ####
# 2.1. population check ####
ggplot(pop.check.all) +
   aes(x = date, y = value) +
   geom_line(colour = "#112446") +
  labs(x = "year", y = "population", title = "Population check") +
   theme_minimal() +
   facet_wrap(vars(Scenario))
  
# 2.2. all cases ####
ggplot(all.cases.all) +
  aes(x = date, y = value, colour = variable) +
  geom_line() +
  scale_color_hue(direction = 1) +
  labs(x = "year", y = "population", title = "All Cases") +
  theme_minimal() +
  facet_wrap(vars(Scenario), scales = "free")
  
# 2.3. clinical cases ####
ggplot(clinical.cases.all) +
  aes(x = date, y = value, colour = variable) +
  geom_line() +
  scale_color_hue(direction = 1) +
  labs(x = "year", y = "population", title = "Clinical Cases") +
  theme_minimal() +
  facet_wrap(vars(Scenario), scales = "free")

# 2.4. treated cases ####
ggplot(treated.cases.all) +
  aes(x = date, y = value, colour = variable) +
  geom_line() +
  scale_color_hue(direction = 1) +
  labs(x = "year", y = "population", title = "Treated Cases") +
  theme_minimal() +
  facet_wrap(vars(Scenario), scales = "free")

# 2.5. Incidence ####
# 2.5.1. All monthly incidence ####
ggplot(incidence.all) +
  aes(x = date, y = value, colour = variable) +
  geom_line() +
  scale_color_hue(direction = 1) +
  labs(x = "year", y = "population", title = "Incidence") +
  theme_minimal() +
  facet_wrap(vars(Scenario), scales = "free")

# 2.5.2. Baseline and reverse monthly incidence ####
baseline.intervention_plot <- incidence.all %>%
  filter(Scenario %in% c("baseline", "reverse 1","reverse 2")) %>% 
  ggplot() +
  aes(x = date, y = value, colour = Scenario) +
  geom_line() +
  scale_color_manual(
    values = c(baseline = "#2A778E",
               `reverse 1` = "#404385",
               `reverse 2` = "#440154")
  ) +
  theme_minimal() +
  labs(title = "Monthly Falciparum Malaria Incidence: Reverse Scenarios vs Baseline", x = "date", y =("population"))

baseline.intervention_plot

# 2.5.3. Baseline and intervention monthly incidence ####
baseline.reverse_plot <- incidence.all %>%
  filter(Scenario %in% c("baseline", "scenario 1", "scenario 2", "scenario 3")) %>% 
  ggplot() +
  aes(x = date, y = value, colour = Scenario) +
  geom_line() +
  scale_color_manual(
    values = c(baseline = "#2A778E",
               `scenario 1` = "#27A882",
               `scenario 2` = "#7BD04F",
               `scenario 3` = "#FF9800")
  ) +
  theme_minimal() +
  labs(title = "Monthly Falciparum Malaria Incidence: Intervention Scenarios vs Baseline", x = "date", y =("population"))

baseline.reverse_plot

# 2.5.4. combine plot ####
grid.arrange(baseline.intervention_plot, baseline.reverse_plot)

# 2.6. deaths ####
ggplot(death.all) +
  aes(x = date, y = value, colour = variable) +
  geom_line() +
  scale_color_hue(direction = 1) +
  labs(x = "year", y = "population", title = "Deaths") +
  theme_minimal() +
  facet_wrap(vars(Scenario), scales = "free")

## 3. Contour Plot for Model fitting ####
# 3.1. contour plot ####
par(mfrow = c(1,2))
contour(NLL.grid, nlevels = 250,x = S0.grid, y = ppn.grid, 
        xlab = "Initial susceptible individuals", ylab = "Person per net", 
        main = "Contour Plot: Grid Search for Maximum Likelihood")
points(x = best.fit.S0, y = best.fit.ppn, pch= 20, col = "red", cex = 2)
points(x = S0.fit, y = ppn.fit, pch= 20, col = "green", cex = 2)

# 3.2. heat plot ####
image(S0.grid, ppn.grid, -log(NLL.grid), col = heat.colors(10000),
      xlab="Initial susceptible individuals",ylab="Person per net", 
      main = "Heat Map: Grid Search for Maximum Likelihood")
points(x = best.fit.S0, y = best.fit.ppn, pch= 20, col = "red", cex = 2)
points(x = S0.fit, y = ppn.fit, pch= 20, col = "green", cex = 2)

# 3.3. daily incidence and data ####
plot_incidence.fit

## 4. Economic Sensitivity Analyses ####
# 4.1 No Behaviour Change ####
par(mfrow = c(2,3))
image(netcoverage.grid, mdacoverage.grid, nbc.nhb, col = hcl.colors(4, palette = "YlOrRd"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Net Health Benefit: No Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -nbc.cost, col = hcl.colors(2, palette = "greens"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Cost: No Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -nbc.nhb.cost.accepted, col = hcl.colors(10, palette = "blues"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "NHB and Cost: No Behaviour Change")

# 4.2. With Behaviour Change ####
image(netcoverage.grid, mdacoverage.grid, bc.nhb, col = hcl.colors(4, palette = "YlOrRd"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Net Health Benefit: With Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -bc.cost, col = hcl.colors(2, palette = "greens"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "Cost: With Behaviour Change")

image(netcoverage.grid, mdacoverage.grid, -bc.nhb.cost.accepted, col = hcl.colors(10, palette = "blues"),
      xlab="Net Coverage",ylab="MDA Coverage", main = "NHB and Cost: With Behaviour Change")
