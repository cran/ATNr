## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE, echo=FALSE-----------------------------------------------
oldpar <- par()

## -----------------------------------------------------------------------------
library(ATNr)
set.seed(123)
n_species <- 20 # number of species
conn <- 0.3 # connectance
fw <- create_niche_model(n_species, conn)

## -----------------------------------------------------------------------------
n_species <- 20
n_basal <- 5
masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
L <- create_Lmatrix(masses, n_basal)

## -----------------------------------------------------------------------------
fw <- L
fw[fw > 0] <- 1

## -----------------------------------------------------------------------------
# initialisation of the model object. It is possible to create a ode corresponding to 
# Schneider et al. 2016, Delmas et al. 2016 or Binzer et al. 2016:
# 1) Schneider et al. 2016
n_nutrients <- 3
model_unscaled_nuts <- create_model_Unscaled_nuts(n_species, n_basal, n_nutrients, masses, fw)
# 2) Delmas et al. 2016:
model_scaled <- create_model_Scaled(n_species, n_basal, masses, fw)
# 3) Binzer et al. 2016
model_unscaled <- create_model_Unscaled(n_species, n_basal, masses, fw)

## -----------------------------------------------------------------------------
# updating the hill coefficient in the Unscaled_nuts model:
model_unscaled_nuts$q <- 1.4
# Changing the assimilation efficiencies of all species to 0.5 in the Scaled model:
model_scaled$e = rep(0.5, model_scaled$nb_s)
# print the different fields that can be updated and their values:
str(model_unscaled_nuts)

## -----------------------------------------------------------------------------
# for a model created by create_model_Unscaled_nuts():
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, L)
model_unscaled_nuts$initialisations()
# for a model created by create_model_Scaled():
model_scaled <- initialise_default_Scaled(model_scaled)
model_scaled$initialisations()
# for a model created by create_model_Unscaled():
model_unscaled <- initialise_default_Unscaled(model_unscaled)
model_unscaled$initialisations()

## ----wrappers-----------------------------------------------------------------
biomasses <- masses ^ -0.75 * 1e1 # starting biomasses
biomasses <- append(runif(3, 20, 30), biomasses) # nutrient concentration
# defining the desired integration time
times <- seq(0, 1500, 5)
sol <- lsoda_wrapper(times, biomasses, model_unscaled_nuts)

## ----deSolve------------------------------------------------------------------
# running simulations for the Schneider model
sol <- deSolve::lsoda(
  biomasses,
  times,
  function(t, y, params) {
    return(list(params$ODE(y, t)))
  },
  model_unscaled_nuts
)

## ----plot_odeweb, fig.width=4, fig.height=3, fig.align='center'---------------
par(mar = c(4, 4, 1, 1))
plot_odeweb(sol, model_unscaled_nuts$nb_s)

## -----------------------------------------------------------------------------
# function to plot the fw
show_fw <- function(mat, title = NULL) {
  par(mar = c(.5, .5, 2, .5))
  S <- nrow(mat)
  mat <- mat[nrow(mat):1, ]
  mat <- t(mat)
  image(mat, col = c("goldenrod", "steelblue"),
        frame = FALSE, axes = FALSE)
  title(title)
  grid(nx = S, ny = S, lty = 1, col = adjustcolor("grey20", alpha.f = .1))
}

## -----------------------------------------------------------------------------
S <- 50 # number of species
C <- 0.2 # connectance
fw <- create_niche_model(S, C)

## -----------------------------------------------------------------------------
# number of species and body masses
n_species <- 20
n_basal <- 5
# body mass of species. Here we assume two specific rules for basal and non basal species
masses <- c(sort(10^runif(n_basal, 1, 3)), sort(10^runif(n_species - n_basal, 2, 6)))
L <- create_Lmatrix(masses, n_basal, Ropt = 100, gamma = 2, th = 0.01)

## ---- fig.width=4, fig.height=4, fig.align='center'---------------------------
# boolean version
fw <- L > 0
# 0/1 version:
fw <- L
fw[fw > 0] <- 1
show_fw(fw, title = "L-matrix model food web")

## -----------------------------------------------------------------------------
set.seed(12)
# 1) define number of species, their body masses, and the structure of the
# community
n_species <- 50
n_basal <- 20
n_nut <- 2
# body mass of species
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
# 2) create the food web
# create the L matrix
L <- create_Lmatrix(masses, n_basal, Ropt = 50, gamma = 2, th = 0.01)
# create the 0/1 version of the food web
fw <- L
fw[fw > 0] <- 1
# 3) create the model
model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
# 4) define the temperature gradient and initial conditions
temperatures <- seq(4, 22, by = 2)
extinctions <- rep(NA, length(temperatures))
# defining biomasses
biomasses <- runif(n_species + n_nut, 2, 3)
# 5) define the desired integration time.
times <- seq(0, 100000, 100)
# 6) and loop over temperature to run the population dynamics
i <- 0
for (t in temperatures){
  # initialise the model parameters for the specific temperature
  # Here, no key parameters (numbers of species or species' body masses) are modified
  # Therefore, it is not needed to create a new model object
  # TO reinitialise the different parameters is enough
  model <- initialise_default_Unscaled_nuts(model, L, temperature = t)
  model$q <- 1.2
  model$S <- rep(10, n_nut)
  model$initialisations()
  # running simulations for the Schneider model:
  sol <- lsoda_wrapper(times, biomasses, model, verbose = FALSE)
  # retrieve the number of species that went extinct before the end of the
  # simulation excluding here the 3 first columns: first is time, 2nd and 3rd
  # are nutrients
  i <- i + 1
  extinctions[i] <- sum(sol[nrow(sol), 4:ncol(sol)] < 1e-6)
}

## ---- fig.width=4, fig.height=3, fig.align='center'---------------------------
plot(temperatures, extinctions,
     pch = 20, cex = 0.5, ylim = c(0,50), frame = FALSE,
     ylab = "Number of Extinctions", xlab = "Temperature (°C)")
lines(temperatures, extinctions, col = 'blue')

## ----binzer example-----------------------------------------------------------
# set.seed(142)

# number of species
S <- 30 

# vector containing the predator prey body mass ratios to test
scaling <- 10 ^ seq(-1, 4, by = .5)

# vectors to store the results
persistence0 <- c()
persistence40 <- c()

# create the studied food web
fw <- create_niche_model(S = S, C = 0.1)
# calculating trophic levels
TL = TroLev(fw)
biomasses <- runif(S, 2, 3)

# run a loop over the different pred-prey body mass ratios
for (scal in scaling) {
  # update species body masses following the specific body mass ratio
  masses <- 0.01 * scal ^ (TL - 1)
  
  # create the models with parameters corresponding to 0 and 40 degrees Celcius
  mod0 <- create_model_Unscaled(nb_s = S,
                              nb_b = sum(colSums(fw) == 0),
                              BM = masses,
                              fw = fw)
  mod0 <- initialise_default_Unscaled(mod0, temperature = 0)
  mod0$c <- rep(0, mod0$nb_s - mod0$nb_b)
  mod0$alpha <- diag(mod0$nb_b)
  mod0$initialisations()
  
  mod40 <- create_model_Unscaled(nb_s = S,
                               nb_b = sum(colSums(fw) == 0),
                               BM = masses,
                               fw = fw)
  mod40 <- initialise_default_Unscaled(mod40, temperature = 40)
  mod40$c <- rep(0, mod40$nb_s - mod40$nb_b)
  mod40$alpha <- diag(mod40$nb_b)
  mod40$initialisations()
  
  times <- seq(1, 1e9, by = 1e7)
  
  # run the model corresponding to the 0 degree conditions
  sol <- lsoda_wrapper(times, biomasses, mod0, verbose = FALSE)
  persistence0 <- append(persistence0, sum(sol[nrow(sol), -1] > mod0$ext) / S)
  # run the model corresponding to the 40 degrees conditions
  sol <- lsoda_wrapper(times, biomasses, mod40, verbose = FALSE)
  persistence40 <- append(persistence40, sum(sol[nrow(sol), -1] > mod40$ext) / S)
}

## ----binzer example plot, fig.width=6, fig.height=4, fig.align='center'-------
plot(log10(scaling), persistence40,
     xlab = expression("Body mass ratio between TL"[i + 1]* " and TL"[i]),
     ylab = "Persistence",
     ylim = c(0, 1),
     frame = FALSE, axes = FALSE, type = 'l', col = "red")
lines(log10(scaling), persistence0, col = "blue")
axis(2, at = seq(0, 1, by = .1), labels = seq(0, 1, by = .1))
axis(1, at = seq(-1, 4, by = 1), labels = 10 ^ seq(-1, 4, by = 1))
legend(0.1, 0.9, legend = c("40 \u00B0C", "0 \u00B0C"), fill = c("red", "blue"))

## ----delmas 1-----------------------------------------------------------------
set.seed(1234)
S <- 10
fw <- NULL
TL <- NULL
fw <- create_niche_model(S, C = .15)
TL <- TroLev(fw)

masses <- 0.01 * 100 ^ (TL - 1) #body mass of species

mod <- create_model_Scaled(nb_s = S, 
                           nb_b = sum(colSums(fw) == 0),
                           BM = masses,
                           fw = fw)
mod <- initialise_default_Scaled(mod)
mod$initialisations()
times <- seq(0, 500, by = 2)
biomasses <- runif(S, 2, 3) # starting biomasses

## ----delmas 2-----------------------------------------------------------------
mod$K <- 1
sol1 <- lsoda_wrapper(times, biomasses, mod, verbose = FALSE)
mod$K <- 10
sol10 <- lsoda_wrapper(times, biomasses, mod, verbose = FALSE)

## ----delmas 3, fig.width=6, fig.height=6, fig.align='center'------------------
par(mfrow = c(2, 1))
plot_odeweb(sol1, S)
title("Carrying capacity = 1")
plot_odeweb(sol10, S)
title("Carrying capacity = 10")

## ----mistake 1----------------------------------------------------------------
set.seed(1234)
nb_s <- 20
nb_b <- 5
nb_n <- 2
masses <- sort(10 ^ runif(nb_s, 2, 6)) #body mass of species
biomasses = runif(nb_s + nb_n, 2, 3)
L <- create_Lmatrix(masses, nb_b, Ropt = 50)
L[, 1:nb_b] <- 0
fw <- L
fw[fw > 0] <- 1
model_unscaled_nuts <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, masses, fw)
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, L)
nb_s <- 30 #this does not change the model parameter
model_unscaled_nuts$nb_s #this is the model parameter

## -----------------------------------------------------------------------------
times <- seq(0, 15000, 150)
model_unscaled_nuts$nb_s = 40
model_unscaled_nuts$initialisations()
# this will return an error :
# sol <- lsoda_wrapper(times, biomasses, model_schneider)

## ----mistake 2, fig.width=6---------------------------------------------------
set.seed(1234)
nb_s <- 20
nb_b <- 5
nb_n <- 2
masses <- sort(10 ^ runif(nb_s, 2, 6)) #body mass of species
biomasses = runif(nb_s + nb_n, 2, 3)
L <- create_Lmatrix(masses, nb_b, Ropt = 50)
L[, 1:nb_b] <- 0
fw <- L
fw[fw > 0] <- 1
model_unscaled_nuts <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, masses, fw)
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, L)
model_unscaled_nuts$BM <- sqrt(model_unscaled_nuts$BM) # we change body masses within the model
model_unscaled_nuts$initialisations()
sol <- lsoda_wrapper(seq(1, 5000, 50), biomasses, model_unscaled_nuts)
par(mar = c(4, 4, 1, 1))
plot_odeweb(sol, model_unscaled_nuts$nb_s)

## -----------------------------------------------------------------------------
nb_s <- 30
nb_n <- 2
masses <- sort(10 ^ runif(nb_s, 2, 6)) #body mass of species
biomasses <- runif(nb_s + nb_n, 2, 3)
L <- create_Lmatrix(masses, nb_b, Ropt = 50)
L[, 1:nb_b] <- 0
fw <- L
fw[fw > 0] <- 1
# create a new object:
model_unscaled_nuts <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, masses, fw)
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, L)
# safely run the integration:
model_unscaled_nuts$initialisations()
sol <- lsoda_wrapper(times, biomasses, model_unscaled_nuts)

## ----mistake 3----------------------------------------------------------------
nb_s <- 30
masses <- sort(10 ^ runif(nb_s, 2, 6)) #body mass of species
L <- create_Lmatrix(masses, nb_b, Ropt = 50)
L[, 1:nb_b] <- 0
fw <- L
fw[fw > 0] <- 1
# create a new object:
model_unscaled_nuts <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, masses, fw)
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, L)
model_unscaled_nuts$initialisations() # commenting this line will lead to wrong results.
sol <- lsoda_wrapper(times, biomasses, model_unscaled_nuts)

## ----mistake 4----------------------------------------------------------------
nb_s <- 30
nb_n <- 2
masses <- sort(10 ^ runif(nb_s, 2, 6)) #body mass of species
biomasses <- runif(nb_s + nb_n, 2, 3)
L <- create_Lmatrix(masses, nb_b, Ropt = 50)
L[, 1:nb_b] <- 0
fw <- L
fw[fw > 0] <- 1
# create a new object:
model_1 <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, masses, fw)
model_1 <- initialise_default_Unscaled_nuts(model_1, L)

# trying to create a new model that is similar to model_1
model_2 = model_1

## -----------------------------------------------------------------------------
model_1$q = 1.8
# this also updated the value in model_2:
model_2$q

## ----restore par, include=FALSE, echo=FALSE-----------------------------------
par(oldpar)

