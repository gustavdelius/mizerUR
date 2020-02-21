context("getUREncounter")
library(mizer)
dynamics = list(detritus = "detritus_dynamics", carrion = "carrion_dynamics")
dynamics_params <- list(detritus_external = 0, detritus_proportion = 0,
                        carrion_external = 0)
params <- NS_params
params@species_params$rho_detritus <- 1
params@species_params$rho_carrion <- 2
params <- setUR(params, dynamics, dynamics_params)

test_that("getUREncounter returns zero when given zero n_other", {
    encounter <- getUREncounter(params)
    zero <- params@metab
    zero[] <- 0
    expect_identical(encounter, zero)
})
