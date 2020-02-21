context("setUR")
library(mizer)
dynamics = list(detritus = "detritus_dynamics", carrion = "carrion_dynamics")
dynamics_params <- list(detritus_external = 0, detritus_proportion = 0,
                        carrion_external = 0)
params <- NS_params
no_sp <- nrow(params@species_params)
params@species_params$rho_detritus <- 1:no_sp
params@species_params$rho_carrion <- no_sp:1
params <- setUR(params, dynamics, dynamics_params)
rho <- params@other_params$UR$rho
sim <- project(params, t_max = 1, dt = 1)

test_that("setUR works", {
    expect_identical(params@initial_n_other,
                   list(UR = c(detritus = 0, carrion = 0)))
    params@initial_n_other <- list(UR = c(detritus = 10, carrion = 20))
    expect_identical(sim@n_other[[2]], list(UR = c(detritus = 0, carrion = 0)))
    expect_identical(params@other_encounter, list(UR = "getUREncounter"))
})

test_that("setUR sets rho", {
    expect_equal(rho[2, 1, 1], 2 * params@w[1]^params@n)
    expect_equal(rho[2, 2, 1], (no_sp - 1) * params@w[1]^params@n)
    rho[2, 1, 1] <- 5
    params <- setUR(params, dynamics, dynamics_params, rho)
    expect_equal(rho[2, 1, 1], 5)
})

test_that("setUR works with single component", {
    dynamics = list(detritus = "detritus_dynamics")
    params <- setUR(params, dynamics, dynamics_params)
    sim <- project(params, t_max = 1, dt = 1)
    expect_identical(sim@n_other[[2]], list(UR = c(detritus = 0)))
})

test_that("setUR throws errors", {
    dynamics_no_names = list("detritus_dynamics", "carrion_dynamics")
    expect_error(setUR(params, dynamics_no_names, dynamics_params),
                 "must be a named list")
    dynamics_typo = list(detritus = "detitus_dynamics",
                         carrion = "carrion_dynamics")
    expect_error(setUR(params, dynamics_typo, dynamics_params),
                 "detitus_dynamics is not defined")
    dynamics_wrong = list(detritus = "sum",
                          carrion = "carrion_dynamics")
    expect_error(setUR(params, dynamics_wrong, dynamics_params),
                 "at least the following arguments")
    expect_error(setUR(params, dynamics),
                 "The following parameters needed")
    dynamics_params_NA <- list(detritus_external = 0, detritus_proportion = NA,
                               carrion_external = 0)
    expect_error(setUR(params, dynamics, dynamics_params_NA),
                 "'detritus_dynamics' are NA ")
    rho_wrong_names <- rho
    dimnames(rho_wrong_names)[[2]][[1]] <- "test"
    expect_error(setUR(params, dynamics, dynamics_params, rho_wrong_names),
                 "not equal to resource_names")
})
