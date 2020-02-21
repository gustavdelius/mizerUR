#' Calculate the biomasses of unstructured resources at next time step
#'
#' This function is called in the `project()` loop at each time step to
#' advance the biomasses of all unstructured resources by one time step.
#' As a user you do not need to deal with this function directly.
#'
#' @param params A [MizerParams] object
#' @param n A matrix of current species abundances (species x size)
#' @param n_pp A vector of current plankton abundance by size
#' @param n_other List of abundances of other dynamic components
#' @param rates A list of rates as returned by [getRates()]
#' @param t Current time
#' @param dt Time step size
#' @param ... Unused
#'
#' @return A named list with the biomasses of all unstructured resources at the
#'   next time step.
#' @md
#' @export
UR_dynamics <- function(params, n, n_pp, n_other, rates, t, dt, ...) {
    UR_dynamics <- params@other_params$UR$dynamics
    UR_dynamics_fns <- lapply(UR_dynamics, get)
    for (res in names(UR_dynamics)) {
        n_other$UR[[res]] <-
            UR_dynamics_fns[[res]](
                params,
                n = n,
                n_pp = n_pp,
                n_other = n_other,
                rates = rates,
                t = t,
                dt = dt
            )
    }
    return(n_other$UR)
}


#' Detritus dynamics
#'
#' Calculates the detritus biomass at the next timestep from the current
#' detritus biomass.
#'
#' The equation for the time evolution of the detritus biomass \eqn{B} is
#' assumed to be of the form
#' \deqn{dB/dt = \tt{inflow} - \tt{consumption} * B + \tt{external}}{dB/dt = inflow - consumption * B + external}
#' where
#' * `inflow` comes from feces, calculated as a proportion
#'   `detritus_proportion` of the biomass consumed by all consumers.
#' * `consumption` is by detritivorous species, where the encounter rate is
#'   specified by `rho[, "detritus", ]`.
#' * `external` is an influx from external sources. It can be negative in which
#'   case it represents a loss to external sources.
#'
#' This equation is solved analytically to
#' \deqn{B(t+dt) = B(t)\exp(-\tt{consumption} \cdot dt)
#'   +\frac{\tt{inflow} + \tt{external}}{\tt{consumption}}
#'   (1-\exp(-\tt{consumption} \cdot dt)).}{B(t+dt)
#'   = B(t) exp(-consumption * dt)
#'   +(inflow + external)/(consumption) * (1 - exp(-consumption * dt)).}
#' This avoids the stability problems that would arise if we used the Euler
#' method to solve the equation numerically.
#'
#' @inheritDotParams UR_dynamics
#' @param detritus_external Rate of change from external sources
#' @param detritus_proportion Proportion of consumption by fish that flows into
#'   the detritus component.
#' @param ... Unused
#'
#' @return A single number giving the biomass of detritus at next time step
#' @export
#' @family resource dynamics functions
#' @md
detritus_dynamics <-
    function(params, n, n_pp, n_other, rates, t, dt,
             detritus_external = params@other_params$UR$detritus_external,
             detritus_proportion = params@other_params$UR$detritus_proportion,
             ...) {

    consumption <- getConsumptionByFish(params, "detritus", n, rates)
    inflow <-
        detritus_proportion *
          sum((rates$feeding_level * params@intake_max * n) %*% params@dw)

    if (consumption) {
        et <- exp(-consumption * dt)
        return(n_other[["UR"]][["detritus"]] * et +
                   (inflow  + detritus_external) / consumption  * (1 - et))
    }
    return(n_other[["UR"]][["detritus"]] + (inflow  + detritus_external) * dt)
}

#' Carrion dynamics
#'
#' Calculates the biomass of carrion (dead animals) at the next timestep from
#' the current biomass.
#'
#' The equation for the time evolution of the carrion biomass \eqn{B} is
#' assumed to be of the form
#' \deqn{dB/dt = inflow - consumption * B + external}
#' where
#' * `inflow` comes from
#'     + Discards from fishing.
#'     + Animals killed by fishing gear.
#'     + Animals that have died by causes other than predation.
#' * `consumption` is by scavenger species, where the encounter rate is
#'   specified by `rho[, "carrion", ]`.
#' * `external` is an influx from external sources. It can be negative in which
#'   case it represents a loss to external sources.
#'
#' This equation is solved analytically to
#' \deqn{B(t+dt) = B(t)\exp(-\tt{consumption} \cdot dt)
#'   +\frac{\tt{inflow} + \tt{external}}{\tt{consumption}}
#'   (1-\exp(-\tt{consumption} \cdot dt)).}{B(t+dt)
#'   = B(t) exp(-consumption * dt)
#'   +(inflow + external)/(consumption) * (1 - exp(-consumption * dt)).}
#' This avoids the stability problems that would arise if we used the Euler
#' method to solve the equation numerically.
#'
#' @inheritParams UR_dynamics
#' @param carrion_external External inflow rate of carrion biomass
#' @param ... Unused
#'
#' @return A single number giving the biomass of carrion at next time step
#' @export
#' @family resource dynamics functions
#' @md
carrion_dynamics <-
    function(params, n, n_pp, n_other, rates, t, dt,
             carrion_external = params@other_params$UR$carrion_external,
             ...) {

        inflow <-
            # still need to be written
            0

        consumption <- getConsumptionByFish(params, "carrion", n, rates)
        if (consumption) {
            et <- exp(-consumption * dt)
            return(n_other[["UR"]][["carrion"]] * et +
                       (inflow  + carrion_external) / consumption  * (1 - et))
        }
        return(n_other[["UR"]][["carrion"]] + (inflow  + carrion_external) * dt)
}

#' Get rate of consumption of resource by fish
#'
#' This function is used by the dynamics functions, for example
#' `detritus_dynamics()` and `carrion_dynamics()`, in their calculation of the
#' resource biomasses at the next time step. Not exported.
#'
#' @param params The MizerParams object
#' @param resource A string with the name of the resource for which the
#'   consumption is to be calculated
#' @param n A matrix with the current size spectrum of fish
#' @param rates A list with the rates
#'
#' @return A number giving the rate at which the resource is consumed by all
#'   the fish in units of grams/year
#' @md
getConsumptionByFish <- function(params, resource, n, rates) {
    sum((params@other_params$UR$rho[, resource, ] * n *
             (1 - rates$feeding_level)) %*% params@dw)
}
