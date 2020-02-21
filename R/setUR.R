#' Set up unstructured resource components
#'
#' Allows you to add any number of unstructured resource components to an
#' existing MizerParams object. Such unstructured components are appropriate
#' whenever the predation on these components is not size based. Examples
#' include detritus as a resource for detritivores, carrion as a resource for
#' scavengers, or macroflora on which fish can graze.
#'
#' @section Setting resource dynamics:
#'
#' During a simulation using [project()], the biomasses of the resources are
#' updated at each time step by calling the functions specified in the
#' `dynamics` list. This list should have one named entry for each unstructured
#' resource component, giving the name of the function as a string.
#'
#' Mizer provides two example functions that you can use to model resource
#' dynamics: [detritus_dynamics()] and [carrion_dynamics()], but you can easily
#' implement others by following those templates.
#'
#' As you can see in the documentation of these functions, their arguments are:
#' the `MizerParams` object `params`, the current fish size spectra `n`, the
#' current plankton spectrum `n_pp`, the abundances of any other components
#' 'n_other', the current rates calculated by the [getRates()] function `rates`,
#' the current time `t` and the size of the time step `dt`.
#'
#' The other arguments to the resource dynamics functions are model parameters,
#' like for example growth rates. These need to be provided in the
#' `dynamics_params` argument which is a named list. One model parameter that
#' should always be present in this list is the rate of change due to external
#' causes. This should be given a name of the form `resource_external` where
#' `resource` should be replaced by the name of the resource, see for example
#' `detritus_external` in [detritus_dynamics()].
#'
#' When writing your own resource dynamics functions, you can choose any names
#' for your other model parameters, but you must make sure not to use the same
#' name in the function for another resource component. One way to ensure this
#' is to prefix all parameter names with your resource name.
#'
#' The dynamics for a resource should always have a loss term accounting for
#' the consumption of the resource. The contribution to this loss arising from
#' consumption by fish should be calculated with
#' \code{\link{getConsumptionByFish}}.
#'
#' @section Setting resource encounter rate:
#' The resource encounter rate \eqn{\rho_{id}(w)} (units 1/year) determines the
#' rate at which an individual of species \eqn{i} encounters biomass of resource
#' \eqn{d}, so that the contribution from all unstructured resources to the
#' total encounter rate is
#' \deqn{E_{u.i}(w) = \sum_d\rho_{id}(w) B_d,}
#' where \eqn{B_d} is the biomass of the d-th unstructured resource component.
#'
#' Resource consumption is subject to satiation in the same way as other food,
#' so that a consumer only consumes a fraction \eqn{1-f_i(w)} of the encountered
#' resource biomass, where \eqn{f_i(w)} is the feeding level.
#'
#' If the \code{rho} array is not supplied, then the resource encounter rate is
#' set to a power law
#' \deqn{\rho_{id}(w) = \rho_{id} w^n.}
#' The coefficients \eqn{\rho_{id}} are parameters in the
#' \code{params@species_params} dataframe. For example if there is a resource
#' called "detritus" then the species_params data frame needs to have a column
#' called \code{rho_detritus} and similarly for each other resource.
#'
#' If the \code{rho} array is supplied, the ordering of the entries in the array
#' is important. The order of the species in the first array dimension needs to
#' be the same as that in the species parameter dataframe. The order of the
#' resources in the second array dimension must be the same as in the list of
#' resource dynamics. The third dimension is the size dimension.
#'
#' @param params A MizerParams object
#' @param dynamics A named list of functions that determine the
#'   dynamics of the unstructured resources by calculating their biomasses at
#'   the next time step from the current state. Details are described in the
#'   section "Setting resource dynamics".
#' @param dynamics_params A named list of parameters needed by the
#'   \code{dynamics} functions. An empty list if no parameters are
#'   needed.
#' @param rho Optional. An array (species x resource x size)
#'   holding the rate at which a fish of a particular species and of a
#'   particular size feeds on each resource. Described in the section "Setting
#'   resource encounter rate".
#'
#' @return A MizerParams object
#' @export
#' @md
#' @family functions for setting parameters
setUR <- function(params, dynamics,
                          dynamics_params = list(),
                          rho = NULL) {
    # Check arguments ----
    assert_that(is(params, "MizerParams"),
                is.list(dynamics),
                length(dynamics) > 0,
                is.list(dynamics_params))
    no_sp <- nrow(params@species_params)
    no_res <- length(dynamics)
    resource_names = names(dynamics)
    if (is.null(resource_names)) {
        stop("The `dynamics` list must be a named list.")
    }
    # Check that dynamics functions have the correct arguments and that
    # all necessary dynamic parameters are supplied
    for (d in dynamics) {
        fun <- get0(d)
        if (is.null(fun) || !is.function(fun)) {
            stop("The dynamics function ", d, " is not defined.")
        }
        arg <- names(formals(d))
        required <- c("params", "n", "n_pp", "n_other", "rates", "t", "dt")
        if (!identical(arg[1:7], required)) {
            stop("The dynamics function '", d, "' needs to have at least the ",
                 "following arguments: ", required)
        }
        # lop off the required arguments
        arg <- arg[!(arg %in% c(required, "..."))]
        missing <- !(arg %in% names(dynamics_params))
        if (any(missing)) {
            stop("The following parameters needed for the dynamics function '",
                 d, "' are missing in `dynamics_params`: ",
                 arg[missing])
        }
        # Check that there are no missing values
        if (anyNA(dynamics_params[arg])) {
            stop("Some parameters for the dynamics function '", d,
                 "' are NA in `dynamics_params`.")
        }
    }

    # Default for rho ----
    if (is.null(rho)) {
        rho <-
            array(NA,
                  dim = c(no_sp, no_res, length(params@w)),
                  dimnames = list(sp = params@species_params$species,
                                  res = resource_names,
                                  w = signif(params@w, 3)))
        # Use columns in species_params
        for (res in resource_names) {
            rho_var <- paste0("rho_", res)
            if (!(rho_var %in% names(params@species_params))) {
                stop("The species_params data frame needs a column ", rho_var)
            }
            rho[, res, ] <-
                outer(params@species_params[[rho_var]], params@w^params@n)
        }
    }

    # Check rho ----
    assert_that(is.array(rho),
                length(dim(rho)) == 3)
    if (no_sp != dim(rho)[[1]]) {
        stop("The first dimension of the rho argument should equal the number ",
             "of species.")
    }
    if (no_res != dim(rho)[[2]]) {
        stop("The second dimension of the rho argument should equal the number ",
             "of resources.")
    }
    if (is.character(dimnames(rho)[["res"]])) {
        assert_that(are_equal(dimnames(rho)[["res"]], resource_names))
    }
    if (length(params@w) != dim(rho)[[3]]) {
        stop("The third dimension of the rho array should have one entry for ",
             "every consumer size.")
    }
    assert_that(!anyNA(rho),
                all(rho >= 0))

    # Set slots ----
    params@other_dynamics$UR <- "UR_dynamics"
    params@other_encounter$UR <- "getUREncounter"
    params@other_params$UR <- dynamics_params
    params@other_params$UR$dynamics <- dynamics
    params@other_params$UR$rho <- rho

    params@initial_n_other$UR <- rep(0, no_res)
    names(params@initial_n_other$UR) <- resource_names

    # Set linecolour and linetype if necessary ----
    # Colour-blind-friendly palette
    # From http://dr-k-lo.blogspot.co.uk/2013/07/a-color-blind-friendly-palette-for-r.html
    cbbPalette <- c("#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00",
                    "#CC79A7", "#F0E442")
    unused <- setdiff(cbbPalette, params@linecolour)
    colourless <- setdiff(resource_names, names(params@linecolour))
    linecolour <- rep(unused, length.out = length(colourless))
    names(linecolour) <- colourless
    params@linecolour <- c(params@linecolour, linecolour)

    linetype <- rep("solid", length.out = length(colourless))
    names(linetype) <- colourless
    params@linetype <- c(params@linetype, linetype)

    return(params)
}
