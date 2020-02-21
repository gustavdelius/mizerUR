#' Get encounter rate for unstructured resource
#'
#' Calculates the rate \eqn{E_i(w)} at which a predator of species \eqn{i} and
#' weight \eqn{w} encounters resources (grams/year).
#'
#' @section Resource encounter:
#' In addition to the contribution from predation on fish prey and plankton,
#' the food encounter rate may have a contribution from unstructured resource
#' components. This takes the form
#' \deqn{E_{u.i} = \sum_d \rho_{id}(w) B_d.}
#' where \eqn{B_d} is the biomass of the d-th unstructured resource component
#' and \eqn{\rho_{id}(w)} is a parameter that therefore determines the rate at
#' which a predator of species \eqn{i} and size \eqn{w} encounters biomass from
#' the d-th unstructured resource component. This is set with
#' \code{\link{setResourceEncounter}}.
#'
#' The encounter rate is multiplied by \eqn{1-f_0} to obtain the consumption rate,
#' where \eqn{f_0} is the feeding level calculated with \code{\link{getFeedingLevel}}.
#' This is used by the \code{\link{project}} function for performing simulations.
#'
#' The function returns values also for sizes outside the size-range of the
#' species. These values should not be used, as they are meaningless.
#'
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param n_other A list of abundances of other dynamical components
#'
#' @return A two dimensional array (predator species x predator size)
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' # TODO
#' }
getUREncounter <- function(params, n = params@initial_n,
                           n_pp = params@initial_n_pp,
                           n_other = params@initial_n_other,
                           ...) {

    # Create a vector of the right dimension and with the right dimnames
    encounter <- params@metab
    encounter[] <- 0
    # Add contributions from unstructured resources
    # Can't use rowSums or colSums unfortunately because
    # the resource index that we want to sum over is the middle index.
    for (u in seq_along(n_other$UR)) {
        encounter[] <- encounter +
            params@other_params$UR$rho[, u, ] * n_other$UR[[u]]
    }
    return(encounter)
}
