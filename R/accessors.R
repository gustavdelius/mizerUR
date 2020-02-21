#' Access initial values for unstructured resources
#'
#' @param params The MizerParams object
#' @return Named vector with initial biomasses of unstructured resources
#' @export
URInitial <- function(params) params@initial_n_other[["UR"]]

#' Setter for initial values for unstructured resources
#'
#' @param params The MizerParams object
#' @param value Vector with initial biomasses of unstructured resources.
#' @return MizerParams object
#' @export
`URInitial<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                is.numeric(value),
                all(value >= 0),
                !anyNA(value),
                length(value) == length(params@initial_n_other[["UR"]]))
    if (is.null(names(value))) {
       names(value) <- names(params@initial_n_other[["UR"]])
    }
    if (names(value) != names(params@initial_n_other[["UR"]])) {
        stop("The `value` argument must be a named list of values with ",
             "the names: ", names(params@initial_n_other[["UR"]]))
    }
    params@initial_n_other[["UR"]] <- value
    validObject(params)
    params
}

#' Access simulation results for unstructured resources
#'
#' @param sim The MizerSim object
#' @return A named matrix with rows corresponding to the times in the simulation
#'   and columns corresponging to the resources.
#' @export
UR <- function(sim) {
    matrix(unlist(lapply(sim@n_other, function(x) x[["UR"]])),
           nrow = length(sim@n_other), byrow = TRUE,
           dimnames = list(t = names(sim@n_other),
                           res = names(sim@n_other[[1]][["UR"]])))
}
