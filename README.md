
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mizerUR

<!-- badges: start -->

<!-- badges: end -->

This is an extension package for the mizer package
<https://sizespectrum.org/mizer/>. Mizer has a size-structured
planktonic resource. With this extension you can also model any number
of unstructured resource components. Such unstructured components are
appropriate whenever the predation on these components is not size
based.

Possible applications include the modelling of detritus as a resource
for detritivores, carrion as a resource for scavengers, or macroflora on
which fish can graze.

## Installation

You can install mizerUR from
[GitHub](https://github.com/gustavdelius/mizerUR) with:

``` r
devtools::install_github("gustavdelius/mizerUR")
```

## Example

The mizerUR package facilitates adding unstructured resource components
to an existing MizerParams object. Here we will illustrate the process
with an example MizerParams object that ships with the mizer package.

I will provide a more realistic example here once it is ready and
tested.

``` r
library(mizer)
library(mizerUR)
params <- NS_params
params@species_params$species
#>  [1] "Sprat"   "Sandeel" "N.pout"  "Herring" "Dab"     "Whiting" "Sole"   
#>  [8] "Gurnard" "Plaice"  "Haddock" "Cod"     "Saithe"
```

Let us introduce a detritus component. The mizerUR package provides a
function `detritus_dynamics()` that we will use to model the detritus.
We set this up by creating a named list

``` r
dynamics <- list(detritus = "detritus_dynamics")
```

The `detritus_dynamics()` function uses two parameters

  - `detritus_external` gives the rate of change of detritus from
    external sources in grams/year,
  - `detritus_proportion` gives the proportion of the biomass consumed
    by fish that flows into the detritus component.

We need to give values for these parameters in a named list

``` r
dynamics_params <- list(`detritus_external` = 0,
                        `detritus_proportion` = 0.1)
```

Finally we need to specify the rate in proportion/year at which fish
encounter detritus biomass. We use an allometric rate for this, which
means that we only need to specify the rate for a fish of 1 gram and
rates at other sizes will be scaled appropriately. This rate is species
specific, so it goes into a new column in the `species_params` data
frame. Let us set some (unrealistic made-up) values.

``` r
params@species_params$rho_detritus <- 
    c(0, 0, 0, 0, 1, 5, 10, 0, 20 ,10, 12, 5) * 1e-15
```

We can now add the detritus component with `setUR()`.

``` r
params <- setUR(params, dynamics, dynamics_params)
```

See the help page for `setUR()` for more details.

The initial value of the resource biomass is `URInitial(params)`. We set
it as follows:

``` r
URInitial(params) <- c(detritus = 10^12)
```

We can now run the model in the usual way using the `project()` function
of the mizer package.

``` r
sim <- project(params, progress_bar = FALSE)
```

We can extract the detritus biomass from the sim object with `UR(sim)`
which will return an array with one row for each saved time in the sim
object. We can use this in a plot for example.

``` r
library(ggplot2)
ggplot(melt(UR(sim))) +
    geom_line(aes(x = t, y = value)) +
    scale_y_log10() +
    xlab("Time [years]") +
    ylab("Detritus biomass [g]")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />
