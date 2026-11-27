# evolvability

`evolvability` is an R-package for calculating evolvability parameters from G-matrices (see Hansen and Houle 2008  <https://doi.org/10.1111/j.1420-9101.2008.01573.x> and Bolstad et al. 2014 <https://doi.org/10.1098/rstb.2013.0255>) and for fitting phylogenetic comparative models that link the rate of evolution of a trait to the state of another evolving trait (see Hansen et al. 2021 Systematic Biology <https://doi.org/10.1093/sysbio/syab079>).

## Installation

```r
# From CRAN
install.packages("evolvability")

# The development version from GitHub
devtools::install_github(repo = "GHBolstad/evolvability", build_vignettes = TRUE)
```

## Example

```r
library(evolvability)

# An additive-genetic variance matrix and 100 random selection gradients
G <- matrix(c(1, 1, 0, 1, 2, 1, 0, 1, 2), ncol = 3) / 10
Beta <- randomBeta(100, 3)

# Evolvability, respondability, conditional evolvability, autonomy and
# integration along the selection gradients
summary(evolvabilityBeta(G, Beta))

# Average evolvability parameters across all directions
evolvabilityMeans(G)
```

## Vignettes

The package includes two vignettes: one on analyzing rates of evolution
(`vignette("Analyzing_rates_of_evolution")`) and one on fitting phylogenetic
mixed models (`vignette("Phylogenetic_mixed_model")`).
