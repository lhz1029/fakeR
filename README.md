# fakeR
An R package to generate fake data given an existing dataset. To download the package,  

## Description 
Generates fake data from a dataset of different variable types.
The package contains the functions simulate_dataset and simulate_dataset_ts 
to simulate time-independent and time-dependent data. It randomly samples 
character and factor variables from contingency tables and  numeric and 
ordered factors from a multivariate normal distribution. It currently supports the 
simulation of stationary and zero-inflated count time series. 
## License
CC0
## Imports
mvtnorm,
polycor,
VGAM,
stats
## Suggests
knitr,
rmarkdown,
testthat
# VignetteBuilder
knitr
