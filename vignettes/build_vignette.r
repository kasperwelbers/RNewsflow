## Run each time after building

devtools::install(build_vignettes = TRUE)
devtools::build_vignettes()

vignette('RNewsflow')
