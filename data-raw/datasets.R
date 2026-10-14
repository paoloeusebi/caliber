# datasets.R  --  Builds the package data sets in data/ from the CSV files in
# data-raw/. Run from the package root:  Rscript data-raw/datasets.R
#
# psqi_sleep: sim_psqi_actigraphy.csv joined with eta_true from
#             sim_eta_true.csv by id (simulate_sleep_data(seed = 20261003)).
# pf_gait:    gait_analysis.csv, written by inst/scripts/gait_data.R.

psqi_sleep <- utils::read.csv("data-raw/sim_psqi_actigraphy.csv")
eta <- utils::read.csv("data-raw/sim_eta_true.csv")
stopifnot(nrow(psqi_sleep) == 1000, setequal(psqi_sleep$id, eta$id))
psqi_sleep$eta_true <- eta$eta_true[match(psqi_sleep$id, eta$id)]

pf_gait <- utils::read.csv("data-raw/gait_analysis.csv")
stopifnot(nrow(pf_gait) == 1426)

usethis::use_data(psqi_sleep, pf_gait, overwrite = TRUE, compress = "xz")
