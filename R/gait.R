# ---- The NHANES gait-speed example ------------------------------------------------

#' The NHANES 2001-2002 gait-speed example
#'
#' Gait speed (m/s, timed 20-ft walk, 6.096 m / walk time) and seven
#' lower-body mobility items of the physical functioning questionnaire, scored
#' 0 = no difficulty, 1 = some, 2 = much, 3 = unable (refused and don't know
#' missing), for persons aged 60+ with a valid walk time and at least one item
#' answered (n = 1,426). Built from the public-use files MSX_B, PFQ_B and
#' DEMO_B by `inst/scripts/gait_data.R` and shipped as [pf_gait]. Analyses
#' are unweighted.
#'
#' Items C1 (walking a quarter mile) and C2 (walking up 10 steps) are not asked
#' of people who need special equipment to walk (PFQ055 = 1, 103 people) and
#' are missing for them. Two sensitivity variants are available:
#' `exclude_aids = TRUE` drops the 161 walking-aid users (needs equipment, a
#' device used in the walk test, or could not walk unassisted);
#' `skip_unable = TRUE` scores the skipped C1 and C2 as 3 = unable.
#'
#' @param exclude_aids drop walking-aid users
#' @param skip_unable score the skipped walking items as unable
#' @param file path to an analysis CSV file with the columns of [pf_gait];
#'   `NULL` (default) uses the package data set [pf_gait]
#' @return list(dat, Y, speed, n, tag, exclude_aids, skip_unable, n_skip_recoded);
#'   `tag` is "", "_noaids" or "_skipunable"
#' @export
gait_data <- function(exclude_aids = FALSE, skip_unable = FALSE,
                      file = NULL) {
  if (exclude_aids && skip_unable) stop("Choose exclude_aids or skip_unable, not both.")
  dat <- if (is.null(file)) caliber::pf_gait else utils::read.csv(file)
  tag <- ""
  if (exclude_aids) { dat <- dat[dat$walking_aid == 0, ]; tag <- "_noaids" }
  n_skip_recoded <- 0
  if (skip_unable) {
    skipped <- dat$needs_equipment == 1 & is.na(dat$C1) & is.na(dat$C2)
    dat$C1[skipped] <- 3; dat$C2[skipped] <- 3
    n_skip_recoded <- sum(skipped); tag <- "_skipunable"
  }
  items <- paste0("C", 1:7)
  list(dat = dat, Y = as.matrix(dat[items]), speed = dat$gait_speed, n = nrow(dat),
       tag = tag, exclude_aids = exclude_aids, skip_unable = skip_unable,
       n_skip_recoded = n_skip_recoded)
}

#' Labels of the gait example
#' @export
gait_item_labels <- c("C1 Walking a quarter mile", "C2 Walking up 10 steps",
                      "C3 Walking between rooms", "C4 Standing up from armless chair",
                      "C5 Getting in and out of bed", "C6 Stooping, crouching, kneeling",
                      "C7 Standing for long periods")

#' @rdname gait_item_labels
#' @export
gait_category_labels <- c("0 No difficulty", "1 Some difficulty", "2 Much difficulty", "3 Unable")
