# =============================================================================
# gait_data.R  --  Analysis data set for the real-data example: timed 20-ft
#                  walk (gait speed) read on seven lower-body mobility items
#
# Source files (NHANES 2001-2002, public use), in $NHANES_DIR (default
# "data/gait"; in the project folder: data/gait/):
#   MSX_B.xpt   muscle strength exam; MSXW20TM = time (s) for the 20-ft walk
#   PFQ_B.xpt   physical functioning questionnaire; PFQ060B..S difficulty items
#   DEMO_B.xpt  demographics; RIDAGEYR = age in years (needed for ages 60+)
#
# Digital / performance measure: gait speed (m/s) = 6.096 m / MSXW20TM.
# PRO: seven lower-body mobility items, 1 = no difficulty, 2 = some, 3 = much,
# 4 = unable, recoded to 0..3. Refused (7) and don't know (9) are missing.
# PFQ060B and PFQ060C are not asked of people who need special equipment to
# walk (PFQ055 = 1); those skips stay missing.
#
# Sample: ages 60+ (everyone is asked the items), a valid walk time and at
# least one of the seven items answered. Walk times under 3 s (over 2 m/s at
# "usual pace") are treated as recording errors.
#
# Without DEMO_B.xpt the script stops, unless GAIT_ALLOW_NO_AGE=1: then it keeps
# everyone who answered the items (ages 60+ plus the 50-59 who reported a
# limitation) and writes gait_analysis_PROVISIONAL.csv instead.
#
# Run:  NHANES_DIR=<folder with the .xpt files> Rscript gait_data.R
# Output in $CALIBER_OUT/data (default output/data): gait_analysis.csv,
# gait_sample_flow.csv, gait_item_counts.csv. Copy gait_analysis.csv to
# data-raw/ and run data-raw/datasets.R to update the package data set
# gait (gait_data()).
# =============================================================================

gait_dir <- Sys.getenv("NHANES_DIR", "data/gait")
out_dir <- Sys.getenv("CALIBER_OUT", "output")
gait_items <- c(walk_quarter_mile = "PFQ060B", climb_10_steps = "PFQ060C",
                walk_between_rooms = "PFQ060H", rise_from_chair = "PFQ060I",
                in_out_of_bed = "PFQ060J", stoop_crouch_kneel = "PFQ060D",
                stand_long_periods = "PFQ060M")

msx <- foreign::read.xport(file.path(gait_dir, "MSX_B.xpt"))
pfq <- foreign::read.xport(file.path(gait_dir, "PFQ_B.xpt"))
have_age <- file.exists(file.path(gait_dir, "DEMO_B.xpt"))
if (!have_age && Sys.getenv("GAIT_ALLOW_NO_AGE", "0") != "1")
  stop("DEMO_B.xpt is missing: age (RIDAGEYR) is needed for the 60+ restriction.")

d <- merge(msx[c("SEQN", "MSXW20TM", "MSAWDEV", "MSDCWUA")],
           pfq[c("SEQN", "PFQ055", unname(gait_items))], by = "SEQN")
flow <- data.frame(step = "In MSX_B and PFQ_B (walk exam, ages 50+)", n = nrow(d))

if (have_age) {
  demo <- foreign::read.xport(file.path(gait_dir, "DEMO_B.xpt"))
  d <- merge(d, demo[c("SEQN", "RIDAGEYR", "RIAGENDR")], by = "SEQN")
  d <- d[d$RIDAGEYR >= 60, ]
  flow <- rbind(flow, data.frame(step = "Aged 60+", n = nrow(d)))
} else {
  d$RIDAGEYR <- NA; d$RIAGENDR <- NA
}

d <- d[!is.na(d$MSXW20TM), ]
flow <- rbind(flow, data.frame(step = "Walk time recorded", n = nrow(d)))
d <- d[d$MSXW20TM >= 3, ]
flow <- rbind(flow, data.frame(step = "Walk time >= 3 s (<= 2.03 m/s)", n = nrow(d)))

Y <- sapply(gait_items, function(v) { x <- d[[v]]; ifelse(x %in% 1:4, x - 1, NA) })
colnames(Y) <- paste0("C", seq_along(gait_items))
d <- d[rowSums(!is.na(Y)) >= 1, ]; Y <- Y[rowSums(!is.na(Y)) >= 1, , drop = FALSE]
flow <- rbind(flow, data.frame(step = if (have_age) "At least one of the 7 items answered"
                               else "At least one of the 7 items answered (60+ or limited 50-59)",
                               n = nrow(d)))

out <- data.frame(SEQN = d$SEQN, age = d$RIDAGEYR, sex = d$RIAGENDR,
                  walk_time = d$MSXW20TM, gait_speed = 6.096 / d$MSXW20TM,
                  # needs special equipment to walk, used a device in the walk test,
                  # or could not walk unassisted
                  walking_aid = as.integer(d$PFQ055 %in% 1 | d$MSAWDEV %in% 2:4 | d$MSDCWUA %in% 2),
                  Y, n_items_answered = rowSums(!is.na(Y)),
                  # needs special equipment to walk (PFQ055 = 1): not asked PFQ060B, C
                  needs_equipment = as.integer(d$PFQ055 %in% 1))
flow <- rbind(flow, data.frame(step = "  of whom walking-aid users", n = sum(out$walking_aid)),
              data.frame(step = "  of whom all 7 items answered", n = sum(out$n_items_answered == 7)))

counts <- do.call(rbind, lapply(seq_along(gait_items), function(j) {
  x <- Y[, j]
  data.frame(item = colnames(Y)[j], variable = gait_items[[j]], content = names(gait_items)[j],
             none = sum(x == 0, na.rm = TRUE), some = sum(x == 1, na.rm = TRUE),
             much = sum(x == 2, na.rm = TRUE), unable = sum(x == 3, na.rm = TRUE),
             missing = sum(is.na(x)))
}))

tag <- if (have_age) "" else "_PROVISIONAL"
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(out_dir, "data", sprintf("gait_analysis%s.csv", tag)), row.names = FALSE)
write.csv(flow, file.path(out_dir, "data", sprintf("gait_sample_flow%s.csv", tag)), row.names = FALSE)
write.csv(counts, file.path(out_dir, "data", sprintf("gait_item_counts%s.csv", tag)), row.names = FALSE)
print(flow, row.names = FALSE); print(counts, row.names = FALSE)
cat(sprintf("Gait speed (m/s): mean %.2f, SD %.2f, range %.2f to %.2f\n", mean(out$gait_speed),
            sd(out$gait_speed), min(out$gait_speed), max(out$gait_speed)))
print(table(walking_aid = out$walking_aid, items_answered = out$n_items_answered))
print(table(needs_equipment = out$needs_equipment,
            C1_C2_missing = paste0(is.na(out$C1) + 0, is.na(out$C2) + 0)))
