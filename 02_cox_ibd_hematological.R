############################ 02 IBD and incident hematological disease ############################
# Cox models for IBD, ulcerative colitis (UC) and Crohn's disease (CD) and 15 hematological
# outcomes in UK Biobank (Table S10, Figure 2), and two sensitivity analyses: immunosuppressant
# use (Table S11) and blood-count abnormalities at baseline (Table S12).
#
# Input : data/ukb_analysis_cohort.rds    one row per participant, made by 01_cohort_and_outcomes.R
#         data/ukb_immunosuppressant.rds  age at first immunosuppressant prescription, made by 01
#         data/ukb_blood_traits.rds       baseline blood counts, made by 01
#         definitions/outcomes.csv        the 15 outcomes (13 primary + 2 composites)
# Output: output/Table_S10_cox_ibd_hematological.xlsx
#         output/Table_S11_immunosuppressant_sensitivity.xlsx
#         output/Table_S12_baseline_blood_count_sensitivity.xlsx
#
# Design: follow-up from baseline to the outcome, death, loss to follow-up or the last date the
#         outcome is recorded; attained age as time scale; IBD time-varying, first 6 months after
#         diagnosis left out; deaths censored; BH-FDR over the 13 primary outcomes per exposure.
#         Sensitivity analyses (S11, S12): the 13 primary outcomes, same model, Wald P values.

###################################################################################################

library(survival)
library(openxlsx)

data_dir <- "data"
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE)

dat      <- as.data.frame(readRDS(file.path(data_dir, "ukb_analysis_cohort.rds")))
outcomes <- read.csv("definitions/outcomes.csv")

lag_years  <- 0.5    # time left out after the IBD diagnosis
min_events <- 5      # fewer exposed cases than this: not estimable
max_se     <- 1.5    # larger standard error than this: not estimable
lr_cutoff  <- 50     # fewer exposed cases than this: P from a likelihood-ratio test

covariates <- c("sex", "ethnicity", "education", "deprivation", "deprivation_missing",
                "bmi", "smoking", "alcohol", "diet", "physical_activity")


######################## person-time for one exposure and one outcome ########################
# Returns one row per participant and period: unexposed time before the IBD diagnosis,
# exposed time from 6 months after it. Columns tstart/tstop are ages in years.
make_person_time <- function(dat, exposure, outcome) {

  # UC and CD are compared with participants without any IBD
  if (exposure == "UC") dat <- subset(dat, is.na(ibd_subtype) | ibd_subtype == "UC")
  if (exposure == "CD") dat <- subset(dat, is.na(ibd_subtype) | ibd_subtype == "CD")
  ibd <- !is.na(dat$ibd_subtype)

  # drop people with this outcome on or before baseline
  keep <- dat[[paste0("prior_", outcome)]] %in% FALSE
  dat <- dat[keep, ]; ibd <- ibd[keep]

  # end of follow-up and event
  age_event <- dat[[paste0("age_", outcome)]]
  age_exit  <- pmin(age_event, dat$age_death, dat$age_ltfu, dat$age_admin_end,
                    dat[[paste0("age_comp_", outcome)]], na.rm = TRUE)
  event <- !is.na(age_event) & age_event <= age_exit

  # the 6 months after diagnosis are not counted as exposed or unexposed
  gap_start <- ifelse(ibd & !is.na(dat$age_at_ibd_dx), dat$age_at_ibd_dx, Inf)
  gap_end   <- gap_start + lag_years

  # unexposed period: baseline to diagnosis (an event in the gap is censored)
  u_stop <- pmin(age_exit, gap_start)
  u <- u_stop > dat$age_entry
  unexposed <- data.frame(id = dat$f.eid[u], tstart = dat$age_entry[u], tstop = u_stop[u],
                          event = as.integer(event[u] & age_event[u] <= u_stop[u] &
                                             !(is.finite(gap_start[u]) & age_event[u] >= gap_start[u])),
                          exposed = 0, dat[u, covariates])

  # exposed period: 6 months after diagnosis to the end of follow-up
  e_start <- pmax(dat$age_entry, gap_end)
  e <- age_exit > e_start
  exposed <- data.frame(id = dat$f.eid[e], tstart = e_start[e], tstop = age_exit[e],
                        event = as.integer(event[e] & age_event[e] > e_start[e]),
                        exposed = 1, dat[e, covariates])

  rbind(unexposed, exposed)
}


######################## Cox models: 3 exposures x 15 outcomes ########################
fml  <- as.formula(paste("Surv(tstart, tstop, event) ~ exposed +", paste(covariates, collapse = " + ")))
fml0 <- as.formula(paste("Surv(tstart, tstop, event) ~",           paste(covariates, collapse = " + ")))

results <- data.frame()
for (exposure in c("IBD", "UC", "CD")) {
  for (i in 1:nrow(outcomes)) {
    outcome <- outcomes$outcome_id[i]
    pt <- make_person_time(dat, exposure, outcome)

    res <- data.frame(exposure = exposure, outcome = outcomes$label[i], group = outcomes$group[i],
                      primary_outcome = outcomes$primary_outcome[i],
                      cases = sum(pt$event), cases_exposed = sum(pt$event[pt$exposed == 1]),
                      participants = length(unique(pt$id)), person_years = sum(pt$tstop - pt$tstart),
                      HR = NA, LCI = NA, UCI = NA, P = NA)

    if (res$cases_exposed >= min_events) {
      fit <- coxph(fml, data = pt, ties = "efron")
      se  <- sqrt(fit$var[1, 1])
      if (is.finite(se) && se <= max_se) {
        ci <- confint(fit)                                   # Wald 95% CI
        res$HR  <- exp(coef(fit)["exposed"])
        res$LCI <- exp(ci["exposed", 1])
        res$UCI <- exp(ci["exposed", 2])
        res$P   <- summary(fit)$coefficients["exposed", "Pr(>|z|)"]
        if (res$cases_exposed < lr_cutoff) {                 # few cases: likelihood-ratio test
          fit0  <- coxph(fml0, data = pt, ties = "efron")
          res$P <- pchisq(2 * (fit$loglik[2] - fit0$loglik[2]), df = 1, lower.tail = FALSE)
        }
      }
    }
    results <- rbind(results, res)
    print(paste(exposure, outcome, "done"))
  }
}

# BH-FDR over the 13 primary outcomes within each exposure (composites are not tested)
results$FDR_P <- NA
for (exposure in c("IBD", "UC", "CD")) {
  k <- results$exposure == exposure & results$primary_outcome & !is.na(results$P)
  results$FDR_P[k] <- p.adjust(results$P[k], method = "BH")
}


######################## Table S10 ########################
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_hr <- function(hr, lci, uci, n = NULL) {
  ne <- if (is.null(n)) "not estimable" else
    sprintf("not estimable (%d exposed event%s)", n, ifelse(n == 1, "", "s"))
  ifelse(is.na(hr), ne, sprintf("%.2f (%.2f, %.2f)", hr, lci, uci))
}

results$`HR (95% CI)`    <- fmt_hr(results$HR, results$LCI, results$UCI, results$cases_exposed)
results$`P value`        <- fmt_p(results$P)
results$`FDR-adjusted P` <- ifelse(results$primary_outcome, fmt_p(results$FDR_P), "-")   # composites: no FDR

wb <- createWorkbook()
for (ex in c("IBD", "UC", "CD")) {
  tab <- results[results$exposure == ex, ]
  tab <- data.frame(Group = tab$group, Outcome = tab$outcome, Cases = tab$cases,
                    `Cases in exposed person-time` = tab$cases_exposed,
                    Participants = tab$participants, `Person-years` = round(tab$person_years),
                    tab[, c("HR (95% CI)", "P value", "FDR-adjusted P")], check.names = FALSE)
  addWorksheet(wb, ex)
  writeData(wb, ex, tab)
}
addWorksheet(wb, "all estimates")
writeData(wb, "all estimates", results[, 1:13])
saveWorkbook(wb, file.path(out_dir, "Table_S10_cox_ibd_hematological.xlsx"), overwrite = TRUE)


######################## sensitivity analyses: helpers ########################
# The Cox model above, plus optional extra terms; Wald P values.
fit_sensitivity <- function(pt, extra = NULL) {
  res <- data.frame(cases = sum(pt$event), cases_exposed = sum(pt$event[pt$exposed == 1]),
                    participants = length(unique(pt$id)), person_years = sum(pt$tstop - pt$tstart),
                    HR = NA, LCI = NA, UCI = NA, P = NA)
  if (res$cases_exposed >= min_events) {
    f   <- as.formula(paste("Surv(tstart, tstop, event) ~",
                            paste(c("exposed", covariates, extra), collapse = " + ")))
    fit <- coxph(f, data = pt, ties = "efron")
    se  <- sqrt(fit$var[1, 1])
    if (is.finite(se) && se <= max_se) {
      ci <- confint(fit)
      res$HR  <- exp(coef(fit)["exposed"])
      res$LCI <- exp(ci["exposed", 1])
      res$UCI <- exp(ci["exposed", 2])
      res$P   <- summary(fit)$coefficients["exposed", "Pr(>|z|)"]
    }
  }
  res
}

# Cuts each period that contains the age at first prescription into an untreated part
# (no event) and a treated part (keeps the event); on_is = 1 from that age onwards.
split_at_age <- function(pt, age) {
  cut   <- !is.na(age) & age > pt$tstart & age < pt$tstop
  lower <- pt[cut, ]; lower$tstop  <- age[cut]; lower$event <- 0L
  upper <- pt[cut, ]; upper$tstart <- age[cut]
  age   <- c(age[!cut], age[cut], age[cut])
  pt    <- rbind(pt[!cut, ], lower, upper)
  pt$on_is <- as.integer(!is.na(age) & pt$tstart >= age)
  pt
}

# BH-FDR over the 13 primary outcomes within each specification and exposure
add_fdr <- function(res) {
  res$FDR_P <- NA
  for (g in unique(paste(res$specification, res$exposure))) {
    k <- paste(res$specification, res$exposure) == g & !is.na(res$P)
    res$FDR_P[k] <- p.adjust(res$P[k], method = "BH")
  }
  res
}

prim <- results[results$primary_outcome, ]                  # primary analysis, for comparison


######################## Table S11: immunosuppressant use ########################
ims     <- as.data.frame(readRDS(file.path(data_dir, "ukb_immunosuppressant.rds")))
is_age  <- setNames(ims$age_first_is, ims$f.eid)            # age at first recorded prescription
pre_ids <- ims$f.eid[ims$is_before_baseline]                # first prescription before baseline

# specification 3: no prescription before baseline; follow-up ends at the first prescription
dat_naive <- subset(dat, !f.eid %in% pre_ids)
dat_naive$age_ltfu <- pmin(dat_naive$age_ltfu, is_age[dat_naive$f.eid], na.rm = TRUE)
cat(nrow(dat_naive), "participants without a prescription before baseline\n")

specs <- c("Specification 1: time-varying immunosuppressant use",
           "Specification 2: adjustment for pre-baseline use",
           "Specification 3: restriction to participants without pre-baseline use")
s11 <- data.frame()
for (s in 1:3) {
  for (exposure in c("IBD", "UC", "CD")) {
    for (i in which(outcomes$primary_outcome)) {
      outcome <- outcomes$outcome_id[i]
      if (s == 1) {                    # use switches on at the first prescription and stays on
        pt  <- make_person_time(dat, exposure, outcome)
        res <- fit_sensitivity(split_at_age(pt, unname(is_age[pt$id])), "on_is")
      }
      if (s == 2) {                    # use before baseline as a fixed covariate
        pt  <- make_person_time(dat, exposure, outcome)
        pt$is_before_baseline <- as.integer(pt$id %in% pre_ids)
        res <- fit_sensitivity(pt, "is_before_baseline")
      }
      if (s == 3) res <- fit_sensitivity(make_person_time(dat_naive, exposure, outcome))
      s11 <- rbind(s11, data.frame(specification = specs[s], exposure = exposure,
                                   outcome = outcomes$label[i], res))
      print(paste("S11", s, exposure, outcome, "done"))
    }
  }
}
s11 <- add_fdr(s11)

wb <- createWorkbook()
for (ex in c("IBD", "UC", "CD")) {
  p   <- prim[prim$exposure == ex, ]
  tab <- data.frame(Outcome = p$outcome,
                    `Primary: events in exposed person-time` = p$cases_exposed,
                    `Primary: HR (95% CI)` = fmt_hr(p$HR, p$LCI, p$UCI, p$cases_exposed),
                    `Primary: P value` = fmt_p(p$P), check.names = FALSE)
  for (s in 1:3) {
    x  <- s11[s11$specification == specs[s] & s11$exposure == ex, ]
    sp <- paste0("Specification ", s, ": ")
    if (s == 3) tab[[paste0(sp, "events in exposed person-time")]] <- x$cases_exposed
    tab[[paste0(sp, "HR (95% CI)")]] <- fmt_hr(x$HR, x$LCI, x$UCI, x$cases_exposed)
    tab[[paste0(sp, "P value")]]     <- fmt_p(x$P)
  }
  addWorksheet(wb, ex)
  writeData(wb, ex, tab)
}
addWorksheet(wb, "all estimates")
writeData(wb, "all estimates", s11)
saveWorkbook(wb, file.path(out_dir, "Table_S11_immunosuppressant_sensitivity.xlsx"), overwrite = TRUE)


######################## Table S12: blood-count abnormalities at baseline ########################
bc <- as.data.frame(readRDS(file.path(data_dir, "ukb_blood_traits.rds")))

# abnormal = measured and outside the clinical range on any of six indices (anemia by the WHO
# sex-specific hemoglobin threshold, g/dL); a count not measured is not abnormal, and participants
# without a blood-count record (left out of the trait file for a hematological diagnosis before
# baseline) are kept
abnormal <- with(bc, (sex == "Male" & haemoglobin < 13) | (sex == "Female" & haemoglobin < 12) |
                     platelet_count < 150 | platelet_count > 450 |   # thrombocytopenia, thrombocytosis
                     wbc_count < 4 | wbc_count > 11 |                # leukopenia, leukocytosis
                     mcv > 100)                                      # macrocytosis
dat_bc <- subset(dat, !f.eid %in% bc$f.eid[abnormal %in% TRUE])
cat(sum(!dat$f.eid %in% bc$f.eid), "participants without a blood-count record kept\n")
cat(nrow(dat) - nrow(dat_bc), "participants excluded,", nrow(dat_bc), "remaining\n")

s12 <- data.frame()
for (exposure in c("IBD", "UC", "CD")) {
  for (i in which(outcomes$primary_outcome)) {
    outcome <- outcomes$outcome_id[i]
    res <- fit_sensitivity(make_person_time(dat_bc, exposure, outcome))
    s12 <- rbind(s12, data.frame(specification = "no baseline abnormality", exposure = exposure,
                                 outcome = outcomes$label[i], res))
    print(paste("S12", exposure, outcome, "done"))
  }
}
s12 <- add_fdr(s12)

wb <- createWorkbook()
for (ex in c("IBD", "UC", "CD")) {
  p   <- prim[prim$exposure == ex, ]
  x   <- s12[s12$exposure == ex, ]
  tab <- data.frame(Outcome = p$outcome,
                    `Primary: events` = p$cases,
                    `Primary: HR (95% CI)` = fmt_hr(p$HR, p$LCI, p$UCI),
                    `Primary: P value` = fmt_p(p$P),
                    `Excluding abnormality: events` = x$cases,
                    `Excluding abnormality: HR (95% CI)` = fmt_hr(x$HR, x$LCI, x$UCI),
                    `Excluding abnormality: P value` = fmt_p(x$P),
                    `Excluding abnormality: FDR-adjusted P` = fmt_p(x$FDR_P), check.names = FALSE)
  addWorksheet(wb, ex)
  writeData(wb, ex, tab)
}
addWorksheet(wb, "all estimates")
writeData(wb, "all estimates", s12)
saveWorkbook(wb, file.path(out_dir, "Table_S12_baseline_blood_count_sensitivity.xlsx"), overwrite = TRUE)

sessionInfo()
