############################ 05 Blood traits and incident hematological disease ############################
# Cox models of each of 65 baseline blood traits with each of 15 incident hematological outcomes
# in UK Biobank (Tables S20 and S21, Figure 5A-B).
#
# Input : data/ukb_analysis_cohort.rds                one row per participant, made by 01_cohort_and_outcomes.R
#         data/ukb_blood_traits.rds                   baseline blood traits, made by 01_cohort_and_outcomes.R
#         definitions/outcomes.csv                    the 15 outcomes
#         definitions/blood_traits_prospective.csv    the 65 traits (35 blood cell indices, 30 biochemical markers)
# Output: output/Table_S20_blood_cell_traits_incident_disease.xlsx
#         output/Table_S21_biochemical_markers_incident_disease.xlsx
#
# Design: one model per trait and outcome: rank-inverse-normal trait + IBD (time-varying, 6-month lag)
#         + trait-by-IBD interaction + covariates; attained age as the time scale; follow-up starts
#         6 months after baseline; deaths censored;
#         BH-FDR over the traits within each outcome. 

###################################################################################################

library(survival)
library(openxlsx)

data_dir <- "data"
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE)

dat      <- as.data.frame(readRDS(file.path(data_dir, "ukb_analysis_cohort.rds")))
bt       <- as.data.frame(readRDS(file.path(data_dir, "ukb_blood_traits.rds")))
outcomes <- read.csv("definitions/outcomes.csv")
traits   <- read.csv("definitions/blood_traits_prospective.csv")

lag_years  <- 0.5    # time left out after the IBD diagnosis
min_cases  <- 250    # fewer cases than this: model not fitted (about 10 per model term)
min_events <- 5      # fewer cases in IBD than this: within-IBD estimate not estimable
max_se     <- 1.5    # larger SE than this: within-IBD estimate not estimable

# outcomes whose follow-up starts 2 years after baseline (a blood count can precede diagnosis for years)
landmark_2y <- c("leukaemia_composite", "myeloid_leukaemia", "lymphoid_leukaemia", "multiple_myeloma")

covariates <- c("sex", "ethnicity", "education", "deprivation", "deprivation_missing",
                "bmi", "smoking", "alcohol", "diet", "physical_activity")

# z = transformed trait, exposed = IBD person-time
fml <- as.formula(paste("Surv(tstart, tstop, event) ~ z + exposed + z:exposed +",
                        paste(covariates, collapse = " + ")))

# baseline traits (NA for participants not in the trait file, e.g. a hematological diagnosis before baseline)
dat[, traits$trait_id] <- bt[match(dat$f.eid, bt$f.eid), traits$trait_id]
rm(bt)


######################## helpers ########################
# rank-based inverse normal transformation, average ties, NA kept
rint <- function(x) {
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  out[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / sum(ok))
  out
}

# person-time for one outcome: unexposed time before the IBD diagnosis, exposed time from
# 6 months after it; columns tstart/tstop are ages in years, `row` points back to `dat`
make_person_time <- function(dat, outcome) {
  dat$row <- seq_len(nrow(dat))

  # drop people with this outcome on or before baseline
  dat <- dat[dat[[paste0("prior_", outcome)]] %in% FALSE, ]
  ibd <- !is.na(dat$ibd_subtype)

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
  u <- which(u_stop > dat$age_entry)
  unexposed <- data.frame(row = dat$row[u], age_entry = dat$age_entry[u],
                          tstart = dat$age_entry[u], tstop = u_stop[u],
                          event = as.integer(event[u] & age_event[u] <= u_stop[u] &
                                             !(is.finite(gap_start[u]) & age_event[u] >= gap_start[u])),
                          exposed = 0, dat[u, covariates])

  # exposed period: 6 months after diagnosis to the end of follow-up
  e_start <- pmax(dat$age_entry, gap_end)
  e <- which(age_exit > e_start)
  exposed <- data.frame(row = dat$row[e], age_entry = dat$age_entry[e],
                        tstart = e_start[e], tstop = age_exit[e],
                        event = as.integer(event[e] & age_event[e] > e_start[e]),
                        exposed = 1, dat[e, covariates])

  rbind(unexposed, exposed)
}

# one trait (looked up in `dat`) and one outcome: returns one row of results
fit_trait <- function(pt, trait, landmark) {

  # transform over the person-time rows of this outcome, then drop missing values
  z  <- rint(dat[[trait]][pt$row])
  pt <- pt[!is.na(z), ]
  pt$z <- z[!is.na(z)]

  # landmark: follow-up starts `landmark` years after the baseline visit
  lm_start <- pt$age_entry + landmark
  keep <- pt$tstop > lm_start
  pt <- pt[keep, ]
  pt$tstart <- pmax(pt$tstart, lm_start[keep])

  res <- data.frame(participants = length(unique(pt$row)), cases = sum(pt$event),
                    cases_ibd = sum(pt$event[pt$exposed == 1]), person_years = sum(pt$tstop - pt$tstart),
                    HR = NA, LCI = NA, UCI = NA, P = NA, HR_ibd = NA,
                    HR_in_ibd = NA, LCI_in_ibd = NA, UCI_in_ibd = NA, P_interaction = NA, note = "")
  if (res$cases < min_cases) {
    res$note <- sprintf("too few events for the full model (%d)", res$cases)
    return(res)
  }

  # coxph may warn of infinite coefficients (a covariate level without cases, or the nucleated red
  # cell traits within IBD); the SE rule below catches the within-IBD case
  fit <- coxph(fml, data = pt, ties = "efron")
  b <- coef(fit); V <- vcov(fit); s <- summary(fit)$coefficients
  res$HR     <- exp(b["z"])
  res$LCI    <- exp(b["z"] - 1.96 * s["z", "se(coef)"])
  res$UCI    <- exp(b["z"] + 1.96 * s["z", "se(coef)"])
  res$P      <- s["z", "Pr(>|z|)"]
  res$HR_ibd <- exp(b["exposed"])

  # trait HR within IBD person-time = trait main effect + interaction
  b_in  <- b["z"] + b["z:exposed"]
  se_in <- sqrt(V["z", "z"] + V["z:exposed", "z:exposed"] + 2 * V["z", "z:exposed"])
  if (is.finite(se_in) && se_in <= max_se && res$cases_ibd >= min_events) {
    res$HR_in_ibd     <- exp(b_in)
    res$LCI_in_ibd    <- exp(b_in - 1.96 * se_in)
    res$UCI_in_ibd    <- exp(b_in + 1.96 * se_in)
    res$P_interaction <- s["z:exposed", "Pr(>|z|)"]
  } else {
    res$note <- sprintf(paste("trait HR estimable; within-IBD estimate and interaction are not",
                              "(SE %.2f on the log scale, %d exposed events)"), se_in, res$cases_ibd)
  }
  res
}


######################## Cox models: 65 traits x 15 outcomes ########################
results <- data.frame()
for (i in 1:nrow(outcomes)) {
  outcome  <- outcomes$outcome_id[i]
  landmark <- ifelse(outcome %in% landmark_2y, 2, 0.5)
  pt <- make_person_time(dat, outcome)                     # built once, shared by all traits
  for (j in 1:nrow(traits)) {
    res <- fit_trait(pt, traits$trait_id[j], landmark)
    results <- rbind(results, data.frame(trait_id = traits$trait_id[j], outcome_id = outcome,
                                         outcome = outcomes$label[i], landmark = landmark, res))
  }
  print(paste(outcome, "done"))
}
rownames(results) <- NULL

# BH-FDR over the traits within each outcome; the two nucleated red cell traits (99% at the
# detection floor) are outside the family
results <- merge(traits, results, by = "trait_id", sort = FALSE)
results$FDR_P <- NA
for (outcome in unique(results$outcome_id)) {
  k <- results$outcome_id == outcome & results$in_fdr_family & !is.na(results$P)
  results$FDR_P[k] <- p.adjust(results$P[k], method = "BH")
}


######################## Tables S20 and S21 ########################
fmt_p   <- function(p, na = "") ifelse(is.na(p), na, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
# an FDR just below 0.05 gets 4 decimals so that it does not print as 0.050
fmt_fdr <- function(p) ifelse(!is.na(p) & p < 0.05 & sprintf("%.3f", p) == "0.050",
                              sprintf("%.4f", p), fmt_p(p, na = "-"))
fmt_ci  <- function(hr, lci, uci) ifelse(is.na(hr), "", sprintf("%.2f (%.2f, %.2f)", hr, lci, uci))
fmt_n  <- function(x) formatC(x, format = "d", big.mark = ",")

# rows: traits in the order of the definitions file, then outcomes by id
results <- results[order(match(results$trait_id, traits$trait_id), results$outcome_id), ]

make_table <- function(tab) {
  data.frame(Trait = tab$label,
             `UK Biobank field` = ifelse(tab$ukb_field == "", "derived ratio", tab$ukb_field),
             `Trait category` = tab$group, Outcome = tab$outcome,
             Events = fmt_n(tab$cases), `Events in IBD` = fmt_n(tab$cases_ibd),
             Participants = fmt_n(tab$participants),
             `HR per SD (95% CI)` = ifelse(is.na(tab$HR), "not estimable", fmt_ci(tab$HR, tab$LCI, tab$UCI)),
             `P value` = fmt_p(tab$P),
             `FDR-adjusted P (within outcome)` = fmt_fdr(tab$FDR_P),
             `HR per SD in IBD (95% CI)` = fmt_ci(tab$HR_in_ibd, tab$LCI_in_ibd, tab$UCI_in_ibd),
             `P value for trait-by-IBD interaction` = fmt_p(tab$P_interaction),
             `Landmark (years)` = as.character(tab$landmark), Note = tab$note, check.names = FALSE)
}

# FDR_P is missing for traits outside the FDR family (in_fdr_family FALSE) and for models not fitted
estimates <- c("trait_id", "label", "class", "outcome_id", "outcome", "landmark", "participants", "cases",
               "cases_ibd", "person_years", "HR", "LCI", "UCI", "P", "FDR_P", "in_fdr_family", "HR_ibd",
               "HR_in_ibd", "LCI_in_ibd", "UCI_in_ibd", "P_interaction", "note")
files <- c(cell    = "Table_S20_blood_cell_traits_incident_disease.xlsx",
           biochem = "Table_S21_biochemical_markers_incident_disease.xlsx")
for (cl in names(files)) {
  tab <- results[results$class == cl, ]
  wb <- createWorkbook()
  addWorksheet(wb, "table")
  writeData(wb, "table", make_table(tab))
  addWorksheet(wb, "all estimates")                        # numeric values, also behind Figure 5A-B
  writeData(wb, "all estimates", tab[, estimates])
  saveWorkbook(wb, file.path(out_dir, files[cl]), overwrite = TRUE)
}

sessionInfo()
