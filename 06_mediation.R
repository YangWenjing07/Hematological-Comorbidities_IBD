############################ 06 Mediation by baseline blood traits ############################
# Single-trait mediation of the associations between IBD and incident hematological disease by
# 65 baseline blood traits in UK Biobank (Table S22, Figure 5C).
#
# Input : data/ukb_analysis_cohort.rds   one row per participant, made by 01_cohort_and_outcomes.R
#         data/ukb_blood_traits.rds      baseline blood traits, made by 01_cohort_and_outcomes.R
#         output/Table_S10_cox_ibd_hematological.xlsx   made by 02 (outcomes with FDR P < 0.05)
#         definitions/outcomes.csv, mediation_traits.csv, mediation_definitional_pairs.csv
# Output: output/Table_S22_mediation_blood_traits.xlsx   (sheet "all pairs": every estimate)
#
# Design: prevalent IBD (on or before baseline) against never IBD; incident IBD and anyone with a
#         hematological diagnosis other than VTE at baseline left out. lm of the rank-inverse-normal
#         trait on IBD; Cox model on attained age with trait, IBD and trait x IBD, from 6 months after
#         baseline. Natural effects (VanderWeele) from 20,000 Monte Carlo draws; BH-FDR without nRBC.

###################################################################################################

library(survival)
library(MASS)
library(openxlsx)
library(parallel)

data_dir <- "data"
out_dir  <- "output"

dat      <- as.data.frame(readRDS(file.path(data_dir, "ukb_analysis_cohort.rds")))
blood    <- as.data.frame(readRDS(file.path(data_dir, "ukb_blood_traits.rds")))
s10      <- read.xlsx(file.path(out_dir, "Table_S10_cox_ibd_hematological.xlsx"),
                      sheet = "all estimates")
outcomes <- read.csv("definitions/outcomes.csv")
traits   <- read.csv("definitions/mediation_traits.csv")
def_pair <- read.csv("definitions/mediation_definitional_pairs.csv")

lag_years  <- 0.5        # exposed time starts 6 months after the IBD diagnosis
min_events <- 250        # fewer outcome events than this: not estimable
n_draws    <- 20000      # Monte Carlo draws
seed       <- 20260906   # reset before every pair
n_cores    <- 1          # forked workers (not on Windows); same results, the seed is per pair
clonal     <- c("leukaemia_composite", "myeloid_leukaemia", "lymphoid_leukaemia", "multiple_myeloma")

covariates <- c("sex", "ethnicity", "education", "deprivation", "deprivation_missing",
                "bmi", "smoking", "alcohol", "diet", "physical_activity")
fml_m <- as.formula(paste("z ~ exposed +", paste(covariates, collapse = " + ")))
fml_o <- as.formula(paste("Surv(tstart, tstop, event) ~ z + exposed + z:exposed +",
                          paste(covariates, collapse = " + ")))

# outcomes whose incidence IBD raised at FDR P < 0.05 (Table S10)
sig <- subset(s10, exposure == "IBD" & !is.na(FDR_P) & FDR_P < 0.05)
sig$outcome_id <- outcomes$outcome_id[match(sig$outcome, outcomes$label)]


######################## helpers ########################
# rank-inverse-normal transform, ties averaged, NA kept
rint <- function(x) {
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  out[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
  out
}

# One row per participant: never IBD from baseline, prevalent IBD from 6 months after diagnosis.
make_person_time <- function(dat, outcome) {
  dat <- subset(dat, exposure_group != "incident")           # IBD must precede the blood draw
  dat <- dat[dat[[paste0("prior_", outcome)]] %in% FALSE, ]  # outcome on or before baseline
  ibd <- dat$exposure_group == "prevalent"

  age_event <- dat[[paste0("age_", outcome)]]
  age_exit  <- pmin(age_event, dat$age_death, dat$age_ltfu, dat$age_admin_end,
                    dat[[paste0("age_comp_", outcome)]], na.rm = TRUE)
  tstart <- ifelse(ibd, pmax(dat$age_entry, dat$age_at_ibd_dx + lag_years), dat$age_entry)
  event  <- !is.na(age_event) & age_event <= age_exit & (!ibd | age_event > tstart)

  pt <- data.frame(id = dat$f.eid, age_entry = dat$age_entry, tstart = tstart, tstop = age_exit,
                   event = as.integer(event), exposed = as.integer(ibd), dat[, covariates])
  pt <- pt[pt$tstop > pt$tstart, ]
  pt[order(pt$exposed), ]    # never-IBD rows first; the row order fixes the Monte Carlo draws
}

# Natural effects for a = 1 vs a* = 0, continuous mediator, rare outcome (VanderWeele 2015)
# b1: IBD on trait; s2: residual variance; t1, t2, t3: IBD, trait and IBD x trait in the Cox model
decompose <- function(b1, m0, s2, t1, t2, t3) {
  log_nde <- t1 + t3 * m0 + t3 * t2 * s2 + 0.5 * t3^2 * s2
  log_nie <- b1 * (t2 + t3)
  nde <- exp(log_nde); nie <- exp(log_nie)
  list(nde = nde, nie = nie, te = nde * nie, pm = nde * (nie - 1) / (nde * nie - 1),
       log_nie = log_nie)
}

# Point estimates from the fitted models; 95% CIs and P from draws of both models' coefficients
mc_decompose <- function(mfit, ofit, m0) {
  set.seed(seed)
  bm <- coef(mfit); bo <- coef(ofit)
  s2 <- summary(mfit)$sigma^2; df <- mfit$df.residual
  est <- decompose(bm["exposed"], m0, s2, bo["exposed"], bo["z"], bo["z:exposed"])

  # draws depend on the eigenvectors of vcov, i.e. on row order and on the BLAS/LAPACK used
  Bm <- mvrnorm(n_draws, bm, vcov(mfit))
  Bo <- mvrnorm(n_draws, bo, vcov(ofit))
  s2_draws <- s2 * df / rchisq(n_draws, df)
  drw <- decompose(Bm[, "exposed"], m0, s2_draws, Bo[, "exposed"], Bo[, "z"], Bo[, "z:exposed"])

  ci <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))
  pr <- mean(drw$log_nie <= 0)
  out <- list(p_nie = max(2 * min(pr, 1 - pr), 1 / n_draws))    # two-sided, never below 1/draws
  for (k in c("nde", "nie", "te", "pm"))
    out[paste0(k, c("", "_lci", "_uci"))] <- c(unname(est[[k]]), ci(drw[[k]]))
  out
}


######################## one trait and one outcome ########################
fit_pair <- function(pt, outcome, trait) {
  # trait is NA for participants not in the trait file (hematological diagnosis at baseline);
  # rank-inverse-normal within the outcome's sample
  z <- rint(blood[match(pt$id, blood$f.eid), trait])
  d <- pt[!is.na(z), ]
  d$z <- z[!is.na(z)]

  # landmark: follow-up starts 6 months (2 years for clonal outcomes) after baseline
  landmark <- ifelse(outcome %in% clonal, 2, 0.5)
  d$tstart <- pmax(d$tstart, d$age_entry + landmark)
  d <- d[d$tstop > d$tstart, ]

  res <- data.frame(outcome_id = outcome, trait_id = trait, landmark = landmark,
                    participants = nrow(d), events = sum(d$event),
                    events_ibd = sum(d$event[d$exposed == 1]),
                    a = NA, a_se = NA, b = NA, b_se = NA, p_interaction = NA,
                    nde = NA, nde_lci = NA, nde_uci = NA, nie = NA, nie_lci = NA, nie_uci = NA,
                    p_nie = NA, te = NA, te_lci = NA, te_uci = NA, pm = NA, pm_lci = NA, pm_uci = NA,
                    estimable = FALSE, note = "")
  if (res$events < min_events) { res$note <- "too few events"; return(res) }

  mfit <- lm(fml_m, data = d)                           # one row per participant
  ofit <- tryCatch(suppressWarnings(coxph(fml_o, data = d, ties = "efron")),
                   error = function(e) NULL)
  if (is.null(ofit)) { res$note <- "outcome model failed"; return(res) }

  sm <- summary(mfit)$coefficients; so <- summary(ofit)$coefficients
  res[c("a", "a_se")] <- sm["exposed", c("Estimate", "Std. Error")]
  res[c("b", "b_se", "p_interaction")] <- c(so["z", c("coef", "se(coef)")],
                                            so["z:exposed", "Pr(>|z|)"])

  # mean trait without IBD over the sample's covariates
  m0 <- mean(predict(mfit, newdata = transform(d, exposed = 0)))
  mc <- mc_decompose(mfit, ofit, m0)

  # direct effect outside its own interval: the IBD x trait term is too unstable to identify it
  if (!is.finite(mc$nde_uci) || mc$nde < mc$nde_lci || mc$nde > mc$nde_uci) {
    res$note <- "direct effect not identifiable"
    return(res)
  }
  res[names(mc)] <- mc
  res$estimable <- TRUE
  res
}


######################## all traits x the FDR-significant outcomes ########################
# about 450 Cox fits, roughly 2-3 hours on one core
results <- data.frame()
for (outcome in sig$outcome_id) {
  pt  <- make_person_time(dat, outcome)
  res <- mclapply(traits$trait_id, function(tr) fit_pair(pt, outcome, tr), mc.cores = n_cores)
  results <- rbind(results, do.call(rbind, res))
  print(paste(outcome, "done:", sum(results$estimable[results$outcome_id == outcome]), "estimable"))
}

# BH-FDR over all estimable pairs (nucleated red cells, 99% at the detection floor, left out)
results$fdr_family <- results$estimable &
  results$trait_id %in% traits$trait_id[traits$in_fdr_family]
results$FDR_P <- NA
results$FDR_P[results$fdr_family] <- p.adjust(results$p_nie[results$fdr_family], method = "BH")

# reported pairs: indirect-effect CI excludes 1, 0 < proportion < 1, its CI excludes 0, FDR < 0.05
results$meets_criteria <- results$fdr_family & results$FDR_P < 0.05 &
  (results$nie_lci > 1 | results$nie_uci < 1) &
  results$pm > 0 & results$pm < 1 & results$pm_lci > 0
results$definitional <- paste(results$outcome_id, results$trait_id) %in%
  paste(def_pair$outcome_id, def_pair$trait_id)

results$outcome <- outcomes$label[match(results$outcome_id, outcomes$outcome_id)]
results$trait   <- traits$label[match(results$trait_id, traits$trait_id)]

cat(sprintf("%d pairs, %d estimable, %d in the FDR family, %d meet all four criteria\n",
            nrow(results), sum(results$estimable), sum(results$fdr_family),
            sum(results$meets_criteria)))
cat(sprintf("Median proportion mediated across estimable pairs: %.1f%%\n",
            100 * median(results$pm[results$estimable])))


######################## Table S22 ########################
fmt_ci <- function(x, lo, hi, digits) {
  f <- paste0("%.", digits, "f")
  sprintf(paste0(f, " (", f, ", ", f, ")"), x, lo, hi)
}

tab <- subset(results, meets_criteria)
tab <- tab[order(tab$outcome_id, -tab$pm), ]
total <- sig[match(tab$outcome_id, sig$outcome_id), ]

tab <- data.frame(Outcome = tab$outcome,
                  `Total effect of IBD, HR (95% CI)` = fmt_ci(total$HR, total$LCI, total$UCI, 2),
                  Mediator = tab$trait,
                  `a, SD (95% CI)` = fmt_ci(tab$a, tab$a - 1.96 * tab$a_se,
                                            tab$a + 1.96 * tab$a_se, 3),
                  `b, log HR per SD (95% CI)` = fmt_ci(tab$b, tab$b - 1.96 * tab$b_se,
                                                       tab$b + 1.96 * tab$b_se, 3),
                  `Natural direct effect, HR (95% CI)` = fmt_ci(tab$nde, tab$nde_lci, tab$nde_uci, 3),
                  `Natural indirect effect, HR (95% CI)` = fmt_ci(tab$nie, tab$nie_lci, tab$nie_uci, 3),
                  `Proportion mediated (95% CI)` = sprintf("%.1f%% (%.1f, %.1f)", 100 * tab$pm,
                                                           100 * tab$pm_lci, 100 * tab$pm_uci),
                  `Definitional overlap` = ifelse(tab$definitional, "Yes", ""), check.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Table S22")
writeData(wb, "Table S22", tab)
addWorksheet(wb, "all pairs")
writeData(wb, "all pairs", results[, c("outcome", "trait",
                                       setdiff(names(results), c("outcome", "trait")))])
saveWorkbook(wb, file.path(out_dir, "Table_S22_mediation_blood_traits.xlsx"), overwrite = TRUE)

sessionInfo()
