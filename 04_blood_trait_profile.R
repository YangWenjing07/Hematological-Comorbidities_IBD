############################ 04 Blood-trait profile of IBD ############################
# Blood cell indices and biochemical markers in IBD, ulcerative colitis (UC) and Crohn's
# disease (CD) compared with people without IBD, in UK Biobank and in the Chinese multicenter
# cohort (Tables S15-S19; the estimates drawn in Figure 4).
#
# Input : data/cn_cohort.xlsx            Chinese multicenter cohort, one row per participant
#         data/ukb_blood_traits.rds      UK Biobank traits and covariates (01_cohort_and_outcomes.R)
#         data/ukb_analysis_cohort.rds   age at IBD diagnosis (disease-duration strata)
#         definitions/blood_trait_crosswalk.csv   65 traits: labels, Chinese columns, plausible ranges
# Output: data/cn_blood_traits.rds       Chinese analysis table, without identifiers
#         output/Table_S15_baseline_cn.xlsx, output/Table_S16_blood_traits_ukb.xlsx,
#         output/Table_S17_blood_traits_cn.xlsx, output/Table_S18_age_duration_sensitivity.xlsx,
#         output/Table_S19_blood_traits_ukb_harmonized.xlsx
#
# Design: rank-based inverse-normal transformed trait regressed on exposure and covariates (age as a
#         3-df natural spline); compared if >= 300 measured, >= 100 per group and <= 30 points apart in
#         % measured (Chinese cohort: within centers enrolling both groups, center adjusted); BH-FDR per
#         trait class and exposure. S18: UKB 40-50 y vs SRRSH >= 40 y; UKB IBD by years since diagnosis
#         (harmonized covariates). S19: UKB refitted with the covariates available in both cohorts.

###################################################################################################

library(openxlsx)

data_dir <- "data"
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE)

traits <- read.csv("definitions/blood_trait_crosswalk.csv", stringsAsFactors = FALSE)
cn_ids <- traits$trait_id[traits$in_cn]          # the 52 traits measured in the Chinese cohort

min_n     <- 300    # fewer measured values than this: not compared
min_group <- 100    # fewer measured values in either group than this: not compared
max_gap   <- 30     # larger difference (percentage points) in the proportion measured: not compared
age_term  <- "splines::ns(age, df = 3)"

ukb_covariates  <- c("sex", "ethnicity", "education", "deprivation", "deprivation_missing",
                     "bmi", "smoking", "alcohol", "diet", "physical_activity")
harm_covariates <- c("sex", "education_h", "bmi_h", "smoking_h", "alcohol_h")   # in both cohorts
cn_covariates   <- c(harm_covariates, "medication_h")


######################## Chinese cohort: analysis table ########################
raw <- read.xlsx(file.path(data_dir, "cn_cohort.xlsx"), sheet = 1)[, -1]   # drop the admission number

num <- function(x) suppressWarnings(as.numeric(x))
int <- function(x) suppressWarnings(as.integer(x))

# results below the detection limit (e.g. "<5") are set to half the limit
parse_value <- function(x) {
  x <- trimws(as.character(x))
  ifelse(grepl("^<\\s*[0-9.]+$", x), num(sub("^<\\s*", "", x)) / 2, num(x))
}

# implausible values set to missing: outside the trait's range or, for traits without a range,
# above 20 times the median (median of the non-zero values when the median is 0)
plausible <- function(x, lower = NA, upper = NA) {
  if (is.na(lower)) {
    med <- median(x, na.rm = TRUE)
    if (!(is.finite(med) && med > 0)) med <- median(x[!is.na(x) & x > 0])
    lower <- -Inf; upper <- if (is.finite(med) && med > 0) 20 * med else Inf
  }
  x[!is.na(x) & (x < lower | x > upper)] <- NA
  x
}

cn <- data.frame(IBD = as.integer(raw$IBD), centre = tolower(raw$source), src_block = raw$src_block)
cn$subtype <- ifelse(cn$IBD == 1, ifelse(is.na(raw$subtype), "IBD-U", raw$subtype), NA)

# hemoglobin is stored in two columns (two extractions) that never overlap
hb1 <- num(raw$Haemoglobin_concentration); hb2 <- num(raw$Hemoglobin_content)
stopifnot(!any(!is.na(hb1) & !is.na(hb2)))
raw$Haemoglobin_concentration <- ifelse(is.na(hb1), hb2, hb1)

for (i in which(traits$cn_column != "")) {
  cn[[traits$trait_id[i]]] <- plausible(parse_value(raw[[traits$cn_column[i]]]),
                                        traits$cn_lower[i], traits$cn_upper[i])
}

cn[cn$centre == "henan", c("lpa", "ast")] <- NA   # Henan: Lp(a) in two units, AST not in U/L

# platelet count not recorded: rebuilt from plateletcrit = platelet count x MPV / 10000
pl <- traits$trait_id == "platelet_count"
rebuild <- is.na(cn$platelet_count) & !is.na(cn$platelet_crit) & !is.na(cn$mpv)
cn$platelet_count[rebuild] <- plausible(cn$platelet_crit[rebuild] * 10000 / cn$mpv[rebuild],
                                        traits$cn_lower[pl], traits$cn_upper[pl])

# derived ratios, defined as in UK Biobank
ratio <- function(a, b) ifelse(b > 0, a / b, NA)
cn$sii <- ratio(cn$neutrophil_count * cn$platelet_count, cn$lymphocyte_count)
cn$nlr <- ratio(cn$neutrophil_count, cn$lymphocyte_count)
cn$plr <- ratio(cn$platelet_count, cn$lymphocyte_count)
cn$lmr <- ratio(cn$lymphocyte_count, cn$monocyte_count)

# covariates, each with an explicit unknown level
cn$age <- plausible(num(raw$recruitment_age), 1, 100)
cn$sex <- factor(ifelse(int(raw$gender) == 1, "Male", "Female"), levels = c("Female", "Male"))
edu <- ifelse(raw$education %in% c("college", "non_college"), raw$education, "unknown")
cn$education_h <- factor(edu, levels = c("non_college", "college", "unknown"))
# smoking and alcohol coded 1 never, 2 former, 3 current
never_ever <- function(x) ifelse(int(x) %in% 1, "never", ifelse(int(x) %in% 2:3, "ever", "unknown"))
cn$smoking_h <- factor(never_ever(raw$smoking), levels = c("never", "ever", "unknown"))
cn$alcohol_h <- factor(never_ever(raw$alcohol), levels = c("never", "ever", "unknown"))
gcs <- int(raw$glucocorticosteroid); ims <- int(raw$immunosuppressive)   # medication flags
med <- ifelse(gcs %in% 1 | ims %in% 1, "yes", ifelse(gcs %in% 0 & ims %in% 0, "no", "unknown"))
cn$medication_h <- factor(med, levels = c("no", "yes", "unknown"))
bmi    <- plausible(num(raw$BMI), 10, 70)
height <- plausible(num(raw$height), 100, 230)
weight <- plausible(num(raw$weight), 20, 250)
cn$bmi_value <- b <- ifelse(is.na(bmi), weight / (height / 100)^2, bmi)
b <- ifelse(is.na(b), "unknown", ifelse(b < 18.5, "Underweight", ifelse(b < 25, "Normal",
            ifelse(b < 30, "Overweight", "Obese"))))
cn$bmi_h <- factor(b, levels = c("Normal", "Underweight", "Overweight", "Obese", "unknown"))

# participants without a plausible age are dropped
keep <- c("IBD", "subtype", "centre", "src_block", "age", "sex", "education_h", "bmi_h", "bmi_value",
          "smoking_h", "alcohol_h", "medication_h", traits$trait_id[traits$cn_column != ""],
          "sii", "nlr", "plr", "lmr")
cn <- cn[!is.na(cn$age), keep]
rownames(cn) <- NULL
saveRDS(cn, file.path(data_dir, "cn_blood_traits.rds"))


######################## Table S15: baseline characteristics, Chinese cohort ########################
fmt_n   <- function(x) formatC(x, big.mark = ",", format = "d")
mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
n_pct   <- function(x, lev) sprintf("%s (%.1f)", fmt_n(sum(x == lev, na.rm = TRUE)),
                                    100 * sum(x == lev, na.rm = TRUE) / length(x))

s15 <- subset(cn, !is.na(wbc_count))                       # participants with a blood count
s15 <- split(s15, s15$IBD)
row15 <- function(label, f) data.frame(a = label, b = f(s15[["0"]]), c = f(s15[["1"]]))
block <- function(title, variable, labels) {      # a heading row, then one row per level (names)
  rows <- row15(title, function(g) "")
  for (lev in names(labels)) {
    rows <- rbind(rows, row15(paste0("    ", labels[[lev]]), function(g) n_pct(g[[variable]], lev)))
  }
  rows
}
tab15 <- rbind(
  row15("Age (years), mean (SD)", function(g) mean_sd(g$age)),
  row15("Aged under 18 years, n (%)", function(g) n_pct(g$age < 18, TRUE)),
  block("Sex, n (%)", "sex", c(Female = "Female", Male = "Male")),
  block("Education, n (%)", "education_h",
        c(college = "College", non_college = "Non-college", unknown = "Unknown")),
  row15("Body mass index (kg/m\u00b2), mean (SD)", function(g) mean_sd(g$bmi_value)),
  block("Body mass index (kg/m\u00b2), n (%)", "bmi_h",
        c(Underweight = "Underweight (<18.5)", Normal = "Normal (18.5-24.9)",
          Overweight = "Overweight (25-29.9)", Obese = "Obese (\u226530)", unknown = "Unknown")),
  block("Smoking status, n (%)", "smoking_h",
        c(never = "Never", ever = "Ever (former or current)", unknown = "Unknown")),
  block("Alcohol intake, n (%)", "alcohol_h", c(never = "Never", ever = "Ever", unknown = "Unknown")),
  block("Immunosuppressant use, n (%)", "medication_h", c(no = "No", yes = "Yes", unknown = "Unknown")),
  block("Recruitment center, n (%)", "centre", c(henan = "Henan", srrsh = "SRRSH", xiangya = "Xiangya")))
names(tab15) <- c("Characteristic", sprintf("Controls without IBD (n=%s)", fmt_n(nrow(s15[["0"]]))),
                  sprintf("IBD (n=%s)", fmt_n(nrow(s15[["1"]]))))
write.xlsx(tab15, file.path(out_dir, "Table_S15_baseline_cn.xlsx"), sheetName = "Table S15")


######################## the profiling model ########################
rint <- function(x) qnorm((rank(x) - 0.5) / length(x))     # rank-based inverse-normal transformation

# One trait, one comparison. d has a 0/1 column `exposed`; by_centre = TRUE judges eligibility
# within the centers that measured the trait in >= 100 people in both groups and adds a center term.
fit_trait <- function(d, trait, covariates, by_centre = FALSE) {
  x  <- d[[trait]]
  dd <- d[!is.na(x), c("exposed", "age", if (by_centre) "centre", covariates, trait)]
  res <- data.frame(trait_id = trait, n = nrow(dd), n_exposed = sum(dd$exposed == 1),
                    n_unexposed = sum(dd$exposed == 0), zero_pct = round(100 * mean(dd[[trait]] == 0), 2),
                    beta = NA, se = NA, lci = NA, uci = NA, p = NA, note = "")
  scope <- rep(TRUE, nrow(d))
  if (by_centre) {
    n1 <- table(dd$centre[dd$exposed == 1]); n0 <- table(dd$centre[dd$exposed == 0])
    scope <- d$centre %in% names(n1)[n1 >= min_group & n0 >= min_group]
  }
  measured1 <- 100 * mean(!is.na(x[scope & d$exposed == 1]))
  measured0 <- 100 * mean(!is.na(x[scope & d$exposed == 0]))
  if (!any(scope)) {
    res$note <- "Not compared: no center with 100 measured values in both groups"
  } else if (res$n < min_n || (!by_centre && min(res$n_exposed, res$n_unexposed) < min_group)) {
    res$note <- sprintf("Not compared: fewer than 100 measured values per group (n=%s; IBD %s, controls %s)",
                        fmt_n(res$n), fmt_n(res$n_exposed), fmt_n(res$n_unexposed))
  } else if (abs(measured1 - measured0) > max_gap) {
    res$note <- sprintf("Not compared: measured in %.0f%% of patients with IBD vs %.0f%% of controls",
                        measured1, measured0)
  } else if (by_centre && length(unique(dd$centre)) < 2) {
    res$note <- "Not estimable: measured in one center only, no center term possible"
  }
  if (res$note != "") return(res)
  dd$z <- rint(dd[[trait]])
  terms <- c(age_term, covariates, if (by_centre) "centre")
  fit <- lm(as.formula(paste("z ~ exposed +", paste(terms, collapse = " + "))), data = dd)
  s <- summary(fit)$coefficients["exposed", ]
  res$beta <- s[["Estimate"]]; res$se <- s[["Std. Error"]]; res$p <- s[["Pr(>|t|)"]]
  res$lci  <- res$beta - qnorm(0.975) * res$se                 # Wald 95% CI
  res$uci  <- res$beta + qnorm(0.975) * res$se
  res
}

run_traits <- function(d, trait_ids, covariates, model, exposure, by_centre = FALSE) {
  out <- data.frame()
  for (t in trait_ids) out <- rbind(out, fit_trait(d, t, covariates, by_centre))
  print(paste(model, exposure, "done"))
  data.frame(model = model, exposure = exposure, out)
}


######################## IBD, UC and CD vs no IBD (Tables S16, S17, S19) ########################
ukb <- as.data.frame(readRDS(file.path(data_dir, "ukb_blood_traits.rds")))
ukb <- subset(ukb, profile_group %in% c("never", "prevalent"))   # diagnosed during follow-up: not used

results <- data.frame()
for (exposure in c("IBD", "UC", "CD")) {
  # UK Biobank, primary and harmonized covariates; "IBD" subtype = IBD unclassified
  subtypes <- if (exposure == "IBD") c("UC", "CD", "IBD") else exposure
  d <- subset(ukb, profile_group == "never" | ibd_subtype %in% subtypes)
  d$exposed <- as.integer(d$profile_group == "prevalent")
  results <- rbind(results, run_traits(d, traits$trait_id, ukb_covariates, "UKB primary", exposure),
                   run_traits(d, traits$trait_id, harm_covariates, "UKB harmonized", exposure))
  # Chinese cohort, all centers
  d <- if (exposure == "IBD") cn else subset(cn, IBD == 0 | subtype == exposure)
  d$exposed <- d$IBD
  d$centre  <- factor(d$centre)
  results <- rbind(results, run_traits(d, cn_ids, cn_covariates, "China all centers", exposure,
                                       by_centre = TRUE))
}


######################## Table S18: age-restricted and duration comparisons ########################
# the ages where the two cohorts overlap: UK Biobank 40-50 years, SRRSH 40 years or older
d <- subset(ukb, age >= 40 & age <= 50)
d$exposed <- as.integer(d$profile_group == "prevalent")
results <- rbind(results, run_traits(d, traits$trait_id, ukb_covariates, "UKB aged 40-50", "IBD"))
d <- subset(cn, centre == "srrsh" & age >= 40)
d$exposed <- d$IBD
results <- rbind(results, run_traits(d, cn_ids, cn_covariates, "SRRSH aged 40+", "IBD"))

# UK Biobank IBD by years since diagnosis at baseline, each stratum vs never IBD
base <- readRDS(file.path(data_dir, "ukb_analysis_cohort.rds"))
duration <- ukb$age - base$age_at_ibd_dx[match(ukb$f.eid, base$f.eid)]
rm(base)
strata <- c("<5 years", "5-15 years", ">=15 years")
ukb$stratum <- cut(duration, c(0, 5, 15, Inf), labels = strata, right = FALSE)
for (s in strata) {
  d <- subset(ukb, profile_group == "never" | (profile_group == "prevalent" & stratum %in% s))
  d$exposed <- as.integer(d$profile_group == "prevalent")
  model <- paste("UKB duration", s)
  results <- rbind(results, run_traits(d, traits$trait_id, harm_covariates, model, "IBD"))
}


######################## BH-FDR ########################
k <- match(results$trait_id, traits$trait_id)
info <- c("label", "class", "category", "category_order")
results[info] <- traits[k, info]
results$in_family <- !is.na(results$beta) & traits$in_fdr_family[k]   # not: nucleated red cells
results$FDR_P <- NA
family <- paste(results$model, results$exposure, results$class)
for (f in unique(family)) {
  i <- family == f & results$in_family
  results$FDR_P[i] <- p.adjust(results$p[i], method = "BH")
}
results$significant <- !is.na(results$FDR_P) & results$FDR_P < 0.05

# UK Biobank FDR-significant traits: direction in the Chinese cohort, and under harmonized covariates
pick <- function(model, exposure, ids) {
  r <- results[results$model == model & results$exposure == exposure, ]
  r[match(ids, r$trait_id), ]
}
for (exposure in c("IBD", "UC", "CD")) {
  u  <- pick("UKB primary", exposure, traits$trait_id)
  ch <- pick("China all centers", exposure, traits$trait_id)    # NA where not measured in China
  h  <- pick("UKB harmonized", exposure, traits$trait_id)
  k  <- u$significant & !is.na(ch$beta)
  same <- sign(u$beta[k]) == sign(ch$beta[k]); sig <- ch$significant[k]
  cat(sprintf("%s: %d of %d same direction in China", exposure, sum(same), sum(k)),
      sprintf("(%d FDR-significant there; %d significant, opposite)\n", sum(same & sig), sum(!same & sig)))
  lost <- u$label[u$significant & !h$significant]
  cat(sprintf("     harmonized covariates (all traits): %d sign changes; no longer FDR-significant: %s\n",
              sum(sign(h$beta) != sign(u$beta), na.rm = TRUE), paste(lost, collapse = ", ")))
}


######################## Tables S16-S19 ########################
fmt_p    <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_beta <- function(b, l, u) ifelse(is.na(b), "", sprintf("%.3f (%.3f, %.3f)", b, l, u))
beta_col <- "\u03b2 (95% CI), SD"

profile_table <- function(model) {
  r <- results[results$model == model & !is.na(results$beta), ]          # traits compared
  r <- r[order(match(r$exposure, c("IBD", "UC", "CD")), r$category_order, -abs(r$beta)), ]
  tab <- data.frame(Exposure = r$exposure, `Trait category` = r$category, Trait = r$label,
                    `N analyzed` = fmt_n(r$n), `N exposed` = fmt_n(r$n_exposed),
                    b = fmt_beta(r$beta, r$lci, r$uci), `P value` = fmt_p(r$p),
                    `FDR-adjusted P` = ifelse(r$in_family, fmt_p(r$FDR_P), "-"),
                    `Zero values, %` = sprintf("%.1f", r$zero_pct), check.names = FALSE)
  names(tab)[6] <- beta_col
  tab
}

save_table <- function(sheets, models, file) {        # the table(s), plus the numeric estimates
  wb <- createWorkbook()
  for (s in names(sheets)) { addWorksheet(wb, s); writeData(wb, s, sheets[[s]]) }
  addWorksheet(wb, "all estimates")
  writeData(wb, "all estimates", results[results$model %in% models, ])
  saveWorkbook(wb, file.path(out_dir, file), overwrite = TRUE)
}
save_table(list("Table S16" = profile_table("UKB primary")), "UKB primary",
           "Table_S16_blood_traits_ukb.xlsx")
save_table(list("Table S17" = profile_table("China all centers")[, -9]), "China all centers",
           "Table_S17_blood_traits_cn.xlsx")          # zero values: UK Biobank only
save_table(list("Table S19" = profile_table("UKB harmonized")[, -9]), "UKB harmonized",
           "Table_S19_blood_traits_ukb_harmonized.xlsx")

# Table S18: traits in the order of Table S16 (IBD); "same direction" where estimable in both
ref <- results[results$model == "UKB primary" & results$exposure == "IBD", ]
ref <- ref[order(ref$category_order, -abs(ref$beta)), ]
same_sign <- function(a, b) ifelse(is.na(a) | is.na(b), "", ifelse(sign(a) == sign(b), "Yes", "No"))
add_cols <- function(tab, r, prefix, n_analyzed = TRUE) {
  if (n_analyzed) tab[[paste(prefix, "N analyzed")]] <- fmt_n(r$n)
  tab[[paste(prefix, "N with IBD")]] <- fmt_n(r$n_exposed)
  tab[[paste(prefix, beta_col)]] <- fmt_beta(r$beta, r$lci, r$uci)
  tab[[paste(prefix, "FDR-adjusted P")]] <- fmt_p(r$FDR_P)
  tab
}

ids <- ref$trait_id[ref$trait_id %in% cn_ids]
u <- pick("UKB aged 40-50", "IBD", ids); ch <- pick("SRRSH aged 40+", "IBD", ids)
sig <- pick("UKB primary", "IBD", ids)$significant
age_tab <- data.frame(`Trait category` = u$category, Trait = u$label, check.names = FALSE)
age_tab <- add_cols(add_cols(age_tab, u, "UK Biobank 40-50 y:"), ch, "SRRSH >=40 y:")
age_tab$`FDR-significant in the primary UK Biobank analysis` <- ifelse(sig, "Yes", "No")
age_tab$`Same direction` <- same <- same_sign(u$beta, ch$beta)
age_tab$Note <- ch$note
cat(sprintf("Age-restricted: %d of %d UKB FDR-significant traits\n", sum(same[sig] == "Yes"),
            sum(same[sig] != "")))

dur <- lapply(paste("UKB duration", strata), function(m) pick(m, "IBD", ref$trait_id))
dur_tab <- data.frame(`Trait category` = ref$category, Trait = ref$label, check.names = FALSE)
for (j in 1:3) dur_tab <- add_cols(dur_tab, dur[[j]], paste0(strata[j], ":"), n_analyzed = FALSE)
dur_tab$`Same direction, <5 vs >=15 years` <- same <- same_sign(dur[[1]]$beta, dur[[3]]$beta)
notes <- sapply(dur, function(r) r$note)                        # one column per stratum
dur_tab$Note <- apply(notes, 1, function(x) paste(unique(x[x != ""]), collapse = "; "))
cat(sprintf("Disease duration, <5 vs >=15 years: %d of %d traits\n", sum(same == "Yes"), sum(same != "")))

save_table(list("age-restricted" = age_tab, "disease duration" = dur_tab),
           c("UKB aged 40-50", "SRRSH aged 40+", paste("UKB duration", strata)),
           "Table_S18_age_duration_sensitivity.xlsx")

sessionInfo()
