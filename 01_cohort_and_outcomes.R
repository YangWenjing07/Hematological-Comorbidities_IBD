############################ 01 UK Biobank cohort, outcomes and blood traits ############################
#
# Input : data/ukb_analysis_dataset.Rdata       IBD_ICD: dates of baseline visit, birth, IBD diagnosis and
#                                               death; IBD subtype, covariates, blood traits, outcome dates
#         data/ukb_townsend_ltfu.Rdata          rq: Townsend index (f.189.0.0), date lost to follow-up
#                                               (f.191.0.0)
#         data/ukb_hospital_inpatient.Rdata     IBD: inpatient ICD-10 codes and dates (f.41270, f.41280),
#                                               ICD-9 codes and dates (f.41271, f.41281)
#         data/ukb_cancer_registry_hema.Rdata   hema_cancer: cancer-registry records ($Non_Hodgkins_lymphoma)
#         data/ukb_gp_immunosuppressants.Rdata  df1: primary-care immunosuppressant prescriptions (ATC L04A)
#         definitions/outcomes.csv, outcome_codes.csv (Table S1), covariates.csv, blood_traits.csv (Table S2)
# Output: data/ukb_analysis_cohort.rds          one row per participant
#         data/ukb_immunosuppressant.rds        one row per participant with a prescription
#         data/ukb_blood_traits.rds             participants without a hematological outcome at baseline
#
# Design: participants missing from the current data release (no baseline visit date) or with an IBD
#         diagnosis dated before birth are excluded; IBD diagnosed on or before the baseline visit is
#         prevalent, later IBD is incident; each outcome is dated by its earliest record in any source
#         and followed only while its sources are complete.

###################################################################################################

data_dir <- "data"

years    <- function(from, to) as.numeric(to - from) / 365.25
num_date <- function(x) as.Date(as.vector(x), origin = "1970-01-01")
load_object <- function(file, object) {
  e <- new.env()
  load(file.path(data_dir, file), envir = e)
  e[[object]]
}

codes      <- read.csv("definitions/outcome_codes.csv", colClasses = "character")
codes$included <- codes$included == "TRUE"
covariates <- read.csv("definitions/covariates.csv", colClasses = "character")
traits     <- read.csv("definitions/blood_traits.csv", colClasses = "character")
outcome_ids <- c(read.csv("definitions/outcomes.csv")$outcome_id,
                 "multiple_myeloma_dx")   # myeloma in hospital records alone, a check on the registry


######################## hospital inpatient records ########################
# Earliest date of an included code per participant and outcome, over every array position of the
# ICD-10 and ICD-9 fields; participants carrying a code listed as not included are flagged.
hosp <- load_object("ukb_hospital_inpatient.Rdata", "IBD")
scan <- subset(codes, grepl("hospital", source))
scan_ids <- unique(scan$outcome_id)

first_hosp    <- matrix(Inf,   nrow(hosp), length(scan_ids), dimnames = list(NULL, scan_ids))
excluded_code <- matrix(FALSE, nrow(hosp), length(scan_ids), dimnames = list(NULL, scan_ids))
for (system in c("ICD10", "ICD9")) {
  sys <- subset(scan, code_system == system)
  code_cols <- grep(ifelse(system == "ICD10", "^f\\.41270\\.", "^f\\.41271\\."), names(hosp), value = TRUE)
  for (col in code_cols) {
    rows <- which(hosp[[col]] %in% sys$code)
    if (length(rows) == 0) next
    code <- hosp[[col]][rows]
    date <- as.numeric(as.Date(hosp[[sub("^f\\.4127", "f.4128", col)]][rows]))   # matching date field
    for (o in scan_ids) {
      inc <- code %in% sys$code[sys$outcome_id == o & sys$included] & !is.na(date)
      first_hosp[rows[inc], o] <- pmin(first_hosp[rows[inc], o], date[inc])
      exc <- code %in% sys$code[sys$outcome_id == o & !sys$included]
      excluded_code[rows[exc], o] <- TRUE
    }
  }
}
first_hosp[is.infinite(first_hosp)] <- NA
hospital_end <- num_date(max(first_hosp, na.rm = TRUE))   # last date in these records
hosp_id <- hosp$f.eid
rm(hosp); invisible(gc())
cat("hospital inpatient records end on", format(hospital_end), "\n")


######################## outcome dates ########################
ukb <- load_object("ukb_analysis_dataset.Rdata", "IBD_ICD")
h <- match(ukb$f.eid, hosp_id)
hosp_date <- function(o) num_date(first_hosp[h, o])
dates <- list()

# hospital outcomes: the earlier of the inpatient records above (to 30 Sep 2021) and the hospital
# diagnosis date provided in the analysis dataset (records to 31 Oct 2022)
provided <- c(myeloid_leukaemia          = "Myeloid_leukemia_diag_date",
              myeloproliferative_disease = "Myelodysplastic_syndrome_diag_date",
              iron_deficiency_anaemia    = "Iron_deficiency_diag_date",
              pernicious_anaemia         = "Pernicious_anemia_diag_date",
              haemolytic_anaemia         = "hemolytic_anemias_diag_date",
              aplastic_anaemia           = "Aplastic_Anemias_diag_date",
              thrombocytopenia           = "Thrombocytopenia_diag_date",
              coagulation_defects        = "Coagulation_defects_diag_date",
              multiple_myeloma_dx        = "Multiple_Myeloma_diag_date")
for (o in names(provided)) dates[[o]] <- pmin(ukb[[provided[o]]], hosp_date(o), na.rm = TRUE)

# the provided myeloproliferative date (phecode 200) also counted codes that Table S1 leaves out:
# drop participants found only there whose hospital record carries such a code
only_excluded <- !is.na(ukb$Myelodysplastic_syndrome_diag_date) &
                 is.na(hosp_date("myeloproliferative_disease")) &
                 excluded_code[h, "myeloproliferative_disease"]
dates$myeloproliferative_disease[only_excluded] <- NA

# non-Hodgkin lymphoma: cancer-registry records with a Table S1 code, or hospital records
reg <- load_object("ukb_cancer_registry_hema.Rdata", "hema_cancer")$Non_Hodgkins_lymphoma
reg <- data.frame(f.eid = reg$f.eid, date = as.Date(reg$Non_Hodgkins_lymphoma_diag_date),
                  code = reg$Non_Hodgkins_lymphoma_diag_ICD)
reg <- subset(reg, !is.na(date))
registry_end <- max(reg$date)                                     # last date in the registry records
nhl <- subset(reg, code %in% codes$code[codes$outcome_id == "non_hodgkin_lymphoma" & codes$included])
dates$non_hodgkin_lymphoma <- pmin(nhl$date[match(ukb$f.eid, nhl$f.eid)],
                                   hosp_date("non_hodgkin_lymphoma"), na.rm = TRUE)

# venous thromboembolism: hospital inpatient records only
dates$vte <- hosp_date("vte")

# cancer-registry outcomes, as dated in the analysis dataset
dates$lymphoid_leukaemia <- ukb$Lymphoid_Leukaemia_cancer.date
dates$hodgkin_disease    <- ukb$Hodgkin_Disease_cancer.date
dates$multiple_myeloma   <- ukb$Multiple_Myeloma_cancer.date

# composites
dates$leukaemia_composite <- pmin(dates$myeloid_leukaemia, dates$lymphoid_leukaemia, na.rm = TRUE)
dates$lymphoma_composite  <- pmin(dates$hodgkin_disease, dates$non_hodgkin_lymphoma, na.rm = TRUE)
dates <- dates[outcome_ids]

# follow-up for each outcome ends when its sources stop being complete
complete_to <- rep(as.Date("2022-10-31"), length(outcome_ids))   # outcome dates provided in the dataset end here
names(complete_to) <- outcome_ids
complete_to[c("non_hodgkin_lymphoma", "lymphoma_composite")] <- registry_end
complete_to["vte"] <- hospital_end
rm(first_hosp, excluded_code, reg, nhl)


######################## cohort and exposure ########################
age_at_ibd_dx <- years(ukb$birthday.x, ukb$IBD_diag_date)
no_baseline   <- is.na(ukb$f.53.0.0)
dx_pre_birth  <- !no_baseline & !is.na(age_at_ibd_dx) & age_at_ibd_dx < 0
cat("participants:", nrow(ukb), "| not in the current release (no baseline date):", sum(no_baseline),
    "| IBD diagnosis dated before birth:", sum(dx_pre_birth), "\n")
keep  <- !no_baseline & !dx_pre_birth
ukb   <- ukb[keep, ]
dates <- lapply(dates, function(x) x[keep])

dat <- data.frame(f.eid = ukb$f.eid, baseline_date = ukb$f.53.0.0, birth_date = ukb$birthday.x,
                  ibd_subtype = ukb$IBD_subtype, ibd_diag_date = ukb$IBD_diag_date)
dat$exposure_group <- ifelse(is.na(dat$ibd_diag_date), "never",
                             ifelse(dat$ibd_diag_date <= dat$baseline_date, "prevalent", "incident"))
dat$ibd_subtype[dat$exposure_group == "never"] <- NA
age_at_ibd_dx <- years(dat$birth_date, dat$ibd_diag_date)
dat$flag_veo_ibd <- !is.na(age_at_ibd_dx) & age_at_ibd_dx < 5    # very early onset: flagged, not excluded
print(table(dat$exposure_group, dat$ibd_subtype, useNA = "ifany"))


######################## follow-up on the age scale ########################
rq <- load_object("ukb_townsend_ltfu.Rdata", "rq")
r  <- match(dat$f.eid, as.character(rq$f.eid))

dat$age_entry       <- years(dat$birth_date, dat$baseline_date)
dat$age_admin_end   <- years(dat$birth_date, as.Date("2022-12-31"))   # administrative end of follow-up
dat$age_death       <- years(dat$birth_date, ukb$death.date)
dat$age_ltfu        <- years(dat$birth_date, as.Date(rq$f.191.0.0[r]))
dat$age_switch      <- age_at_ibd_dx + 0.5     # exposed person-time starts 6 months after diagnosis
dat$age_dx_unlagged <- age_at_ibd_dx
dat$age_at_ibd_dx   <- age_at_ibd_dx


######################## covariates ########################
# categorical covariates, with the levels of definitions/covariates.csv (reference level first)
for (i in which(covariates$source_column != "" & covariates$levels != "")) {
  lev <- strsplit(covariates$levels[i], "|", fixed = TRUE)[[1]]
  lab <- if (covariates$labels[i] == "") lev else strsplit(covariates$labels[i], "|", fixed = TRUE)[[1]]
  x <- factor(as.character(ukb[[covariates$source_column[i]]]), levels = lev, labels = lab)
  stopifnot(!anyNA(x))
  dat[[covariates$covariate[i]]] <- x
}

# Townsend index: missing values set to the median, with a missing indicator
townsend <- rq$f.189.0.0[r]
dat$deprivation         <- ifelse(is.na(townsend), median(townsend, na.rm = TRUE), townsend)
dat$deprivation_missing <- factor(as.integer(is.na(townsend)), levels = 0:1)

dat$recruitment_age_int <- ukb$recruitment_Age
dat$is_baseline_flag    <- ukb$IS_all     # immunosuppressant use at baseline (primary care or self-report)


######################## outcome landmarks ########################
# date, age at the end of complete ascertainment, age at first record, record on or before baseline
for (o in outcome_ids) dat[[paste0("dt_", o)]] <- dates[[o]]
for (o in outcome_ids) dat[[paste0("age_comp_", o)]] <- years(dat$birth_date, complete_to[[o]])
for (o in outcome_ids) {
  dat[[paste0("age_", o)]]   <- years(dat$birth_date, dates[[o]])
  dat[[paste0("prior_", o)]] <- !is.na(dates[[o]]) & dates[[o]] <= dat$baseline_date
}

counts <- data.frame(outcome = outcome_ids,
                     dated = sapply(dates, function(x) sum(!is.na(x))),
                     prior = sapply(outcome_ids, function(o) sum(dat[[paste0("prior_", o)]])),
                     complete_to = format(complete_to), row.names = NULL)
counts$incident <- counts$dated - counts$prior
print(counts)

dat <- dat[, c("f.eid", "baseline_date", "birth_date", "ibd_subtype", "ibd_diag_date", "exposure_group",
               "flag_veo_ibd", "age_entry", "age_admin_end", "age_death", "age_ltfu", "age_switch",
               "age_dx_unlagged", "age_at_ibd_dx", covariates$covariate, "recruitment_age_int",
               "is_baseline_flag", paste0("dt_", outcome_ids), paste0("age_comp_", outcome_ids),
               as.vector(rbind(paste0("age_", outcome_ids), paste0("prior_", outcome_ids))))]


######################## immunosuppressant prescriptions ########################
# One record per participant and prescription date (drug names are not used); dates after the
# end of follow-up or before birth are dropped.
gp <- load_object("ukb_gp_immunosuppressants.Rdata", "df1")
gp <- unique(data.frame(f.eid = as.character(gp$f.eid), date = num_date(unclass(gp$issue_date))))
b  <- match(gp$f.eid, dat$f.eid)
gp <- gp[which(!is.na(b) & gp$date <= as.Date("2022-12-31") & gp$date >= dat$birth_date[b]), ]

first <- tapply(as.numeric(gp$date), gp$f.eid, min)
last  <- tapply(as.numeric(gp$date), gp$f.eid, max)
b <- match(names(first), dat$f.eid)
is_history <- data.frame(f.eid = names(first), first_is_date = num_date(first), last_is_date = num_date(last),
                         n_script_dates = as.vector(tapply(gp$date, gp$f.eid, length)))
days_from_baseline <- as.numeric(is_history$first_is_date - dat$baseline_date[b])
is_history$age_first_is       <- dat$age_entry[b] + days_from_baseline / 365.25
is_history$is_before_baseline <- is_history$first_is_date <= dat$baseline_date[b]
is_history$years_first_is_from_baseline <- days_from_baseline / 365.25
cat("participants with an immunosuppressant prescription:", nrow(is_history),
    "| first before baseline:", sum(is_history$is_before_baseline), "\n")


######################## blood traits ########################
# 61 assayed traits as measured (each analysis transforms within its own sample) and 4 ratios
blood <- data.frame(f.eid = dat$f.eid, profile_group = dat$exposure_group, ibd_subtype = dat$ibd_subtype,
                    age = dat$age_entry, dat[, covariates$covariate])

# covariates that the Chinese cohort can also supply (Table S4)
recode <- function(x, from, to) factor(to[match(as.character(x), from)], levels = unique(to))
blood$education_h <- recode(blood$education, c("NonCollege", "College", "Missing"),
                            c("non_college", "college", "unknown"))
blood$bmi_h       <- recode(blood$bmi, c("Normal", "Underweight", "Overweight", "Obese", "Missing"),
                            c("Normal", "Underweight", "Overweight", "Obese", "unknown"))
blood$smoking_h   <- recode(blood$smoking, c("1", "0", "miss"), c("never", "ever", "unknown"))
blood$alcohol_h   <- recode(blood$alcohol, c("1", "0", "miss"),
                            c("within_guideline", "above_guideline", "unknown"))
blood$is_baseline_flag <- dat$is_baseline_flag

measured <- subset(traits, ukb_field != "")
for (k in 1:nrow(measured)) blood[[measured$trait_id[k]]] <- ukb[[measured$ukb_field[k]]]
lym <- blood$lymphocyte_count
blood$sii <- ifelse(lym > 0, blood$neutrophil_count * blood$platelet_count / lym, NA)
blood$nlr <- ifelse(lym > 0, blood$neutrophil_count / lym, NA)
blood$plr <- ifelse(lym > 0, blood$platelet_count / lym, NA)
blood$lmr <- ifelse(blood$monocyte_count > 0, lym / blood$monocyte_count, NA)
stopifnot(all(traits$trait_id %in% names(blood)))

# drop participants with any hematological outcome (venous thromboembolism aside) on or before baseline
any_prior <- Reduce(`|`, dat[paste0("prior_", setdiff(outcome_ids, "vte"))])
print(table(group = blood$profile_group, excluded = any_prior))
blood <- blood[!any_prior, ]


######################## save ########################
rownames(dat) <- NULL; rownames(blood) <- NULL
saveRDS(dat,        file.path(data_dir, "ukb_analysis_cohort.rds"))
saveRDS(is_history, file.path(data_dir, "ukb_immunosuppressant.rds"))
saveRDS(blood,      file.path(data_dir, "ukb_blood_traits.rds"))
cat("cohort:", nrow(dat), "x", ncol(dat), "| blood traits:", nrow(blood), "x", ncol(blood), "\n")

sessionInfo()
