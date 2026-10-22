############################ 03 Mendelian randomization and genetic correlation ############################
# Two-sample Mendelian randomization (MR) and LD score regression (LDSC) of genetic liability to IBD,
# Crohn's disease (CD) and ulcerative colitis (UC) with 15 hematological outcomes (Tables S6-S8, S13, S14).
#
# Input : definitions/mr_instruments.csv   instruments with their exposure GWAS statistics
#         definitions/gwas_sources.csv     the 15 outcome GWAS (Table S5)
#         definitions/ldsc_rg_without_log.csv  rg for myeloid leukemia (ldsc.py log not kept)
#         data/gwas/<file>                 outcome GWAS as downloaded (file names in gwas_sources.csv): FinnGen R10,
#                                          or GWAS Catalog harmonised (columns rsid, effect_allele, other_allele,
#                                          beta, standard_error, p_value)
#         data/gwas/IBD17_EUR.tsv.gz       exposure GWAS, de Lange 2017 European, and CD17_EUR, UC17_EUR:
#                                          tab-separated SNP A1 A2 BETA SE P, A1 = effect allele,
#                                          rsIDs from 1000 Genomes
#         data/ldsc/logs/*.log             ldsc.py --rg logs written by R/03_ldsc_commands.sh
#         data/ldsc/eur_w_ld_chr/, data/ldsc/w_hm3.snplist   LDSC European LD scores and HapMap3 SNP list
# Output: output/Table_S6-S8_instruments.xlsx, output/Table_S13_genetic_correlation_mr.xlsx,
#         output/Table_S14_mr_sensitivity.xlsx
#
# Design: instruments P < 5e-8, clumped at r2 < 0.001 within 10,000 kb (1000 Genomes EUR); IVW, MR-Egger,
#         weighted median, weighted mode, simple mode, MR-Egger intercept, Cochran's Q and MR-PRESSO;
#         rg from LDSC with constrained intercepts (ldsc.py v1.0.1, or an R reimplementation of it for
#         the outcomes with rg_method = R); BH-FDR over the 15 outcomes within each exposure.
#  TwoSampleMR 0.7.9, MRPRESSO 1.0, data.table 1.17.8, R.utils 2.13.0, openxlsx 4.2.8.1
###################################################################################################

library(TwoSampleMR)
library(MRPRESSO)
library(data.table)
library(openxlsx)

out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)

instruments <- read.csv("definitions/mr_instruments.csv")
sources     <- read.csv("definitions/gwas_sources.csv")
exposures   <- c("IBD", "CD", "UC")
exposure_n  <- c(IBD = 59957, CD = 40266, UC = 45975)   # de Lange 2017 sample sizes, for LDSC

seed       <- 20260904   # bootstrap of weighted median/mode and MR-PRESSO
run_presso <- TRUE       # MR-PRESSO takes about 4 hours on one core (FALSE skips it)
presso_nb  <- 3000       # MR-PRESSO NbDistribution, must exceed (no. instruments) / 0.05

fmt_p  <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_or <- function(b, lo, hi) sprintf("%.2f (%.2f, %.2f)", exp(b), exp(lo), exp(hi))
fmt_4  <- function(x) ifelse(is.na(x), "-", sprintf("%.4f", x))


######################## instruments (Tables S6-S8) ########################
# Selected with TwoSampleMR::clump_data(clump_r2 = 0.001) on the variants with P < 5e-8; the effect allele
# is the allele the de Lange 2017 effect refers to (Allele2 in the IIBDGC files).
wb <- createWorkbook()
for (k in 1:3) {
  iv <- instruments[instruments$exposure == exposures[k], ]
  tab <- data.frame(SNP = iv$SNP, Chr = iv$chr, Pos = iv$pos, `Other allele` = iv$other_allele,
                    `Effect allele` = iv$effect_allele, Effect = iv$beta, SE = iv$se, `P value` = iv$p,
                    Exposure = iv$exposure, check.names = FALSE)
  addWorksheet(wb, paste0("Table S", 5 + k))
  writeData(wb, paste0("Table S", 5 + k), tab)
}
# instrument strength (Table S5): F = (beta / se)^2 per variant
f_stat <- data.frame()
for (ex in exposures) {
  f <- with(instruments[instruments$exposure == ex, ], (beta / se)^2)
  f_stat <- rbind(f_stat, data.frame(Exposure = ex, `No. instruments` = length(f),
                                     `Mean F-statistic` = round(mean(f), 1), `Minimum F-statistic` = round(min(f), 1),
                                     `% with F<10` = round(100 * mean(f < 10), 1), check.names = FALSE))
}
addWorksheet(wb, "F statistics")
writeData(wb, "F statistics", f_stat)
saveWorkbook(wb, file.path(out_dir, "Table_S6-S8_instruments.xlsx"), overwrite = TRUE)


######################## reading summary statistics ########################
# One GWAS file, restricted to the given SNPs; columns SNP, ref, alt (effect allele), beta, se, p.
read_gwas <- function(file, snps) {
  d <- fread(file.path("data/gwas", file), showProgress = FALSE)
  if ("rsids" %in% names(d)) {
    setnames(d, c("rsids", "sebeta", "pval"), c("SNP", "se", "p"))                       # FinnGen R10
  } else {
    setnames(d, c("rsid", "other_allele", "effect_allele", "standard_error", "p_value"),  # GWAS Catalog
             c("SNP", "ref", "alt", "se", "p"))
  }
  d[SNP %chin% snps, .(SNP, ref, alt, beta, se, p)]
}


######################## MR for one exposure and one outcome ########################
mr_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_simple_mode")
short      <- c("Inverse variance weighted" = "ivw", "MR Egger" = "egger", "Weighted median" = "wmedian",
                "Weighted mode" = "wmode", "Simple mode" = "smode")

run_mr <- function(iv, gwas) {
  exp_dat <- format_data(iv, type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se",
                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p")
  out_dat <- format_data(as.data.frame(gwas[SNP %chin% iv$SNP]), type = "outcome", snp_col = "SNP",
                         beta_col = "beta", se_col = "se", effect_allele_col = "alt", other_allele_col = "ref",
                         pval_col = "p")
  # alleles assumed to be on the forward strand; palindromic variants kept
  dat <- harmonise_data(exp_dat, out_dat, action = 1)
  dat <- dat[dat$mr_keep, ]

  set.seed(seed)
  est   <- mr(dat, method_list = mr_methods)
  egger <- mr_pleiotropy_test(dat)
  het   <- mr_heterogeneity(dat, method_list = "mr_ivw")

  res <- data.frame(nsnp = nrow(dat))
  for (k in 1:nrow(est)) {
    m <- short[[est$method[k]]]
    res[[paste0(m, "_b")]] <- est$b[k]; res[[paste0(m, "_se")]] <- est$se[k]; res[[paste0(m, "_p")]] <- est$pval[k]
  }
  res$egger_intercept    <- egger$egger_intercept
  res$egger_intercept_se <- egger$se
  res$egger_intercept_p  <- egger$pval
  res$Q <- het$Q; res$Q_df <- het$Q_df; res$Q_p <- het$Q_pval

  # MR-PRESSO: global test; outlier test and outlier-corrected IVW only when the global test is significant
  res$presso_global_p <- NA; res$presso_outliers <- NA; res$presso_b <- NA; res$presso_se <- NA; res$presso_p <- NA
  if (run_presso) {
    pr <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                    SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
                    NbDistribution = presso_nb, SignifThreshold = 0.05, seed = seed)
    main <- pr$`Main MR results`
    res$presso_global_p <- as.character(pr$`MR-PRESSO results`$`Global Test`$Pvalue)   # "<1/NbDistribution" if 0
    res$presso_outliers <- 0
    if (!is.na(main$`Causal Estimate`[2])) {
      res$presso_outliers <- length(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
      res$presso_b <- main$`Causal Estimate`[2]; res$presso_se <- main$Sd[2]; res$presso_p <- main$`P-value`[2]
    }
  }
  res
}


######################## MR: 3 exposures x 15 outcomes ########################
mr_results <- data.frame()
for (i in 1:nrow(sources)) {
  gwas <- read_gwas(sources$file[i], instruments$SNP)
  for (ex in exposures) {
    res <- run_mr(instruments[instruments$exposure == ex, ], gwas)
    mr_results <- rbind(mr_results, data.frame(exposure = ex, outcome_id = sources$outcome_id[i], res))
    print(paste(ex, sources$outcome_id[i], "done"))
  }
}

# BH-FDR over the 15 outcomes within each exposure
mr_results$ivw_fdr <- NA
for (ex in exposures) {
  k <- mr_results$exposure == ex
  mr_results$ivw_fdr[k] <- p.adjust(mr_results$ivw_p[k], method = "BH")
}


######################## genetic correlation: ldsc.py logs ########################
# The summary line of each ldsc.py --rg log: p1 p2 rg se z p ...
read_ldsc_log <- function(f) {
  L <- readLines(f)
  h <- grep("^p1\\s+p2\\s+rg", L)
  if (length(h) == 0) return(NULL)
  v <- strsplit(trimws(L[max(h) + 1]), "\\s+")[[1]]
  id <- sub("(\\.gz)?\\.sumstats\\.gz$", "", basename(v[2]))
  id <- sub(".*(GCST[0-9]+).*", "\\1", id)              # GWAS Catalog outcomes by accession
  data.frame(exposure = sub("17.*", "", basename(v[1])), outcome_id = id,
             rg = suppressWarnings(as.numeric(v[3])), rg_se = suppressWarnings(as.numeric(v[4])),
             rg_p = suppressWarnings(as.numeric(v[6])))
}
logs <- list.files("data/ldsc/logs", pattern = "\\.log$", full.names = TRUE)
rg <- do.call(rbind, lapply(logs, read_ldsc_log))
rg <- rg[rg$outcome_id %in% sources$outcome_id[sources$rg_method == "ldsc.py"], ]
# myeloid leukemia: the ldsc.py log of the original run was not kept; rg, SE and P as recorded from that run
no_log <- read.csv("definitions/ldsc_rg_without_log.csv")
rg <- rbind(rg, no_log[!paste(no_log$exposure, no_log$outcome_id) %in% paste(rg$exposure, rg$outcome_id), ])


######################## genetic correlation: constrained LDSC in R ########################
# For the outcomes with rg_method = R: an R reimplementation of the ldsc.py estimator with the h2 intercept
# fixed at 1 and the genetic covariance intercept at 0, i.e. weighted regressions through the origin on LD
# scores (one weight update; no genetic covariance term in the weights), rg = gcov / sqrt(h2_1 h2_2),
# SE from a 200-block jackknife in genome order.
ld_dir <- "data/ldsc/eur_w_ld_chr"
ld  <- rbindlist(lapply(1:22, function(k) fread(file.path(ld_dir, paste0(k, ".l2.ldscore.gz")))[, .(CHR, SNP, BP, L2)]))
M   <- sum(sapply(1:22, function(k) as.numeric(readLines(file.path(ld_dir, paste0(k, ".l2.M_5_50")))[1])))
hm3 <- fread("data/ldsc/w_hm3.snplist", col.names = c("SNP", "H1", "H2"))
comp <- c(A = "T", T = "A", C = "G", G = "C")

# munge_sumstats.py QC: SNVs, no strand-ambiguous pairs, unique rsIDs, HapMap3 SNPs with matching alleles
munge <- function(d, N) {   # d: SNP, A1 (effect allele), A2, Z
  d <- d[is.finite(Z) & A1 %chin% names(comp) & A2 %chin% names(comp)]
  d <- d[!(A2 == comp[A1])]
  d <- d[!(SNP %chin% d[duplicated(SNP), SNP])]
  d <- merge(d, hm3, by = "SNP")
  d <- d[(A1 == H1 & A2 == H2) | (A1 == H2 & A2 == H1) | (comp[A1] == H1 & comp[A2] == H2) |
         (comp[A1] == H2 & comp[A2] == H1), .(SNP, A1, A2, Z, N = as.numeric(N))]
  merge(d, ld, by = "SNP")[order(CHR, BP)]
}
chi2_ok <- function(z, N) z^2 < max(0.001 * N, 80)     # chi-square filter max(0.001 N, 80), as ldsc.py uses for h2
blocks  <- function(n) pmin(200L, 1L + floor(seq_len(n) * 200L / (n + 1)))
jack <- function(x, y, w, blk) {                       # slope through the origin, delete-one-block jackknife
  sxy <- w * x * y; sxx <- w * x * x
  ps <- (sum(sxy) - as.numeric(rowsum(sxy, blk))) / (sum(sxx) - as.numeric(rowsum(sxx, blk)))
  ps <- ps[is.finite(ps)]; g <- length(ps)
  c(est = sum(sxy) / sum(sxx), se = sqrt((g - 1) / g * sum((ps - mean(ps))^2)))
}
h2_ldsc <- function(d) {
  d <- d[chi2_ok(Z, N[1])]
  x <- d$L2 * d$N[1] / M; y <- d$Z^2 - 1; wld <- pmax(d$L2, 1)
  w0 <- 1 / wld
  h <- max(min(sum(w0 * x * y) / sum(w0 * x * x), 1), 0)          # starting value for the weights
  jack(x, y, 1 / (wld * (1 + h * x)^2), blocks(nrow(d)))
}
gcov_ldsc <- function(d1, d2, h1, h2) {
  m <- merge(d1[, .(SNP, A1, A2, Z1 = Z, N1 = N, CHR, BP, L2)], d2[, .(SNP, A1, A2, Z2 = Z, N2 = N)], by = "SNP")
  m <- m[(A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y)]
  m[, Z2 := fifelse(A1.x == A1.y, Z2, -Z2)]               # align to the exposure's effect allele
  m <- m[chi2_ok(Z1, N1[1]) & chi2_ok(Z2, N2[1])]
  setorder(m, CHR, BP)
  x  <- m$L2 * sqrt(m$N1[1] * m$N2[1]) / M
  v1 <- 1 + max(h1, 0) * m$L2 * m$N1[1] / M
  v2 <- 1 + max(h2, 0) * m$L2 * m$N2[1] / M
  jack(x, m$Z1 * m$Z2, 1 / (pmax(m$L2, 1) * v1 * v2), blocks(nrow(m)))
}
rg_ldsc <- function(d1, d2) {
  h1 <- h2_ldsc(d1)[["est"]]; h2 <- h2_ldsc(d2)[["est"]]
  g  <- gcov_ldsc(d1, d2, h1, h2)
  den <- sqrt(max(h1, 0) * max(h2, 0))
  data.frame(rg = g[["est"]] / den, rg_se = g[["se"]] / den, rg_p = 2 * pnorm(-abs(g[["est"]] / g[["se"]])))
}

exposure_ld <- list()
for (ex in exposures) {
  d <- fread(file.path("data/gwas", paste0(ex, "17_EUR.tsv.gz")))
  exposure_ld[[ex]] <- munge(d[, .(SNP, A1, A2, Z = BETA / SE)], exposure_n[[ex]])
}
for (i in which(sources$rg_method == "R")) {
  d <- read_gwas(sources$file[i], hm3$SNP)
  outcome_ld <- munge(d[, .(SNP, A1 = alt, A2 = ref, Z = beta / se)], sources$n[i])
  for (ex in exposures) {
    rg <- rbind(rg, data.frame(exposure = ex, outcome_id = sources$outcome_id[i],
                               rg_ldsc(exposure_ld[[ex]], outcome_ld)))
  }
}

# BH-FDR over the 15 outcomes within each exposure
rg$rg_fdr <- NA
for (ex in exposures) {
  k <- rg$exposure == ex & !is.na(rg$rg_p)
  rg$rg_fdr[k] <- p.adjust(rg$rg_p[k], method = "BH")
}


######################## Table S13: rg and IVW estimate ########################
res <- merge(mr_results, rg, by = c("exposure", "outcome_id"), all.x = TRUE)
res <- merge(res, sources[, c("outcome_id", "outcome", "source")], by = "outcome_id")
res$gwas   <- ifelse(res$source == "GWAS Catalog", paste(res$source, res$outcome_id), res$source)
res$ivw_or <- fmt_or(res$ivw_b, res$ivw_b - 1.96 * res$ivw_se, res$ivw_b + 1.96 * res$ivw_se)

s13 <- res[order(match(res$exposure, exposures), res$rg_p), ]
tab13 <- data.frame(Exposure = s13$exposure, Outcome = s13$outcome, `Outcome GWAS` = s13$gwas,
                    `LDSC regression: genetic correlation rg` = fmt_4(s13$rg), `rg SE` = fmt_4(s13$rg_se),
                    `rg P value` = fmt_p(s13$rg_p), `rg FDR-adjusted P (15-outcome panel)` = fmt_p(s13$rg_fdr),
                    `Mendelian randomization: no. instruments` = s13$nsnp, `MR IVW OR (95% CI)` = s13$ivw_or,
                    `MR IVW P value` = fmt_p(s13$ivw_p), `MR FDR-adjusted P (15-outcome panel)` = fmt_p(s13$ivw_fdr),
                    check.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Table S13")
writeData(wb, "Table S13", tab13)
addWorksheet(wb, "all estimates")
writeData(wb, "all estimates", s13[, c("exposure", "outcome_id", "outcome", "rg", "rg_se", "rg_p", "rg_fdr",
                                       "nsnp", "ivw_b", "ivw_se", "ivw_p", "ivw_fdr")])
saveWorkbook(wb, file.path(out_dir, "Table_S13_genetic_correlation_mr.xlsx"), overwrite = TRUE)


######################## Table S14: MR sensitivity analyses ########################
s14 <- res[order(match(res$exposure, exposures), res$ivw_p), ]
or_ci <- function(m) {
  b <- s14[[paste0(m, "_b")]]; se <- s14[[paste0(m, "_se")]]
  fmt_or(b, b - 1.96 * se, b + 1.96 * se)
}

# MR-PRESSO outlier-corrected OR: t-based CI on (variants kept - 1) df, as in its regression
global_p  <- suppressWarnings(as.numeric(sub("^<", "", s14$presso_global_p)))
corrected <- !is.na(s14$presso_b)
no_estimate <- ifelse(is.na(global_p), "not run", ifelse(global_p < 0.05, "no significant outlier", "not applicable"))
t_crit <- qt(0.975, s14$nsnp - s14$presso_outliers - 1)
presso_or <- ifelse(corrected, fmt_or(s14$presso_b, s14$presso_b - t_crit * s14$presso_se,
                                      s14$presso_b + t_crit * s14$presso_se), no_estimate)

tab14 <- data.frame(Exposure = s14$exposure, Outcome = s14$outcome, `Outcome GWAS` = s14$gwas, `No. SNPs` = s14$nsnp,
                    `IVW OR (95% CI)` = s14$ivw_or, `IVW P value` = fmt_p(s14$ivw_p),
                    `IVW FDR-adjusted P (panel)` = fmt_p(s14$ivw_fdr),
                    `Weighted median OR (95% CI)` = or_ci("wmedian"), `Weighted median P value` = fmt_p(s14$wmedian_p),
                    `Weighted mode OR (95% CI)` = or_ci("wmode"), `Weighted mode P value` = fmt_p(s14$wmode_p),
                    `MR-Egger slope OR (95% CI)` = or_ci("egger"), `MR-Egger slope P value` = fmt_p(s14$egger_p),
                    `MR-Egger intercept` = sprintf("%.4f", s14$egger_intercept),
                    `MR-Egger intercept SE` = sprintf("%.4f", s14$egger_intercept_se),
                    `MR-Egger intercept P value` = fmt_p(s14$egger_intercept_p),
                    `Cochran's Q` = sprintf("%.1f", s14$Q), `Cochran's Q df` = s14$Q_df,
                    `Cochran's Q P value` = fmt_p(s14$Q_p), `MR-PRESSO global P value` = fmt_p(global_p),
                    `MR-PRESSO outliers removed` = ifelse(is.na(s14$presso_outliers), "-", s14$presso_outliers),
                    `MR-PRESSO outlier-corrected OR (95% CI)` = presso_or,
                    `MR-PRESSO outlier-corrected P value` = ifelse(corrected, fmt_p(s14$presso_p), no_estimate),
                    check.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Table S14")
writeData(wb, "Table S14", tab14)
addWorksheet(wb, "all estimates")
writeData(wb, "all estimates", s14[, c("exposure", "outcome_id", "outcome", names(mr_results)[-(1:2)])])
saveWorkbook(wb, file.path(out_dir, "Table_S14_mr_sensitivity.xlsx"), overwrite = TRUE)

sessionInfo()
