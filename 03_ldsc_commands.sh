#!/usr/bin/env bash
############################ 03 LD score regression with ldsc.py ############################
# Genetic correlation of IBD, CD and UC with the 13 outcomes marked rg_method = ldsc.py in
# definitions/gwas_sources.csv (the other two are estimated in R/03_mr_ldsc.R).
#
# Input : data/gwas/   exposure and outcome GWAS, tab-separated, as described in R/03_mr_ldsc.R
#         data/ldsc/eur_w_ld_chr/, data/ldsc/w_hm3.snplist   European LD scores and HapMap3 SNP list,
#                      https://data.broadinstitute.org/alkesgroup/LDSCORE/ (eur_w_ld_chr.tar.bz2, w_hm3.snplist.bz2)
# Output: data/ldsc/logs/ldsc.<exposure>.<outcome>.log, read by R/03_mr_ldsc.R
#
# Design: munge_sumstats.py restricted to HapMap3 SNPs, then ldsc.py --rg with the intercepts constrained
#         (--no-intercept). Run before R/03_mr_ldsc.R, from the repository root.
# ldsc v1.0.1 (https://github.com/bulik/ldsc), Python 2.7; set LDSC to the ldsc folder
###################################################################################################
set -e
LDSC=$(cd "${LDSC:-ldsc}" && pwd)   # ldsc folder, resolved before changing directory
cd data/ldsc
mkdir -p logs


######################## munge input: SNP A1 A2 Z N P, A1 = effect allele ########################
# arguments: file in data/gwas, total sample size, column names of rsID, effect allele, other allele, beta, SE, P
# (N only scales the h2 estimates; rg does not depend on it)
prep () {
  gzip -dc "../gwas/$1" | awk -F '\t' -v N="$2" -v cols="$3" '
    BEGIN { OFS = "\t"; OFMT = "%.10g"; split(cols, k, ",") }
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "SNP", "A1", "A2", "Z", "N", "P"; next }
    ($c[k[5]] + 0) > 0 { print $c[k[1]], $c[k[2]], $c[k[3]], $c[k[4]] / $c[k[5]], N, $c[k[6]] }'
}
# --a1 A2 --a2 A1 swaps the allele labels of every file alike, so the sign of rg is unaffected
munge () {
  python "$LDSC/munge_sumstats.py" --sumstats "$1.txt" --out "$1" --chunksize 500000 \
    --merge-alleles w_hm3.snplist --snp SNP --a1 A2 --a2 A1 --ignore CHR,BP,BETA,SE
}


######################## exposures: de Lange 2017, European ########################
prep IBD17_EUR.tsv.gz 59957 SNP,A1,A2,BETA,SE,P > IBD17_EUR.txt; munge IBD17_EUR
prep CD17_EUR.tsv.gz  40266 SNP,A1,A2,BETA,SE,P > CD17_EUR.txt;  munge CD17_EUR
prep UC17_EUR.tsv.gz  45975 SNP,A1,A2,BETA,SE,P > UC17_EUR.txt;  munge UC17_EUR


######################## outcomes: FinnGen R10 and GWAS Catalog ########################
OUTCOMES=""
while read id n; do
  prep "$id.gz" "$n" rsids,alt,ref,beta,sebeta,pval > "$id.txt"; munge "$id"
  OUTCOMES="$OUTCOMES $id"
done <<'END'
finngen_R10_D3_ANAEMIA_IRONDEF 408837
finngen_R10_D3_THROMBOCYTOPENIANAS 407814
finngen_R10_CD2_PRIMARY_LYMPHOID_HEMATOPOIETIC_EXALLC 332813
finngen_R10_CD2_LYMPHOID_LEUKAEMIA_EXALLC 326267
finngen_R10_D3_HAEMOLYTICANAEMIA 412181
finngen_R10_D3_ANAEMIA 124542
finngen_R10_D3_COAGDEF_PURPUR_HAEMORRHAGIC 412181
finngen_R10_CD2_MULTIPLE_MYELOMA_PLASMA_CELL_EXALLC 325987
finngen_R10_CD2_HODGKIN_LYMPHOMA_EXALLC 325496
finngen_R10_C3_NONHODGKIN_EXALLC 315265
finngen_R10_CD2_MYELOID_LEUKAEMIA_EXALLC 324650
finngen_R10_MYELOPROF_NONCML 411923
END
prep GCST90018794.h.tsv.gz 473500 rsid,effect_allele,other_allele,beta,standard_error,p_value > GCST90018794.txt
munge GCST90018794
OUTCOMES="$OUTCOMES GCST90018794"


######################## genetic correlation, intercepts constrained ########################
for e in IBD CD UC; do
  for o in $OUTCOMES; do
    python "$LDSC/ldsc.py" --rg "${e}17_EUR.sumstats.gz,$o.sumstats.gz" --ref-ld-chr eur_w_ld_chr/ \
      --w-ld-chr eur_w_ld_chr/ --no-intercept --out "logs/ldsc.$e.$o"
  done
done
