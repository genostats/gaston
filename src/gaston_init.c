#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  generated with tools::package_native_routine_registration_skeleton
*/

/* .C calls */
extern void qfc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP gg_AIREML1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREML1_contrast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREML1_logit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREML1_logit_nofix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREML1_nofix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREMLn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREMLn_contrast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREMLn_logit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREMLn_logit_nofix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_AIREMLn_nofix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_alleles_duplicated(SEXP, SEXP);
extern SEXP gg_alleles_recoding(SEXP);
extern SEXP gg_as_matrix4(SEXP);
extern SEXP gg_bind_inds2(SEXP, SEXP);
extern SEXP gg_bind_snps(SEXP);
extern SEXP gg_diago_full_likelihood1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_full_likelihood1_nocovar(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_full_likelihood2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_full_likelihood2_nocovar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_likelihood1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_likelihood1_nocovar(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_likelihood2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_diago_likelihood2_nocovar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_duplicated_remove(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_extract_inds_bool(SEXP, SEXP);
extern SEXP gg_extract_inds_indices(SEXP, SEXP);
extern SEXP gg_extract_snps_bool(SEXP, SEXP);
extern SEXP gg_extract_snps_indices(SEXP, SEXP);
extern SEXP gg_fit_diago(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_fit_diago_nocovar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_geno_stats(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_geno_stats_inds(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_geno_stats_snps(SEXP, SEXP, SEXP);
extern SEXP gg_GWAS_lm_quanti(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_GWAS_lmm_lrt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_GWAS_lmm_score_f(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_GWAS_lmm_wald(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_GWAS_logit_wald_f(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_GWAS_logitmm_wald_f(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_hwe(SEXP, SEXP, SEXP);
extern SEXP gg_hwe_chi(SEXP, SEXP, SEXP);
extern SEXP gg_invert_snp_coding(SEXP, SEXP);
extern SEXP gg_Kinship_pw(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_Kinship_w(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_LD(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_LD_chunk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_LD_chunk_p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_LD_p(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_ld_thin_left(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_ld_thin_random(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_ld_thin_right(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_logp_thinning(SEXP, SEXP);
extern SEXP gg_m4_as_scaled_matrix_mu_sigma(SEXP, SEXP, SEXP);
extern SEXP gg_m4_as_scaled_matrix_p(SEXP, SEXP);
extern SEXP gg_m4_as012(SEXP);
extern SEXP gg_m4_loading_to_pc_ms(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_m4_loading_to_pc_p(SEXP, SEXP, SEXP);
extern SEXP gg_m4_pc_to_loading_ms(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_m4_pc_to_loading_p(SEXP, SEXP, SEXP);
extern SEXP gg_manhattan_thinning(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_ninds(SEXP);
extern SEXP gg_nsnps(SEXP);
extern SEXP gg_pre_likelihood(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_pre_likelihood_nofix(SEXP, SEXP, SEXP);
extern SEXP gg_random_ortho(SEXP);
extern SEXP gg_re_likelihood(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_re_likelihood_nofix(SEXP, SEXP, SEXP);
extern SEXP gg_read_bed_file(SEXP, SEXP, SEXP);
extern SEXP gg_read_vcf_filtered(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_read_vcf2(SEXP, SEXP, SEXP);
extern SEXP gg_set_snp_to_na(SEXP, SEXP);
extern SEXP gg_snp_hz_to_na(SEXP, SEXP);
extern SEXP gg_SNPmatch(SEXP, SEXP);
extern SEXP gg_which_duplicated_chr_pos(SEXP, SEXP);
extern SEXP gg_which_duplicated_chr_pos_alleles(SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_which_duplicated_id(SEXP);
extern SEXP gg_which_duplicated_id_chr_pos(SEXP, SEXP, SEXP);
extern SEXP gg_which_duplicated_id_chr_pos_alleles(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gg_write_bed_file(SEXP, SEXP);
extern SEXP isnullptr(SEXP);

static const R_CMethodDef CEntries[] = {
    {"qfc", (DL_FUNC) &qfc, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"gg_AIREML1",                             (DL_FUNC) &gg_AIREML1,                             15},
    {"gg_AIREML1_contrast",                    (DL_FUNC) &gg_AIREML1_contrast,                    15},
    {"gg_AIREML1_logit",                       (DL_FUNC) &gg_AIREML1_logit,                       14},
    {"gg_AIREML1_logit_nofix",                 (DL_FUNC) &gg_AIREML1_logit_nofix,                 11},
    {"gg_AIREML1_nofix",                       (DL_FUNC) &gg_AIREML1_nofix,                       14},
    {"gg_AIREMLn",                             (DL_FUNC) &gg_AIREMLn,                             15},
    {"gg_AIREMLn_contrast",                    (DL_FUNC) &gg_AIREMLn_contrast,                    15},
    {"gg_AIREMLn_logit",                       (DL_FUNC) &gg_AIREMLn_logit,                       14},
    {"gg_AIREMLn_logit_nofix",                 (DL_FUNC) &gg_AIREMLn_logit_nofix,                 11},
    {"gg_AIREMLn_nofix",                       (DL_FUNC) &gg_AIREMLn_nofix,                       14},
    {"gg_alleles_duplicated",                  (DL_FUNC) &gg_alleles_duplicated,                   2},
    {"gg_alleles_recoding",                    (DL_FUNC) &gg_alleles_recoding,                     1},
    {"gg_as_matrix4",                          (DL_FUNC) &gg_as_matrix4,                           1},
    {"gg_bind_inds2",                          (DL_FUNC) &gg_bind_inds2,                           2},
    {"gg_bind_snps",                           (DL_FUNC) &gg_bind_snps,                            1},
    {"gg_diago_full_likelihood1",              (DL_FUNC) &gg_diago_full_likelihood1,               6},
    {"gg_diago_full_likelihood1_nocovar",      (DL_FUNC) &gg_diago_full_likelihood1_nocovar,       5},
    {"gg_diago_full_likelihood2",              (DL_FUNC) &gg_diago_full_likelihood2,               7},
    {"gg_diago_full_likelihood2_nocovar",      (DL_FUNC) &gg_diago_full_likelihood2_nocovar,       6},
    {"gg_diago_likelihood1",                   (DL_FUNC) &gg_diago_likelihood1,                    6},
    {"gg_diago_likelihood1_nocovar",           (DL_FUNC) &gg_diago_likelihood1_nocovar,            5},
    {"gg_diago_likelihood2",                   (DL_FUNC) &gg_diago_likelihood2,                    7},
    {"gg_diago_likelihood2_nocovar",           (DL_FUNC) &gg_diago_likelihood2_nocovar,            6},
    {"gg_duplicated_remove",                   (DL_FUNC) &gg_duplicated_remove,                    7},
    {"gg_extract_inds_bool",                   (DL_FUNC) &gg_extract_inds_bool,                    2},
    {"gg_extract_inds_indices",                (DL_FUNC) &gg_extract_inds_indices,                 2},
    {"gg_extract_snps_bool",                   (DL_FUNC) &gg_extract_snps_bool,                    2},
    {"gg_extract_snps_indices",                (DL_FUNC) &gg_extract_snps_indices,                 2},
    {"gg_fit_diago",                           (DL_FUNC) &gg_fit_diago,                           10},
    {"gg_fit_diago_nocovar",                   (DL_FUNC) &gg_fit_diago_nocovar,                    9},
    {"gg_geno_stats",                          (DL_FUNC) &gg_geno_stats,                           5},
    {"gg_geno_stats_inds",                     (DL_FUNC) &gg_geno_stats_inds,                      4},
    {"gg_geno_stats_snps",                     (DL_FUNC) &gg_geno_stats_snps,                      3},
    {"gg_GWAS_lm_quanti",                      (DL_FUNC) &gg_GWAS_lm_quanti,                       6},
    {"gg_GWAS_lmm_lrt",                        (DL_FUNC) &gg_GWAS_lmm_lrt,                        10},
    {"gg_GWAS_lmm_score_f",                    (DL_FUNC) &gg_GWAS_lmm_score_f,                     6},
    {"gg_GWAS_lmm_wald",                       (DL_FUNC) &gg_GWAS_lmm_wald,                       10},
    {"gg_GWAS_logit_wald_f",                   (DL_FUNC) &gg_GWAS_logit_wald_f,                    7},
    {"gg_GWAS_logitmm_wald_f",                 (DL_FUNC) &gg_GWAS_logitmm_wald_f,                  8},
    {"gg_hwe",                                 (DL_FUNC) &gg_hwe,                                  3},
    {"gg_hwe_chi",                             (DL_FUNC) &gg_hwe_chi,                              3},
    {"gg_invert_snp_coding",                   (DL_FUNC) &gg_invert_snp_coding,                    2},
    {"gg_Kinship_pw",                          (DL_FUNC) &gg_Kinship_pw,                           5},
    {"gg_Kinship_w",                           (DL_FUNC) &gg_Kinship_w,                            5},
    {"gg_LD",                                  (DL_FUNC) &gg_LD,                                   5},
    {"gg_LD_chunk",                            (DL_FUNC) &gg_LD_chunk,                             7},
    {"gg_LD_chunk_p",                          (DL_FUNC) &gg_LD_chunk_p,                           6},
    {"gg_LD_p",                                (DL_FUNC) &gg_LD_p,                                 4},
    {"gg_ld_thin_left",                        (DL_FUNC) &gg_ld_thin_left,                        10},
    {"gg_ld_thin_random",                      (DL_FUNC) &gg_ld_thin_random,                      10},
    {"gg_ld_thin_right",                       (DL_FUNC) &gg_ld_thin_right,                       10},
    {"gg_logp_thinning",                       (DL_FUNC) &gg_logp_thinning,                        2},
    {"gg_m4_as_scaled_matrix_mu_sigma",        (DL_FUNC) &gg_m4_as_scaled_matrix_mu_sigma,         3},
    {"gg_m4_as_scaled_matrix_p",               (DL_FUNC) &gg_m4_as_scaled_matrix_p,                2},
    {"gg_m4_as012",                            (DL_FUNC) &gg_m4_as012,                             1},
    {"gg_m4_loading_to_pc_ms",                 (DL_FUNC) &gg_m4_loading_to_pc_ms,                  4},
    {"gg_m4_loading_to_pc_p",                  (DL_FUNC) &gg_m4_loading_to_pc_p,                   3},
    {"gg_m4_pc_to_loading_ms",                 (DL_FUNC) &gg_m4_pc_to_loading_ms,                  4},
    {"gg_m4_pc_to_loading_p",                  (DL_FUNC) &gg_m4_pc_to_loading_p,                   3},
    {"gg_manhattan_thinning",                  (DL_FUNC) &gg_manhattan_thinning,                   4},
    {"gg_ninds",                               (DL_FUNC) &gg_ninds,                                1},
    {"gg_nsnps",                               (DL_FUNC) &gg_nsnps,                                1},
    {"gg_pre_likelihood",                      (DL_FUNC) &gg_pre_likelihood,                       4},
    {"gg_pre_likelihood_nofix",                (DL_FUNC) &gg_pre_likelihood_nofix,                 3},
    {"gg_random_ortho",                        (DL_FUNC) &gg_random_ortho,                         1},
    {"gg_re_likelihood",                       (DL_FUNC) &gg_re_likelihood,                        4},
    {"gg_re_likelihood_nofix",                 (DL_FUNC) &gg_re_likelihood_nofix,                  3},
    {"gg_read_bed_file",                       (DL_FUNC) &gg_read_bed_file,                        3},
    {"gg_read_vcf_filtered",                   (DL_FUNC) &gg_read_vcf_filtered,                    4},
    {"gg_read_vcf2",                           (DL_FUNC) &gg_read_vcf2,                            3},
    {"gg_set_snp_to_na",                       (DL_FUNC) &gg_set_snp_to_na,                        2},
    {"gg_snp_hz_to_na",                        (DL_FUNC) &gg_snp_hz_to_na,                         2},
    {"gg_SNPmatch",                            (DL_FUNC) &gg_SNPmatch,                             2},
    {"gg_which_duplicated_chr_pos",            (DL_FUNC) &gg_which_duplicated_chr_pos,             2},
    {"gg_which_duplicated_chr_pos_alleles",    (DL_FUNC) &gg_which_duplicated_chr_pos_alleles,     4},
    {"gg_which_duplicated_id",                 (DL_FUNC) &gg_which_duplicated_id,                  1},
    {"gg_which_duplicated_id_chr_pos",         (DL_FUNC) &gg_which_duplicated_id_chr_pos,          3},
    {"gg_which_duplicated_id_chr_pos_alleles", (DL_FUNC) &gg_which_duplicated_id_chr_pos_alleles,  5},
    {"gg_write_bed_file",                      (DL_FUNC) &gg_write_bed_file,                       2},
    {"isnullptr",                              (DL_FUNC) &isnullptr,                               1},
    {NULL, NULL, 0}
};

void R_init_gaston(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



