
select_set <- function(set) {
if (set=="d3a") {
    samples <- list(list("moe_wt_hmc.bed", "moe_d3a_hmc.bed"), list("moe_wt_mc.bed", "moe_d3a_mc.bed"))
    legend <- c("Dnmt3a +/+", "Dnmt3a -/-")
    rows <- 2
    columns <- 1
  } else if (set=="d3a_2") {
    samples <- list(list("moe_d3a_wt_hmc", "moe_d3a_ko_hmc"), list("moe_d3a_wt_mc", "moe_d3a_ko_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_3") {
    samples <- list(list("d3a_wt_hmc", "d3a_ko_hmc"), list("d3a_wt_mc", "d3a_ko_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_4") {
    samples <- list(list("moe_d3a_wt_hmc_rpkm",  "moe_d3a_ko_hmc_rpkm"),
                    list("moe_d3a_wt_mc_rpkm", "moe_d3a_ko_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_hmc") {
    samples <- list(list("moe_d3a_wt_hmc_rpkm", "moe_d3a_ko_hmc_rpkm"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="d3a_mc") {
    samples <- list(list("moe_d3a_wt_mc_rpkm", "moe_d3a_ko_mc_rpkm"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="d3a_bt2") {
    samples <- list(list("moe_d3a_wt_hmc_bt2_rpkm",  "moe_d3a_ko_hmc_bt2_rpkm"),
                    list("moe_d3a_wt_mc_bt2_rpkm", "moe_d3a_ko_mc_bt2_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
    
  } else if (set=="d3xog_hmc") {
    samples <- list(list("omp_hmc_rep1_mean_omp_hmc_rep2", "d3xog_het_hmc_paired_q30", "d3xog_ko_hmc_paired_q30"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="d3xog_mc") {
    samples <- list(list("omp_mc_rmdup",  "d3xog_het_mc_paired_q30", "d3xog_ko_mc_paired_q30"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="d3xog_rmrna") {
    samples <- list(list("omp_rmrna_rep1_protein",
                         "omp_rmrna_rep2_protein",
                         "d3xog_wt_rmrna_blank_protein",  
                         "d3xog_het_rmrna_blank_protein", 
                         "d3xog_ko_rmrna_blank_protein",
                         "d3xog_ko_rmrna_rep2_blank_protein"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="d3a_rlm") {
    samples <- list(list("moe_wt_hmc_rlm", "moe_d3a_hmc_rlm"), list("moe_wt_mc_rlm", "moe_d3a_mc_rlm"))
    rows <- 2
    columns <- 1
  } else if (set=="d3a_rf") {
    samples <- list(list("moe_wt_mc_rf", "moe_d3a_mc_rlm"))
    rows <- 1
    columns <- 1
  } else if (set=="d3a_dnase") {
    samples <- list(list("d3a_het_dnase_sort_q30", "d3a_ko_dnase_sort_q30"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="cells") {
    samples <- list(list("omp_hmedip.bed", "ngn_hmedip.bed", "icam_hmedip.bed"),
                    list("omp_medip.bed", "ngn_medip.bed", "icam_medip.bed"))
    rows <- 2
    columns <- 1
    legend <- c("mOSN", "GBC", "HBC")
  } else if (set=="cells_rlm") {
    samples <- list(list("omp_hmc_rlm", "ngn_hmc_rlm", "icam_hmc_rlm"),
                    list("omp_mc_rlm", "ngn_mc_rlm", "icam_mc_rlm"))
    rows <- 2
    columns <- 1
  } else if (set=="cells_norm") {
    samples <- list(list("omp_hmc", "ngn_hmc", "icam_hmc"),
                    list("omp_mc", "ngn_mc", "icam_mc"))
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 2
    orient <- 2
  } else if (set=="cells_rpkm") {
    samples <- list(list("omp_hmc_120424_rpkm", "ngn_hmc_rpkm", "icam_hmc_rpkm"),
                    list("omp_mc_rpkm", "ngn_mc_rpkm", "icam_mc_rpkm"))
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 2
    orient <- 2
  } else if (set=="cells_hmc") {
    samples <- list(list("omp_hmc_120424_rpkm", "ngn_hmc_rpkm", "icam_hmc_rpkm"))
                    
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 1
    orient <- 2
    
  } else if (set=="cells_hmc_rep") {
    samples <- list(list("omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup", 
                         "ngn_hmc_rep1_q30_rmdup_extend300_mean_ngn_hmc_rep2_q30_rmdup",
                         "icam_hmc_rep1_q30_rmdup_extend300_mean_icam_hmc_rep2_q30_rmdup"))
    
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 1
    orient <- 2
    
  } else if (set=="cells_mc") {
    samples <- list(list("omp_mc_rpkm", "ngn_mc_rpkm", "icam_mc_rpkm"))
    
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 1
    orient <- 2
    
  } else if (set=="cells_norm_hmc") {
    samples <- list("omp_hmc", "ngn_hmc", "icam_hmc")
    columns <- 1
    rows <- 2
    orient <- 2
  } else if (set=="cells_norm_mc") {
    samples <- list("omp_mc", "ngn_mc", "icam_mc")
    columns <- 1
    rows <- 2
    orient <- 2
  } else if (set=="tfo_omp") {
    samples <- list(list("tfo_hmc", "omp_hmc"),
                    list("tfo_mc", "omp_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="tfo_omp_rpkm") {
    samples <- list(list("tfo_hmc_22M", "^omp_hmc_120424_rpkm"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="tfo_omp_ac") {
    samples <- list(list("tfo_hmc_22M", "^omp_hmc_120424_rpkm_ac"))
    rows <- 1
    columns <- 1
    orient <- 2

  } else if (set=="tfo_omp_d3a") {
    samples <- list(list("tfo_hmc", "omp_hmc", "moe_d3a_wt_hmc", "moe_d3a_ko_hmc"),
                    list("tfo_mc", "omp_mc", "moe_d3a_wt_mc", "moe_d3a_ko_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="nuc") {
    samples <- list(list("omp_nuc_0123", "icam_nuc_01234"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="tt3") {
    samples <- list(list("o.tt3.1_hmc_rpkm", "o.tt3.2_hmc_rpkm"),
                    list("o.tt3.1_mc_rpkm", "o.tt3.2_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="tt3_2") {
    samples <- list(list("omp_hmc_rpkm", "o.tt3.2_hmc_rpkm"),
                    list("omp_mc_rpkm", "o.tt3.2_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="tt3_3") {
    samples <- list(list("omp_hmc_120424_rpkm", "ott3_1_hmc_rpkm"),
                    list("omp_mc_rpkm", "ott3_1_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="tt3_3_hmc") {
    samples <- list(list("omp_hmc_120424_rpkm", "ott3_1_hmc_rpkm"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="tt3_3_mc") {
    samples <- list(list("omp_mc_rpkm", "ott3_1_mc_rpkm"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="tt3_rep") {
    samples <- list(list("omp_hmc_rep1_mean_omp_hmc_rep2", "ott3_hmc_rep1_mean_ott3_hmc_rep2"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="tt3_rmrna") {
    samples <- list(list("omp_rmrna_rep1", "omp_rmrna_rep2", 
                         "ott3_rmrna_rep1", "ott3_rmrna_rep2"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="d3xog_rmrna") {
    samples <- list(list("d3xog_wt_rmrna_blank_paired_q30", "d3xog_het_rmrna_blank_paired_q30", 
                         "d3xog_ko_rmrna_blank_paired_q30", "d3xog_ko_rmrna_rep2_blank_paired_q30"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="d3xog_rmrna_v_rmsk") {
    samples <- list(list("d3xog_wt_rmrna_blank_v_rmsk", "d3xog_het_rmrna_blank_v_rmsk", 
                         "d3xog_ko_rmrna_blank_v_rmsk", "d3xog_ko_rmrna_rep2_blank_v_rmsk"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="d3a_mrna") {
    samples <- list(list("moe_d3a_wt_mrna", "moe_d3a_ko_mrna"))
    rows <- 1
    columns <- 1
    orient <- 2
    
  } else if (set=="cells_nuc") {
    samples <- list(list("omp_nuc_0123","ngn_nuc_456", "icam_nuc_01234"))
    rows <- 1
    columns <- 1
    orient <- 2
  
  } else if (set=="d3a_nuc") {
    samples <- list(list("d3xog_wt_nuc_478_p1", "d3xog_ko_nuc_256_p1"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="d3a_nuc2") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup", "d3xog_ko_nuc_256_rmdup"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="d3a_nuc3") {
    samples <- list(list("d3xog_wt_nuc_478", "d3xog_ko_nuc_256"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="d3a_nuc3_sub_naked") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup_sub_omp_naked_s78", "d3xog_ko_nuc_256_rmdup_sub_omp_naked_s78"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="d3a_nuc_extend") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup_q30_extend", "d3xog_ko_nuc_256_rmdup_q30_extend"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="d3a_nuc_extend_sub") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup_q30_extend_insert", "d3xog_ko_nuc_256_rmdup_q30_506M_extend_insert"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="d3a_nuc_dyad") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup_q30_dyad", "d3xog_ko_nuc_256_rmdup_q30_dyad"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="d3a_nuc_dyad_1") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup_q30_dyad_1", "d3xog_ko_nuc_256_rmdup_q30_dyad_1"))
    rows <- 1
    columns <- 1
    orient <- 2 
  } else if (set=="encode_dnase") {
    samples <- list(list("d3a_het_dnase_sort_q30", "wgEncodeUwDnaseCerebrumC57bl6MAdult8wksAlnRep1", "wgEncodeUwDnaseCerebellumC57bl6MAdult8wksAlnRep1", "wgEncodeUwDnaseRetinaC57bl6MAdult1wksAlnRep1",
                    "wgEncodeUwDnaseHeartC57bl6MAdult8wksAlnRep1", "wgEncodeUwDnaseLiverC57bl6MAdult8wksAlnRep1"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="guanz_dorsal") {
    samples <- list(list("cd_hmc", "id_hmc"), list("cd_mc", "id_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="guanz_ventral") {
    samples <- list(list("cv_hmc", "iv_hmc"), list("cv_mc", "iv_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
    
  }

return(list(samples, rows, columns, orient))
}
