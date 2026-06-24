#!/usr/bin/env bash
# =====================================================================
# run_c24med_overnight.sh
#   Sequentially runs 6 c24med projection jobs, each starting only after
#   the previous one's .rds has been saved.  Every job writes a DISTINCT
#   output file (verified no-overwrite).  Continues to the next job even
#   if one returns a non-zero exit (e.g. cosmetic trailing error after
#   the save), logging a clear PASS/FAIL based on whether the .rds exists.
#
#   Jobs (in order):
#     Q1-A  mixed arm, HL=2.64, site retention   -> ..._hl2p64_mixonly.rds
#     Q1-B  mixed arm, HL=5,    site retention   -> ..._hl5_mixonly.rds
#     C     all 8 arms, f=0, HL=2.64, new ret    -> ..._hl2p64_ret1396.rds
#     D     all 8 arms, f=0, HL=5,    new ret    -> ..._hl5_ret1396.rds
#     E     pyr_atn,pyr_cfp_atn, f=1, HL=2.64, new ret -> ..._hl2p64_chem1_ret1396.rds
#     F     pyr_atn,pyr_cfp_atn, f=1, HL=5,    new ret -> ..._hl5_chem1_ret1396.rds
#
#   New retention = 2013.56 * log(2) = 1395.69 d (site value treated as half-life).
# =====================================================================
set -u
cd /c/Users/ag4218/Local/GitHub/malariasimulation || exit 1

RETNEW=1395.69
OUT=dev/outputs
STAMP=$(date +%Y%m%d_%H%M%S)
MASTER=$OUT/c24med_overnight_driver_${STAMP}.log
PFX=mli_c24med_projection_results

log() { echo "[$(date '+%F %T')] $*" | tee -a "$MASTER"; }

run_job () {
  local tag="$1"; local logf="$2"; local expect="$3"; shift 3
  log "START  $tag"
  log "       env: $*"
  env "$@" SWEEP_CORES=18 \
    Rscript -e "devtools::load_all('.', quiet=TRUE); source('dev/c24med_projection_run.R')" \
    > "$logf" 2>&1
  local rc=$?
  if [ -f "$OUT/$expect" ]; then
    local sz; sz=$(du -h "$OUT/$expect" | cut -f1)
    log "DONE   $tag  (exit $rc)  PASS rds saved: $expect ($sz)"
  else
    log "DONE   $tag  (exit $rc)  *** FAIL: expected rds missing: $expect *** (see $logf)"
  fi
}

log "=========== OVERNIGHT DRIVER START (6 jobs) ==========="

run_job "Q1-A mixed HL2.64"  "$OUT/drv_mixed_hl2p64.log"     "${PFX}_hl2p64_mixonly.rds" \
        ANTIMAL_HL_YEARS=2.64 ATN_CHEM_DOSE=0 SWEEP_ARMS=pyr_cfp_mc_atn_cd OUT_SUFFIX=_mixonly

run_job "Q1-B mixed HL5"     "$OUT/drv_mixed_hl5.log"        "${PFX}_hl5_mixonly.rds" \
        ANTIMAL_HL_YEARS=5    ATN_CHEM_DOSE=0 SWEEP_ARMS=pyr_cfp_mc_atn_cd OUT_SUFFIX=_mixonly

run_job "C ret f0 HL2.64"    "$OUT/drv_ret_hl2p64.log"       "${PFX}_hl2p64_ret1396.rds" \
        ANTIMAL_HL_YEARS=2.64 ATN_CHEM_DOSE=0 OUT_SUFFIX=_ret1396 NET_RETENTION_DAYS=$RETNEW

run_job "D ret f0 HL5"       "$OUT/drv_ret_hl5.log"          "${PFX}_hl5_ret1396.rds" \
        ANTIMAL_HL_YEARS=5    ATN_CHEM_DOSE=0 OUT_SUFFIX=_ret1396 NET_RETENTION_DAYS=$RETNEW

run_job "E ret f1 HL2.64"    "$OUT/drv_ret_hl2p64_chem1.log" "${PFX}_hl2p64_chem1_ret1396.rds" \
        ANTIMAL_HL_YEARS=2.64 ATN_CHEM_DOSE=1 SWEEP_ARMS=pyr_atn,pyr_cfp_atn OUT_SUFFIX=_ret1396 NET_RETENTION_DAYS=$RETNEW

run_job "F ret f1 HL5"       "$OUT/drv_ret_hl5_chem1.log"    "${PFX}_hl5_chem1_ret1396.rds" \
        ANTIMAL_HL_YEARS=5    ATN_CHEM_DOSE=1 SWEEP_ARMS=pyr_atn,pyr_cfp_atn OUT_SUFFIX=_ret1396 NET_RETENTION_DAYS=$RETNEW

log "=========== OVERNIGHT DRIVER COMPLETE ==========="
