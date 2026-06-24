#!/usr/bin/env bash
# =====================================================================
# run_c24med_rerun_DE.sh
#   Re-run Job D (ATN arms only, HL=5, ret1396) then Job E (f=1, HL=2.64, ret1396).
#   Run this AFTER Job F (hl5_chem1_ret1396) has finished.
#
#   Job D partial: non-ATN arms already complete in hl2p64_ret1396.rds (Job C).
#   This re-run overwrites hl5_ret1396.rds with ATN-arms-only results,
#   which is safe since the plotting script always reads non-ATN arms from
#   the Job C baseline file, never from Job D.
#
#   Jobs:
#     D-rerun  ATN arms only, f=0, HL=5, new ret  -> ..._hl5_ret1396.rds
#     E        pyr_atn,pyr_cfp_atn, f=1, HL=2.64  -> ..._hl2p64_chem1_ret1396.rds
# =====================================================================
set -u
cd /c/Users/ag4218/Local/GitHub/malariasimulation || exit 1

RETNEW=1395.69
OUT=dev/outputs
STAMP=$(date +%Y%m%d_%H%M%S)
MASTER=$OUT/c24med_rerun_DE_${STAMP}.log
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

log "=========== D-RERUN + E START ==========="

run_job "D-rerun ATN-only HL5" "$OUT/drv_ret_hl5_atn_rerun.log" "${PFX}_hl5_ret1396.rds" \
        ANTIMAL_HL_YEARS=5 ATN_CHEM_DOSE=0 \
        SWEEP_ARMS=atn,pyr_atn,pyr_cfp_atn,pyr_cfp_mc_atn_cd \
        OUT_SUFFIX=_ret1396 NET_RETENTION_DAYS=$RETNEW

run_job "E ret f1 HL2.64"      "$OUT/drv_ret_hl2p64_chem1_rerun.log" "${PFX}_hl2p64_chem1_ret1396.rds" \
        ANTIMAL_HL_YEARS=2.64 ATN_CHEM_DOSE=1 \
        SWEEP_ARMS=pyr_atn,pyr_cfp_atn \
        OUT_SUFFIX=_ret1396 NET_RETENTION_DAYS=$RETNEW

log "=========== D-RERUN + E COMPLETE ==========="
