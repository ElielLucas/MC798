#!/bin/bash

INSTANCIAS=("grf1.txt" "grf2.txt" "grf5.txt" "grf6.txt" "grf19.txt" "grf101.txt" "grf102.txt" "grf103.txt" "grf104.txt" "grf105.txt" "grf106.txt" "grf107.txt" "grf108.txt" "grf109.txt" "grf110.txt" "grf111.txt" "grf112.txt")
SEED=42
ROTULO=295745
MAXTIME=1500

for INST in "${INSTANCIAS[@]}"; do
  LOGFILE="log-${INST%.txt}.txt"
  OUTFILE="saida-${INST%.txt}.txt"

  echo "Rodando instância: $INST" | tee "$OUTFILE"
  
  julia trigger_arc_tsp.jl \
    -f "$INST" \
    -s "$SEED" \
    --maxtime_lb_lp $MAXTIME \
    --maxtime_lb_rlxlag $MAXTIME \
    --maxtime_lb_colgen $MAXTIME \
    --maxtime_ub_lp $MAXTIME \
    --maxtime_ub_rlxlag $MAXTIME \
    --maxtime_ub_colgen $MAXTIME \
    --maxtime_ilp $MAXTIME \
    -r "$ROTULO" \
    -l "$LOGFILE" \
    2>&1 | tee -a "$OUTFILE"
done
