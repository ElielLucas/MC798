#!/bin/bash

INSTANCIAS=("grf113.txt" "grf114.txt" "grf117.txt" "grf118.txt" "grf119.txt" "grf120.txt" "grf121.txt" "grf123.txt" "grf124.txt" "grf125.txt" "grf126.txt" "grf127.txt" "grf129.txt" "grf130.txt" "grf131.txt" "grf132.txt" "grf133.txt" "grf134.txt")
SEED=42
ROTULO=295745
MAXTIME=600

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
