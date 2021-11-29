#!/usr/bin/env bash

set -euo pipefail

export NITER=175
export N=2
export M=3
export SIGMAV2=0.001
export EXP=6
export LANG=en_US
export DIST="gauss"
export DIR_PREFIX="data_${DIST}"

foo()
{
    local beta=$1
    local mse_eea_cache="mse_eea_cache_${N}_${M}.bin"
    local DIR="${DIR_PREFIX}_${N}_${M}/$beta"
    local mse_ia_file="${DIR}/mse_ia_${N}_${M}.txt"
    local mse_st_file="${DIR}/mse_st_${N}_${M}.txt"
    local mse_eea_file="${DIR}/mse_eea_${N}_${M}.txt"
    local mse_png_file="${DIR}/mse_${N}_${M}.png"
    local mse_eps_file="${DIR}/mse_${N}_${M}.eps"

    mkdir -p "$DIR"

   # Por questões de desempenho, guardamos todas as recursões em um arquivo.
   # Se esse arquivo não existir, criamos ele.
   [[ ! -f "$mse_eea_cache" ]] && skewness --writecache "$mse_eea_cache" -N $N -M $M --indmode eea --outmode mse

   # Obtemos o MSE usando EEA
   skewness --readcache "$mse_eea_cache" -N "$N" -M "$M" --indmode eea -n "$NITER" -b "$beta" --sv2 \
   "$SIGMAV2" --outmode mse -d "$DIST" -o "$mse_eea_file"

   # Obtemos o MSE usando IA
   skewness  -N "$N" -M "$M" --indmode ia -n "$NITER" -b "$beta" --sv2 "$SIGMAV2" --outmode mse -d \
   "$DIST" -o "$mse_ia_file"

   # Obtemos o MSE usando ensaios de Monte Carlo
   skemp --filter-length $N --data-length $M --niter $NITER \
       --exp "$EXP" --beta "$beta" --sigmav2 "$SIGMAV2" --dist "$DIST" \
       --mse-file "$mse_st_file" \

   echo "Creating the graphs"

   # Finalmente criamos os gráficos
   sgg --yt db  --xmin 0 --xmax "$NITER" -y "MSE (dB)" -s ":" " -" " --" -l "IA" "EEA" "Empírico" \
       -f  "$mse_ia_file" "$mse_eea_file" "$mse_st_file" -d "$mse_png_file" --eps-file "$mse_eps_file" \
       -c "lightseagreen" "darkorange" "royalblue"
}

foo "0.04"
