#!/usr/bin/env bash

set -euo pipefail

print_options()
{
    echo "DIST=$DIST"
    echo "NITER=$NITER"
    echo "N=$N"
    echo "M=$M"
    echo "SIGMAV2=$SIGMAV2"
    echo "EXP=$EXP"
    echo "LANG=$LANG"
    echo "BETA=$BETA"
    echo "K=$K"
    echo "PDF_START=$PDF_START"
    echo "PDF_END=$PDF_END"
    echo "PDF_SAMPLES=$PDF_SAMPLES"
    echo "KERNEL_EXP=$KERNEL_EXP"
}

experiment()
{
    DIR="data/${DIST}"
    sk_ia_file="${DIR}/sk_ia.txt"
    sk_eea_file="${DIR}/sk_eea.txt"
    sk_st_file="${DIR}/sk_st.txt"
    sk_eps_file="${DIR}/sk.eps"

    mse_ia_file="${DIR}/mse_ia.txt"
    mse_eea_file="${DIR}/mse_eea.txt"
    mse_st_file="${DIR}/mse_st.txt"
    mse_eps_file="${DIR}/mse.eps"

    pdf_file="${DIR}/pdf_st.txt"
    pdf_gauss_file="${DIR}/pdf_gauss.txt"
    pdf_eps_file="${DIR}/pdf.eps"

    sk_eea_cache="sk_eea_cache_${N}_${M}.bin"
    mse_eea_cache="mse_eea_cache_${N}_${M}.bin"

    mkdir -p "$DIR"

    print_options > "${DIR}/README"

    [[ ! -f "$sk_eea_cache" ]] && skewness --writecache "$sk_eea_cache" -N $N -M $M --indmode eea --outmode sk
    skewness --readcache "$sk_eea_cache" -N "$N" -M "$M" --indmode eea -n "$NITER" -b "$BETA" --sv2 "$SIGMAV2" --outmode sk -d "$DIST" -o "$sk_eea_file"
    skewness  -N "$N" -M "$M" --indmode ia -n "$NITER" -b "$BETA" --sv2 "$SIGMAV2" --outmode sk -d "$DIST" -o "$sk_ia_file"


    [[ ! -f "$mse_eea_cache" ]] && skewness --writecache "$mse_eea_cache" -N $N -M $M --indmode eea --outmode mse
    skewness --readcache "$mse_eea_cache" -N "$N" -M "$M" --indmode eea -n "$NITER" -b "$BETA" --sv2 "$SIGMAV2" --outmode mse -d "$DIST" -o "$mse_eea_file"
    skewness  -N "$N" -M "$M" --indmode ia -n "$NITER" -b "$BETA" --sv2 "$SIGMAV2" --outmode mse -d "$DIST" -o "$mse_ia_file"


    skemp --filter-length "$N" --data-length "$M" --niter "$NITER" \
        --exp "$EXP" --beta "$BETA" --sigmav2 "$SIGMAV2" --skewness-file "$sk_st_file" \
        --mse-file "$mse_st_file" \
        --pdf-file "$pdf_file" --pdf-instant "$K" --pdf-start "$PDF_START" --pdf-end "$PDF_END" --kernel-exp "$KERNEL_EXP" \
        --pdf-samples "$PDF_SAMPLES" \
        --dist "$DIST"

    mean_sd=$(awk -F '#'  'FNR == 2 {print $2}'  $pdf_file)
    mean=$(echo "$mean_sd" | cut -d ',' -f1)
    sd=$(echo "$mean_sd" | cut -d ',' -f2)

    printf "mean = %s, σ = %s\n" "$mean" "$sd"

    gauss_emp --mean "$mean" --sd "$sd" --start "$PDF_START" --end "$PDF_END" --samples "$PDF_SAMPLES" \
        --exp 6 --kernel-exp "$KERNEL_EXP" > "$pdf_gauss_file"

    sgg -y "Skewness" -s ":" " -" " --" -l "IA" "EEA" "Empírico" \
        -f  "$sk_ia_file" "$sk_eea_file" "$sk_st_file" --eps-file "$sk_eps_file" \
        -c "lightseagreen" "darkorange" "royalblue"

    sgg -y "MSE" -s ":" " -" " --" -l "IA" "EEA" "Empírico" \
        -f  "$mse_ia_file" "$mse_eea_file" "$mse_st_file" --eps-file "$mse_eps_file" \
        -c "lightseagreen" "darkorange" "royalblue" --yt "db" --yt-files "$mse_ia_file" "$mse_eea_file" "$mse_st_file"

    sgg -y "Função Densidade de Probabilidade" -s " -" " --" -l "Gaussiano" "Empírico" \
        -f   $pdf_gauss_file $pdf_file --eps-file $pdf_eps_file \
        -c  "darkorange" "royalblue"

}
