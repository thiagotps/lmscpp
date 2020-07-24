#!/usr/bin/env bash

SKEWNESS_EXECUTABLE=$1
BASE_RESULTS_DIR=$2

foo()
{
    data_dir=$(mktemp -d --suffix=${DIST})

    sk_ia_file="${data_dir}/${DIST}/sk_ia.txt"
    sk_eea_file="${data_dir}/${DIST}/sk_eea.txt"
    sk_png_file="${data_dir}/${DIST}/sk.png"

    mse_ia_file="${data_dir}/${DIST}/mse_ia.txt"
    mse_eea_file="${data_dir}/${DIST}/mse_eea.txt"
    mse_png_file="${data_dir}/${DIST}/mse.png"

    sk_eea_cache="${data_dir}/${DIST}/sk_eea_cache_${N}_${M}_${DIST}.bin"
    mse_eea_cache="${data_dir}/${DIST}/mse_eea_cache_${N}_${M}_${DIST}.bin"

    mkdir -p ${data_dir}/${DIST}

    [[ ! -f "$mse_eea_cache" ]] && ${SKEWNESS_EXECUTABLE} --writecache "$mse_eea_cache" -N $N -M $M --indmode eea --outmode mse

    ${SKEWNESS_EXECUTABLE} -N $N -M $M --indmode ia -n $NITER -b $BETA --sv2 $SIGMAV2 --outmode mse -d $DIST  -o $mse_ia_file
    [[ $(diff "${BASE_RESULTS_DIR}/${DIST}/$(basename $mse_ia_file)" "$mse_ia_file") ]] \
        && echo "Failed on MSE IA ${DIST}" && exit 1

    ${SKEWNESS_EXECUTABLE} --readcache "$mse_eea_cache" -N $N -M $M --indmode eea -n $NITER -b $BETA --sv2 $SIGMAV2 --outmode mse -d $DIST -o $mse_eea_file
    [[ $(diff "${BASE_RESULTS_DIR}/${DIST}/$(basename ${mse_eea_file})" "$mse_eea_file") ]] \
        && echo "Failed on MSE EEA ${DIST}" && exit 1

    [[ ! -f "$sk_eea_cache" ]] && ${SKEWNESS_EXECUTABLE} --writecache "$sk_eea_cache" -N $N -M $M --indmode eea --outmode sk

    ${SKEWNESS_EXECUTABLE} -N $N -M $M --indmode ia -n $NITER -b $BETA --sv2 $SIGMAV2 --outmode sk  -d $DIST -o $sk_ia_file
    [[ $(diff "${BASE_RESULTS_DIR}/${DIST}/$(basename ${sk_ia_file})" "$sk_ia_file") ]] \
        && echo "Failed on SK IA ${DIST}" && exit 1

    ${SKEWNESS_EXECUTABLE} --readcache "$sk_eea_cache" -N $N -M $M --indmode eea -n $NITER -b $BETA --sv2 $SIGMAV2 --outmode sk -d $DIST -o $sk_eea_file
    [[ $(diff "${BASE_RESULTS_DIR}/${DIST}/$(basename ${sk_eea_file})" "$sk_eea_file") ]] \
        && echo "Failed on SK EEA ${DIST}" && exit 1

    rm -rf $data_dir
}

export NITER=1000
export N=3
export M=2
export SIGMAV2=0.00000000001
export LANG=en_US
export BETA=0.061
export DIST="gauss"

# Gaussian experiment
foo

export BETA=0.005
export DIST="lap"

# Laplacian experiment
foo
