#!/usr/bin/env bash

CLASSICAL_EXECUTABE=$1
BASE_RESULTS_DIR=$2

for line in "1 1 1" "1 2 3" "1 3 19" "1 4 152" "1 5 1341" "2 1 5" "2 2 48"  "2 3 394"  "2 4 3517" "3 1  37"  "3 2 698"
do
    set -- $line
    L=$1
    M=$2
    NUM=$3

    echo "Testing configuration L=$L M=$M"

    RES=$(${CLASSICAL_EXECUTABE}  -L $L -M $M | grep NUMBER_OF_EQUATIONS | cut -d ' ' -f 2)

    [[ "$NUM" != "$RES" ]] && echo "Expected $NUM equations, got $RES" && exit 1
done

nummatrix_file=$(mktemp)
symmatrix_file=$(mktemp)
evol_file=$(mktemp)

${CLASSICAL_EXECUTABE} --writecache -L 1 -M 2 --nummatrix "$nummatrix_file" --symmatrix "$symmatrix_file" --evolution "$evol_file" --niter 1000

[[ $(diff "${BASE_RESULTS_DIR}/12basenummatrix.txt" "$nummatrix_file") ]] && echo "The actual numerical matrix and the base are not equal." && exit 1
[[ $(diff "${BASE_RESULTS_DIR}/12basesymmatrix.txt" "$symmatrix_file") ]] && echo "The actual symbolical matrix and the base are not equal." && exit 1
[[ $(diff "${BASE_RESULTS_DIR}/12baseevol.txt" "$evol_file") ]] && echo "The actual matrix's evolution and the base are not equal." && exit 1

rm -f "$nummatrix_file" "$symmatrix_file" "$evol_file"
${CLASSICAL_EXECUTABE} --readcache -L 1 -M 2 --nummatrix "$nummatrix_file" --symmatrix "$symmatrix_file" --evolution "$evol_file" --niter 1000

[[ $(diff "${BASE_RESULTS_DIR}/12basenummatrix.txt" "$nummatrix_file") ]] && echo "The numerical matrix obtained through the cache differ." && exit 1
[[ $(diff "${BASE_RESULTS_DIR}/12basesymmatrix.txt" "$symmatrix_file") ]] && echo "The symbolical matrix obtained through the cache differ." && exit 1
[[ $(diff "${BASE_RESULTS_DIR}/12baseevol.txt" "$evol_file") ]] && echo "The matrix's evolution obtained through the cache differ." && exit 1

rm -f "$nummatrix_file" "$symmatrix_file" "$evol_file"

exit 0
