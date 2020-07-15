#!/usr/bin/env bash

CLASSICAL_EXECUTABE=$1

for line in "1 1 1" "1 2 3" "1 3 19" "1 4 152" "1 5 1341" "2 1 5" "2 2 48"  "2 3 394"  "2 4 3517" "3 1  37"  "3 2 698"
do
    set -- $line
    L=$1
    M=$2
    NUM=$3

    echo "Testing configuration L=$L M=$M"

    RES=$(exec ${CLASSICAL_EXECUTABE}  -L $L -M $M | grep NUMBER_OF_EQUATIONS | cut -d ' ' -f 2)

    [[ "$NUM" != "$RES" ]] && echo "Expected $NUM equations, got $RES" && exit 1
done

exit 0
