#!/usr/bin/env bash
# GAUSSIAN !!!

set -euo pipefail

export NITER=120
export N=2
export M=2
export SIGMAV2=0.001
export EXP=11
export LANG=en_US
export BETA=0.09
export K=16
export PDF_START=-0.5
export PDF_END=0.5
export PDF_SAMPLES=1000
export KERNEL_EXP=-3
export DIST="gauss"

source ./experiment.sh

experiment
