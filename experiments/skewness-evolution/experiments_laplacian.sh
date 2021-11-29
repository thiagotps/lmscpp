#!/usr/bin/env bash
# LAPLACIAN !!!

set -euo pipefail

export NITER=1000
export N=3
export M=2
export SIGMAV2=0.00000000001
export EXP=6
export LANG=en_US
export BETA=0.005
export K=256
export PDF_START=0.0
export PDF_END=0.20
export PDF_SAMPLES=10000
export KERNEL_EXP=-3
export DIST="lap"

source ./experiment.sh
experiment
