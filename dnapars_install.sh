#!/usr/bin/env bash

PREFIX=$CONDA_PREFIX
mkdir -p $CONDA_PREFIX/bin
export CC=gcc 

$CC -c tribal/dnapars/phylip.c -o tribal/dnapars/phylip.o -fcommon
$CC -c tribal/dnapars/seq.c -o tribal/dnapars/seq.o -fcommon
$CC -c tribal/dnapars/dnapars.c -o tribal/dnapars/dnapars.o -fcommon
$CC tribal/dnapars/seq.o tribal/dnapars/phylip.o tribal/dnapars/dnapars.o -lm -o $PREFIX/bin/dnapars

