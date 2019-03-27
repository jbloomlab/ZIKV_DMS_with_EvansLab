#!/bin/bash

# name of analysis notebook
NB="analysis_nb.ipynb"

# execute the analysis notebook
jupyter nbconvert \
    --to notebook \
    --execute \
    --inplace \
    --ExecutePreprocessor.timeout=-1 \
    $NB

# create markdown notebook output in `results/summary/`
SUMMARYDIR="$RESULTSDIR/summary"
mkdir -p $SUMMARYDIR
jupyter nbconvert \
    --to markdown \
    --output-dir=$SUMMARYDIR \
    $NB

