#!/bin/bash

set -e

RESULTSDIR="results/summary"

mkdir -p $RESULTSDIR

declare -a nbs=(
#                "analysis_notebook.ipynb"
                "plot_growth_data.ipynb"
                )

for nb in "${nbs[@]}"
do
    echo "Running $nb"

    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=-1 \
        $nb

    echo "Converting $nb to Markdown"
    jupyter nbconvert \
        --output-dir $RESULTSDIR \
        --to markdown \
        $nb
done
