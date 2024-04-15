#!/bin/bash

if [[ $# -ne 1 ]]; then
    echo "Illegal number of arguments" >&2
    echo "Must provide the sampling duration (e.g. 10s)"
    exit 2
fi

echo "Sampling for $1"
sudo Scripts/disco.d $1

graphname="out.disco.perf.svg"
echo "Writing $graphname"
stackcollapse.pl out.disco.stacks | flamegraph.pl > $graphname

