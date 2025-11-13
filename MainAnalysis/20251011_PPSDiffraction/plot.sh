#!/usr/bin/env bash
set -e

# Compile loop.C into an executable called 'loop'
g++ -O2 -std=c++17 plot.C $(root-config --cflags --libs) -o plot

# Run it (forward any args)
./plot "$@"
