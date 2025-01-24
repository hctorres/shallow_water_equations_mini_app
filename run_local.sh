#!/bin/bash

set -e

# Build
make

# Run
./main2d.gnu.ex inputs

# Postprocess
python plot_with_yt.py