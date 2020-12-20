#!/bin/sh

strainphlan.py --ifn_samples /tmp/*.markers --output_dir /tmp --print_clades_only --nprocs_main 36 > /tmp/clades.txt

