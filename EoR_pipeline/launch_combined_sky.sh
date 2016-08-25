#!/bin/bash

# USAGE: bash launch_combined_sky.sh <start_channel> <end_channel> <simtype> <"filename">

# where
# <start/end_channel> are the first and last frequncy channel
# (an integer from 0 to max channel count) you want analysing.

# <simtype> is the type of simulation you are passing, usually this should be 3
# so that noise is added. However, there is the option to use this to analyse
# individual elements of the sky. 1 = foregrounds, 2 = cosmological signal.

# <"filename"> should be the filename string only (e.g. fname.fits -> "fname").

# IMPORTANT: OSKAR expects lightcones in units of Kelvin, with dimensions of
# constant FoV*FoV and frequency channel (running from 1 to n).

# Catherine Watkinson (Feb 2016)

# ========================
start_channel=$1
end_channel=$2
simtype=$3
# ========================

# USER: define combined cosmological signal and foregrounds file string
fname1=$4

fname_cs_fg_no=${fname1}

echo "Running OSKAR on a lightcone of the following sky model:"
echo ${fname_cs_fg_no}

sbatch submit_sim.tesla ${fname_cs_fg_no} ${start_channel} ${end_channel} ${simtype}
