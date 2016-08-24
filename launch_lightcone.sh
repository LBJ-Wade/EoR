#!/bin/bash
# Catherine Watkinson (Feb 2016).
# USAGE: bash launch_lightcone.sh <start_channel> <end_channel>

# where <start_channel> is the first frequency slice to be 'observed' in the lightcone (usually 0).
# and <end_channel> is the final frequency slice to be 'observed' in the lightcone.

# USAGE:
# You must provide three lightcones (one containing cosmological signal, one with foregrounds,
# and a third filled with zeros for noise to be added upon).

# IMPORTANT: OSKAR expects lightcones in units of Kelvin, with dimensions of
# constant FoV*FoV and frequency channel (running from 1 to n)

# ========================
start_channel=$1
end_channel=$2
# ========================

# USER: define foreground file string
fnamefg="Models/fg"

# USER: define cosmological signal file string
fnamecs="Models/cs_21cmFast_10deg"

# USER: define noise file string
fnameno="Models/noise"

echo "Running OSKAR on lightcones of following sky models:"
echo $fnamefg
echo $fnamecs
echo $fnameno

''' USER: Submit three parallel realisations of the OSKAR code with jobfile run_sim.tesla
    (this is for use on the Cambs HPC cluster and will be need to be ammended for other clusters)
    running on foregrounds, cosmological signal and noise (top to bottom).
     '''
sbatch submit_sim.tesla $fnamefg $start_channel $end_channel 1
sbatch submit_sim.tesla $fnamecs $start_channel $end_channel 2
sbatch submit_sim.tesla $fnameno $start_channel $end_channel 3
