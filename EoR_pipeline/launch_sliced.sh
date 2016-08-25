#!/bin/bash
# Catherine Watkinson (Feb 2016) edit: Added comments.
# I THINK THIS IS REDUNDANT, WHY PRESLICE WHEN OSKAR WILL DO THE SLICING for
# YOU. THIS IS A USEFUL APPROACH IN PRINCIPAL AS IT ALLOWS DIFFERENT CHANNELS
# TO BE RUN IN PARALLEL.

# Runs the code on a slice-by-slice basis.

# Your three lightcones need to have been presliced into frequency channels
# each one labelled with channel numbers), e.g. for channel 5: fname5.fits
# Provide sky models for the cosmological signal, foregrounds,
# and a third (filled with zeros) for noise to be added upon.

# IMPORTANT: OSKAR expects lightcones in units of Kelvin, with dimensions of
# constant FoV*FoV and frequency channel (running from 1 to n).

# ========================
start_channel=$1
end_channel=$2
# ========================

for i in $(eval echo {${start_channel}..${end_channel}}); do
    echo "Channel ${i}"
    # USER: define foreground file string
    fname1="Models/fg_zeromean_highres"
    # USER: define cosmological signal file string
    fname4="Models/cs_sf_tapered_K_zeromean_highres"
    # USER: define noise file string
    fname44="Models/noise"

    printf -v fname2 "%03d" ${i}
    #cw(redundant) fname3=".fits"

    fnamefg=${fname1}${fname2}
    fnamecs=${fname4}${fname2}
    fnameno=${fname44}${fname2}

    echo "Running OSKAR on a slice-by-slice basis to the following sky models:"
    echo ${fnamefg}
    echo ${fnamecs}
    echo ${fnameno}

    #USER: Submit three parallel realisations of the OSKAR code with jobfile
    #      run_sim.tesla (this is for use on the Cambs HPC cluster and will be
    #      need to be ammended for other clusters) running on foregrounds,
    #      cosmological signal and noise (top to bottom).
    sbatch submit_sim.tesla ${fnamefg} ${i} ${i} 1
    sbatch submit_sim.tesla ${fnamecs} ${i} ${i} 2
    sbatch submit_sim.tesla ${fnameno} ${i} ${i} 3
done
