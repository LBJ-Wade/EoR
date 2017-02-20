#!/bin/bash
# sbatch launch_sims_tidy.sh 140.0 140.5 0.5 'LOFAR' 'NCP' 1200 10 'HBA_CS.tm' 'NCP_below_1Jy' 2
# mode determines whether .osm is needed etc... for diffuse simulations mode 1, pre-existing .osm mode 2, noise mode 3
# ========================
f1=$1
f2=$2
telescope=$4
echo $telescope
field=$5
pix=$6
fov=$7
df=$3
telescope_model=$8
file_name=$9
mode=${10}
echo "mode $mode"
start_channel=0
end_channel_float=$(echo "scale=1;($f2-$f1)/$df" | bc)
end_channel=${end_channel_float%.*} 
echo "start $start_channel"
echo "end $end_channel"
# ========================

for i in $(eval echo {$start_channel..$end_channel}); do
    echo "Channel $i"
    freq=$(echo "scale=1;$f1+$df*$i" | bc)
    echo $freq
    printf -v fname2 "_%05.1fMHz" $freq
    fname3=".fits"

    fname=$file_name$fname2

    echo $fname
    echo $mode

     sbatch run_sim_tidy.tesla $freq $mode $fname $fov $pix $telescope $field $telescope_model
done

