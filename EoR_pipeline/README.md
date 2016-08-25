README: EoR pipeline
--------------------

###### To think about:
Would be nice to have all called programs in a separate function so the choices
appear less overwhelming. This directory would then just consist of launch 
scripts, cleanup and possibly the lightcone generation (although it's likely 
that this will have a launch script associated with it).

Should the name of the visibilities folder (and maybe others) be made a user
variable (the advantage being that you could have several different runs in 
parallel, this would need to source an independent telescope model folder as 
well. Could be a nice simple way of being able to have several different 
observations running at once.)

-------------------------------

In order to run the EoR pipeline you need to install:
`python`, `numpy`, `pillow`, `pyfits` and `OSKAR`.

For more information on OSKAR's installation please visit:
http://www.oerc.ox.ac.uk/~ska/oskar2/

If running outside of the Cambridge HPC, then jobfiles will need to be
adjusted for qsubbing both OSKAR and CASA;
`submit_sim.tesla` and `sumbit_make_images.sandybridge` can be used as guides.

`LD_LIBRARY_PATH` in `sumbit_make_images` jobfile must also be directed to the
location of your version of the EoR pipeline’s `User_defined_variables/` folder,
where all the USER_DEFINED runtime variables are set.

There are several steps to running this pipeline:


> -------------------------------------------------------------------------
> --------- Nb. runtime user defined variables can be adjusted in ---------  
> ---------        User_defined_variables/USER_DEFINED.py         ---------
> -------------------------------------------------------------------------


1. produce a lightcone from a series of simulated cubes (see __??__).
NOTE - The lightcone should have units of Kelvin and the data type should be 
float32.  
—> NEED TO GET/WRITE A LIGHTCONE GENERATING SCRIPT TO CONVERT BOX IN CPC TO
FoF*FoV*FullBandwidthOfObservation. Provide ability to convert as simulation
from mK to K in the process (useful as e.g. 21cmFAST outputs in mK).

2. run ___??__ to get corresponding foreground simulated lightcone 
(see __??__)  
—> NEED TO GET FOREGROUND GENERATION CODE FROM EMMA AND FIND OUT (THEN DOCUMENT)
HOW TO USE IT. IDEALLY WOULD GET HANDS ON A FULL-SKY MODEL FOR THE FOREGROUNDS
AND JUST BUNDLE IT. VIBOR WILL PROVIDE THIS ONCE WE DECIDE ON THE APPROPRIATE
RESOLUTION ETC.

3. run one of the launch_____.sh scripts to get a visibility measurement set 
from OSKAR.

4. turn the measurement set into an image using the 
`submit_make_images.sandybridge` (Nb. If not being used on the Cambs HPC 
cluster, then this will have to be adjusted for the relevant cluster.)  
—> HAVE ADJUSTED THIS SCRIPT TO RUN THE COLLECT IMAGES .sh AFTER IMAGE 
GENERATION. WOULD BE GOOD TO ALSO INCORPORATE A SCRIPT THAT WILL SUM OVER THE 
SEPARATE COMPONENTS OF AN OBSERVATION (i.e. cosmo signal, noise, foregrounds) 
TO CREATE THE COMBINED SKY IMAGE.

5. (optional) Remove foregrounds using GMCA (the alternatives being to use a
different foreground removal algorithm or to exploit the foreground window and
utilise use foreground avoidance)  
--> WOULD NEED TO INTEGRATE GMCA AND/OR SIMILAR, CAN ALREADY BUNDLE 21CMSENSE 
TO PERFORM ROUGH SENSITIVITY CALCULATIONS, WHICH CAN MODEL FOREGROUND AVOIDANCE
OR PERFECT FOREGROUND REMOVAL. I SUGGEST THAT FOR ALL OF THESE EXTERNAL CODES 
THEY ARE PULLED FROM THE RESPECTIVE GITHUBS BY A SCRIPT. WOULD REQUIRE EMMA TO 
MAKE HER GMCA (AND ANY OTHERS SHE HAS) AVAILABLE PUBLICLY.

## SCRIPTS OVERVIEW

### Running OSKAR (to generate a visibility measurement set)

`launch_lightcone.sh`

USAGE: `bash launch_lightcone.sh <start_channel> <end_channel>` 

where `<start/end_channel>` are the first and last frequency channel
(an integer from 0 to max channel count) you want analysing.
The prefix of the filenames under analysis needs to be manually adjusted within
`launch_lightcone.sh`

Shell script to run OSKAR on user defined sky models (e.g. cosmological signal
and foregrounds) independently.
Nb. The output of these runs must ultimately be combined to give the total 
observed sky.

------------------------------------------
`launch_sliced.sh`

USAGE: `bash launch_sliced.sh <start_channel> <end_channel>`  

where `<channel>` is the number of frequency slices in the lightcone.

As above, but simulations have been pre-split into individual frequency slices,
labelled according to slice number (e.g. cs_21cmFast_10deg05) describes
the slice of frequency channel 5 of the cs_21cmFast_10deg lightcone.
The prefix of the filenames under analysis needs to be manually adjusted within
launch_sliced.sh

This is a useful mode as each frequency slice is processed in parallel. However,
it seems unnecessary to go to all the bother of presplitting the input data when
runOSKAR.py can just extract a particular slice. This is presumably useful
if you want to fix the PSF to that of the lowest frequency for all frequencies 
as per Emma FG removal papers.

------------------------------------------
`launch_combined_sky.sh`

USAGE: `bash launch_combined_sky.sh <start_channel> <end_channel> <simtype> 
                                    <"filename">`

`<"filename">` should be the filename string only (e.g. fname.fits -> "fname").
where `<start/end_channel>` are the first and last frequency channel
(an integer from 0 to max channel count) you want analysing.

`<simtype>` is the type of simulation you are passing, usually this should be 3
so that noise is added. However, there is the option to use this to analyse
individual elements of the sky. 1 = foregrounds, 2 = cosmological signal.

This is the most basic way to run OSKAR. It takes a sky model containing both cs
and fg (filename), adds noise and runs OSKAR to produce
visibility set and images for the combined sky.

### Running CASA (returning images from the visibility set)

`submit_make_images.sandybridge`

USAGE: `qsub submit_make_images.sandybridge <start_channel> <end_channel>`  
where `<start_channel>` and `<end_channel>` are the first and last channels 
you wish to image.

This generates images from visibility measurement sets using make_images.py and
then calls collect_images.sh to organise the produced images.

------------------------------------------
`make_images.py`

This converts the visibility set data to images and its running is managed by
submit_make_images.sandybridge

------------------------------------------
`collect_images.sh`

This organises the images produced by make_images into foreground, cosmological
signal and noise folders. Also run by submit_make_images.sandybridge.

### Auxilliary scripts

`cleanup.py`

USAGE: `bash cleanup.py`

This basically removes all output from an OSKAR and make_images run.

**!!!!--------------------- IMPORTANT ---------------------------!!!!**  
Make sure you have moved anything you don't want to lose from the
EoR pipeline folders before using.
