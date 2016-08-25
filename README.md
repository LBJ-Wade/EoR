# Scripts for EoR simulations with the SKA1-low

## 1. Noise cubes
Scripts for simulating noise cubes can be found in the folder 
`noise_cube/`

- `sim_cube.py`: Script to generate noise cubes (simulation + imaging). 
    Simulates noise for a 1000h observation equivalent to repeatedly 
    observing the same field for 5 hours for the same range of hour 
    angles. Uncorrelated random Gaussian noise is added to Stokes-I 
    visibilities with a noise RMS based on antenna element system 
    temperature and effective area values provided by Eloy de Lera 
    Acedo (U. Cambridge) for simulated SKALA antennas. Note these
    values do not take into account of scan angle. Visibilities are
    imaged using a wavelength dependent mask (lambda cut) aimed at
    making the PSF between frequency channels in each of the 8MHz
    cubes as consistent as possible. Imaging is currently performed
    using the OSKAR imager with uniform and natural weighting.

- `plot_cube.py`: Script to plot data the noise cube simulation 
   (used for debugging)

See https://www.overleaf.com/read/nqqdzzdqnpqc for a description of the current noise cube simulation.


## 2. EoR pipeline scripts

EoR pipeline scipts can be found in the folder `EoR_pipeline/`.
