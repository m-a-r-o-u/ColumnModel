# The Atmospheric Column Model

This research model is intendet to investigate the role of radiation on diffusional growth.
Futhermore it served to evaluate the evolution of clouds.

## Usage

The following list of parameters may serve as a starting point, to get things running.
This work is preliminary and imcomplete, if you intend on building upon this work,
you have to dig into the source code.
Create a input file config.yaml:

```
model:
    t_max: 3000
    dt: 0.05
    grid:
        toa: 3000.
        gridlength: 25.
    initial_state:
        ALR: 0.004
        Theta0: 297.2
        p0: 100000.
        cloud_base: 500.
        cloud_roof: 500.
        w: 2.
    radiation:
        sw: false
        lw: false
        data_path: /home/m/Mares.Barekzai/phd/projects/column_model/data/afglus.dat
    particle_source:
        type: twomey
        N_sp: 200
    fluctuations:
        type: markov
        epsilon: 50.e-4
        l: 100.
    collisions:
        type: hall
    sedimentation:
        type: lookup
    advection:
        type: secondfirstorderupwind
        lifetime: 3000.
logger:
    type: netcdf
    dir_name: /project/meteo/scratch/Mares.Barekzai/phd/projects/column_model/
    file_name: dummy
```
