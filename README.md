## About
Revisit the [3D wave tank test cases (TCs)](https://github.com/EAGRE-water-wave-impact-modelling/3D-wave-tank-JCP2022) through the "VP approach" where the weak formulations are now generated automatically based on the time-discretised VP via Firedrake's `derivative()` functionality.

## Codes
- Main file:
    - The new approach: `3D_tank_VP.py`
    - The old approach: `3D_tank.py` (to be updated)
- Shared scripts:
    - Settings of TCn: `settings_TCn.py`
    - Save the job details: `savings.py`

## Developing Log
- :white_check_mark: Change the way by which the wavemaker-related functions update against time (2023.11.9)
- :white_check_mark: Output results into files
    - TC1: `energy.csv`, `.npy` binary files with field data, `readme.txt`
    - TC4: `energy.csv`, `probes.csv` with numerical measurements, `readme.txt`
- :white_check_mark: $\hat{\phi}(z)$ can be switched between 1 and high order Lagrange polynomial based on GLL points.
- :white_check_mark: Interpolation of the experimental data by using the `numpy` method `interp` instead of the `scipy` method `interp1d`.
- TO DO
    - Implement the above changes to `3D_tank.py`.
    
