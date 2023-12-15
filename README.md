## About
Revisit the [3D wave tank test cases (TCs)](https://github.com/EAGRE-water-wave-impact-modelling/3D-wave-tank-JCP2022) through the "VP approach" where the weak formulations are now generated automatically based on the time-discretised VP via Firedrake's `derivative()` functionality.

## Computation Codes
- Main file:
    - The new approach: `3D_tank_VP.py`
    - The old approach: `3D_tank.py` (to be updated)
- Shared files:
    - Settings of TCn: `settings_TCn.py`
    - Save the job details: `savings.py`

| Test Case    | New Approach | Old Approach (to be updated) |
| :---:        |    :----:    |     :----: |
| TC4      |   **`3D_tank_VP.py`** <br/>`settings_TC4.py`, `savings.py`<br/> folder `202002`  |  **`3D_tank.py`** <br/>`settings_TC4.py`, `savings.py`<br/> folder `202002`  |

## Post-Processing Codes
- TC1:
    - `pp_standing_wave_TC1VP_gif.py`
    - `pp_energy_TC1VP.py`
    - `pp_L2norn_TC1VP.py`
- TC4:
    - `pp_probes_TC4VP.py`
    - `pp_wavemaker_TC4VP.py`

## Developing Log
- :white_check_mark: Change the way by which the wavemaker-related functions update against time (2023.11.9)
- :white_check_mark: Output results into files
    - TC1: `energy.csv`, `.npy` binary files with field data, `readme.txt`
    - TC4: `energy.csv`, `probes.csv` with numerical measurements, `readme.txt`
- :white_check_mark: $\hat{\phi}(z)$ can be switched between 1 and high order Lagrange polynomial based on GLL points.
- :white_check_mark: Interpolation of the experimental data by using the `numpy` method `interp` instead of the `scipy` method `interp1d`.
- TO DO
    - Implement the above changes to `3D_tank.py`.
    
