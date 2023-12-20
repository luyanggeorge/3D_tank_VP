## About
Revisit the [3D wave tank test cases (TCs)](https://github.com/EAGRE-water-wave-impact-modelling/3D-wave-tank-JCP2022) through the "VP approach" where the weak formulations are now generated automatically based on the time-discretised VP via Firedrake's `derivative()` functionality. 

To compare the results with the old ones, some changes are also made to the original code, where manually derived weak formulations for the SE and SV schemes are explicitly used in *Firedrake*.

## Computation Codes
- Main file:
    - The new approach: `3D_tank_VP.py`
    - The old approach: `3D_tank.py`
- Shared files:
    - Settings of TCx: `settings_TCx.py`
    - Save the job details: `savings.py`
- Exclusive to the old approach:
    - Explicit weak formulations for SE and SV in `solvers_full.py`

| Test Case | New Approach | Old Approach |
| :---:     |    :----:    |   :----:     |
| TC1       |**`3D_tank_VP.py`** <br/>`settings_TC1.py`, `savings.py` | **`3D_tank.py`** <br/>`settings_TC1.py`, `savings.py`<br/>`solvers_full.py`  |
| TC4       |**`3D_tank_VP.py`** <br/>`settings_TC4.py`, `savings.py`<br/> folder `202002`  |  **`3D_tank.py`** <br/>`settings_TC4.py`, `savings.py` <br/> `solvers_full.py` <br/> folder `202002`  |

## Post-Processing Codes for both approaches
- TC1:
    - `pp_standing_wave_TC1VP_gif.py`
    - `pp_energy_TC1VP.py`
    - `pp_L2norn_TC1VP.py`
- TC4:
    - `pp_probes_TC4VP.py`
    - `pp_wavemaker_TC4VP.py`

## Developing Log
- :white_check_mark: Change the way by which the wavemaker-related functions update against time.
- :white_check_mark: Output results into files
    - TC1: `energy.csv`, a series of `.npy` binary files with field data, `readme.txt`
    - TC4: `energy.csv`, `probes.csv` with numerical measurements, `readme.txt`
- :white_check_mark: $\hat{\phi}(z)$ can be switched between 1 and high order Lagrange polynomial based on GLL points.
- :white_check_mark: Interpolation of the experimental data by using the *NumPy* method `interp` instead of the *SciPy* method [`interp1d`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html) as the latter might be removed in future *SciPy* versions.
- :white_check_mark: The above changes are also implemented into `3D_tank.py`. The Lagrange polynomial $\tilde{\varphi}_i(z)$ is now constructed based on GLL points.
    
