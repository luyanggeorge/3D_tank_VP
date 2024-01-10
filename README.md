## About
Revisit the [3D wave tank test cases (TCs)](https://github.com/EAGRE-water-wave-impact-modelling/3D-wave-tank-JCP2022) through the "VP approach" where the weak formulations are now generated automatically based on the time-discretised VP via Firedrake's `derivative()` functionality. 

To compare the results with the old ones, some changes are also made to the original code, where manually derived weak formulations for the SE and SV schemes are explicitly implemented in *Firedrake*.

## Computation Codes
- Main file:
    - The new approach: `3D_tank_VP.py`
    - The old approach: `3D_tank.py`
- Shared files:
    - Settings of TCx: `settings_TCx.py`
    - Save the job details: `savings.py`
- Exclusive to the old approach:
    - Explicit weak formulations for SE and SV in `solvers_full.py`

## Developing Notes
- Change the way by which the wavemaker-related functions update against time.
- Output results into files during the time-stepping loop
    - TC1: `energy.csv`, a series of `.npy` binary files with field data, `readme.txt`
    - TC4: `energy.csv`, `probes.csv` with numerical measurements, `readme.txt`
- The above changes are also implemented into `3D_tank.py`. In addition, the Lagrange polynomial $\tilde{\varphi}_i(z)$ is now constructed based on GLL points.
- In `3D_tank_VP.py`, $\hat{\phi}(z)$ can be switched between 1 and high order Lagrange polynomial based on GLL points via the flag `hatphi_one`; while the flag `one_ver_ele` decides whether there is only one element or multiple elements in the $z$-direction.

## Log
| Test Case | New Approach | Old Approach (SV-GLL) |
| :---:     |    :----:    |   :----:     |
| TC1       |**`3D_tank_VP.py`** <br/>`settings_TC1.py`, `savings.py` | **`3D_tank.py`** + `solvers_full.py` <br/>`settings_TC1.py`, `savings.py`  |
| TC3       |**`3D_tank_VP.py`** <br/>`settings_TC3.py`, `savings.py`<br/> :clock1: Δt=0.001s TBD (JC/YL); <br/>:white_check_mark: Δt=0.002s Done. 15h(16p-YL). 20230109 | **`3D_tank.py`** + `solvers_full.py` <br/>`settings_TC3.py`, `savings.py` <br/> :clock1: Δt=0.001s: job submitted to HPC <br/> :white_check_mark: Δt=0.002s: Done. 13h(16p-YL). 20230110 |
| TC4       |**`3D_tank_VP.py`** <br/>`settings_TC4.py`, `savings.py`<br/> folder `202002` <br/> :white_check_mark: Done. 28h(16p-YL). 20230108 |  **`3D_tank.py`** + `solvers_full.py`<br/>`settings_TC4.py`, `savings.py` <br/> folder `202002` <br/> :clock1: TBD (OB/YL)  |
