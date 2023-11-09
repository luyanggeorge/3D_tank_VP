## About
Revisit the [3D wave tank test cases (TCs)](https://github.com/EAGRE-water-wave-impact-modelling/3D-wave-tank-JCP2022) through the "VP approach" where the weak formulations are now generated automatically based on the time-discretised VP via Firedrake's `derivative()` functionality.

## Codes
- Main file: `3D_tank_VP.py`
- Setting of TCn: `settings_TCn.py`

## Developing Log
- :white_check_mark: Change the way by which the wavemaker-related functions update against time (2023.11.9)
- :x: Solver fails to converge starting from the very first time iterative:
  ```
  firedrake.exceptions.ConvergenceError: Nonlinear solve failed to converge after 0 nonlinear iterations.
  Reason:
     DIVERGED_FNORM_NAN
  ``` 
- TO DO
    - Output results into files
    - Change the form of $\hat{\phi}(z)$ by using GLL
