## Post-processing code for the TCs

| Test Case| Figure No. | Data Sources | Post-Processing code  |
|  :----:  |  :----:    |    :----:    |      :----:           |
|   TC1    |   Fig.3    |  `energy.csv`, `.npy` field data files | `pp_SWS_TC1VP.py` |
|   TC2    |   Fig.4    |  `.npy` field data files | Fig.4(a),(b): `pp-TC1-convergence-dx.py`<br> Fig.4(c),(d): `pp-TC1-convergence-nz.py`<br> Fig.4(e),(f): `pp-TC1-convergence-dt.py`  |
|   TC2    |   Fig.5, Fig.6    |  `.npy` field data files |  `pp-TC1-advanced-convergence.py` |
|   TC3    |   Fig.7, Fig.8 |  `checkpoints.csv` | `pp_energy_figs_TC3VP.py` |
|   TC4    |   Fig.10    |  `probes.csv`, folder `202002` | `pp_wavemaker_TC4VP.py`  |
|   TC4    |   Fig.12   |  `probes.csv`, folder `202002` | `pp_probes_TC4VP.py`  |
|   TC4    |   Fig.13   |  `probes.csv`, folder `202002` | `FFT_202002.m`  |

**Notes:**
- Figure numbers are subject to change until the paper is published.
