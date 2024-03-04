## Post-processing code for the TCs

| Test Case| Figure No. | Data Sources | Post-Processing code  |
|  :----:  |  :----:    |    :----:    |      :----:           |
|   TC1    |   Fig.3    |  `energy.csv`, `.npy` field data files | `pp_SWS_TC1VP.py` |
|   TC2    |   Fig.5    |  `.npy` field data files | `pp_Atiken_3groups_TC2VP.py` |
|   TC3    |   Figs.6&7 |  `checkpoints.csv` | `pp_energy_figs_TC3VP.py` |
|   TC4    |   Fig.9    |  `probes.csv`, folder `202002` | `pp_wavemaker_TC4VP.py`  |
|   TC4    |   Fig.11   |  `probes.csv`, folder `202002` | `pp_probes_TC4VP.py`  |
|   TC4    |   Fig.12   |  `probes.csv`, folder `202002` | `FFT_202002.m`  |

**Notes:**
- Figure numbers are subjected to changes as the paper is not finalised yet.
- For TC2, which part of the field data to be loaded and processed can be specified via the return value of the function `load_data`, and `weight` should also be modified accordingly.
