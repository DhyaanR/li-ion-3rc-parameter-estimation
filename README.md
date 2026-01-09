## Overview
This repository implements SOC-dependent parameter estimation of a Li-ion
battery using a 3RC Equivalent Circuit Model (ECM).
The model parameters (R0, R1–R3, C1–C3) are identified using nonlinear
least-squares fitting on pulse-discharge segments.

## Methodology
1. Load synthetic pulse-discharge cell data
2. Detect discharge pulse start points using current
3. Segment voltage/current between successive discharge pulses
4. Compute SOC by coulomb integration (SOC start = 100%, end = 0%)
5. Fit a 3RC ECM per segment using nonlinear least squares
6. Assign estimated parameters to corresponding SOC points

## Scripts
- `01_load_synthetic_data.m`
- `02_segment_pulses_and_soc.m`
- `03_LS_fit_all_segments_3rc.m`

## Outputs
- SOC-dependent R0, R1–R3, C1–C3 tables (MATLAB workspace)
- Individual plots of each R and C parameter vs SOC
- Fit quality metrics (RMSE, max error)

## Requirements
- MATLAB
- Optimization Toolbox (for lsqnonlin)

## Author
Dhyaan Radhakrishnamurthy  
M.Eng Electrical & Computer Engineering, University of Ottawa  
LinkedIn: https://www.linkedin.com/in/dhyaan-r-83b0a9198/
