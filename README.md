## Author
Dhyaan Radhakrishnamurthy  
LinkedIn: https://www.linkedin.com/in/dhyaan-r-83b0a9198/

# Li-ion 3RC Parameter Estimation (SOC-Dependent ECM)

High-fidelity SOC-dependent 3RC equivalent circuit model (ECM) parameter estimation for Li-ion cells using pulse test data (e.g., HPPC).
Generates Simulink-ready lookup tables of **R0, R1-C1, R2-C2, R3-C3 vs SOC** and validation plots.

## Highlights
- SOC-dependent 3RC ECM parameter estimation at fine SOC resolution (~1–1.5%)
- Automated preprocessing and pulse segmentation workflow
- Optimization-based fitting in MATLAB (and/or Simulink Design Optimization)
- Validation on held-out pulses / dynamic profiles
- Achieved **≤0.8% max voltage error** on validation data (update with your latest figure)

## Repository Structure
- `data/` raw/processed/sample datasets (raw data excluded by default)
- `src/` preprocessing, estimation, validation modules
- `models/` Simulink + MATLAB ECM models
- `results/` estimated parameters and plots
- `scripts/` runnable entry points

## Quick Start
1. Put your dataset exports in `data/raw/` (kept local, not committed)
2. Run the pipeline:
```matlab
scripts/run_pipeline.m


