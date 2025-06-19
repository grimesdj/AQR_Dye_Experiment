# Rooker Summer 2025 Project README

## Preliminary analysis of dye release lander velocity data
- compare raw/L0 aquadopp velocity with Vector (at approriate bin)
  -- what script does this?
  -- where to find the plots/stats
  -- preliminary observations (R1,R2,R3)
- compare raw/L0 aquadopp with M1-signature 1000 (averaged over appropriate bin range)
  -- ...
- compare all three at bin corresponding to Vector
  -- ...
- [ ] draft readme's for each script worth its weight in salt. 

## Velocity Processing
- Standardizing velocity variables and data files accross instruments.
  -- [ ] Time, Velocity_East, Velocity_North, Amplitude_Minimum, Correltation_Minimum, etc.
  -- [ ] Convert script load/process algorithms to functions, e.g.,
     	```load_AquadoppHR_Release2.m'''
	calls:
     	```load_AQD_data.m'''
     	then does several other operations (qa/qc) and makes plots. Put that into an ```L0''' function?
- [ ] (re) Make L0 files reflecting above changes.
- generate L1 files for subsequent analysis.
  -- [ ] Make 5-min averaged fields for AQD and M1. See ```time_average_and_rotate_sig1000_RDI_matrix_format.m''' for example.
  -- [ ] For Vector, want to use heading from AQD to make ENU velocities. Then generate 15-30 minute spectra. 

