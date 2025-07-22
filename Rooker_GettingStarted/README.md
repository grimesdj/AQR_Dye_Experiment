# Rooker Summer 2025 Project README

## Preliminary analysis of dye release lander velocity data


  - multiplot.m plots East and North velocities for all three instruments at bin equivalent to vector height 
  - aquaplot.m plots ENU, beam velocity, beam amplitude, correlation, and statistics for all bins for a single instrument
 - Plots and figures for velocities are in '../Kelp_data/Summer2025/Rooker/figures/...'
  - Statistics table in '..Kelp_data/Summer2025/Rooker/Stats/...'
  - preliminary observations (R1,R2,R3): 

- compare raw/L0 aquadopp velocity with Vector (at appropriate bin)

- compare raw/L0 aquadopp with M1-signature 1000 (averaged over appropriate bin range)
	- R1: Bin + Time averaged velocities are more than 90 degrees apart
	- R2: Trajectories are similar, with AQD magnitude being smaller 
	- R3: Trajectories are similar but further apart than R2, AQD magnitude is much smaller


- compare all three at bin corresponding to Vector
	- R1: Data from AquaDopp appears to be phase wrapped and is not usable in its current state. 
		Vector and Sig1000 have average velocities nearly perpendicular to each other
	- R2: All three instruments have time averaged currents of similar magnitude and direction.
	- R3: Aquadopp battery starts to die and takes limited samples at regular intervals. AQD and Vector are similar but the Sig1000 is at a slightly different angle and has a greater magnitude.
  
- [ ] draft readme's for each script worth its weight in salt. 




## Velocity Processing
- Standardizing velocity variables and data files across instruments.
  - [x] Time, Velocity_East, Velocity_North, Amplitude_Minimum, Correlation_Minimum, etc.
  - [x] Convert script load/process algorithms to functions, e.g.,
     	```load_AquadoppHR_Release2.m```
	calls:
     	```load_AQD_data.m```
     	then does several other operations (qa/qc) and makes plots. Put that into an ```L0``` function?

  - [x] use jolt filter to figure out what stuff got messed up by the .hr2 and the fix it :) -> ended up using Shcherbinas's magic unwrapping script	


L0:
	- [x] (re) Make L0 files reflecting above changes.

	- [x] add principal axis rotation

			[x] 2 panel plots:
				PCA_X, PCA_Y
				Minamp, Mincor

	- [x] main_sig1k put vars into config struct
	- [ ] descriptions
	- [x] add to standardized variable list (sig 1000) to load and process Burst_AltimeterDistanceAST (or something) and 
		Burst_...TimeOffset 
	- [x] plots comparing mag/direction of M1, AQD, Vector for all releases
L1:
- generate L1 files for subsequent analysis.
 	- [ ] Make 5-min averaged fields for AQD and M1. See ```time_average_and_rotate_sig1000_RDI_matrix_format.m``` for example.
		- for all variables that depend on time, apply an appropriate width Hamming/Hanning filter using conv() or conv2()
		- then sub-sample using the appropriate span (equal to half-width of filter)
 	- [x] For Vector, want to use heading from AQD to make ENU velocities. 
	- [ ]Then generate 15-30 minute spectra. 
 

