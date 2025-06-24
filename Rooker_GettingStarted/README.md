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
	- R3: Aquadopp battery starts to die and takes limited samples at regular intervals. AQD and 	Vector are similar but the Sig1000 is at a slightly different angle and has a greater magnitude.
  
- [ ] draft readme's for each script worth its weight in salt. 




## Velocity Processing
- Standardizing velocity variables and data files across instruments.
  - [ ] Time, Velocity_East, Velocity_North, Amplitude_Minimum, Correltation_Minimum, etc.
  - [ ] Convert script load/process algorithms to functions, e.g.,
     	```load_AquadoppHR_Release2.m```
	calls:
     	```load_AQD_data.m```
     	then does several other operations (qa/qc) and makes plots. Put that into an ```L0``` function?

---------------------------------------------------------------------
What scripts need to be changed?
	[x] load_AQD_data.m
	[ ] make load_Aquadopp_release1 into function call
	[ ] release 2
	[ ] make load_Vector_release1 into function call
	[ ] release 2
	[ ] load_VECTOR_data.m
	[ ] load_and_process_sig1000_to_RDI_matrix_format.m
	[ ] time_average_and_rotate_sig1000_RDI_matrix_format.m
	[ ] estimate_wave_bulk_stats_SIG1000_RDI_matrix_format.m
	[x] rotation (add to multiplot) -> now called ```Vector_rotation.m```

		- find out what each function is supposed to return and then set that to output in the same syntax as the wrapper function
		- Keep load_Aquadopp_release1 and related functions as wrapper and make a new function that those the actual processing
			(e.g.) lAr1 is user input data and make 2 fuction calls.
		- compare processing; can L0 for Vector be the same function as L0 for Vector?
		- Leave .raws alone, but standardize structures for L0 output in above format.
		- double-check that release 2 L0 and raw makes it to the correct folder

		- take a look at extended velocity range mode
		rotation matrix fro vector
---------------------------------------------------------------------

- [ ] (re) Make L0 files reflecting above changes.
- generate L1 files for subsequent analysis.
  - [ ] Make 5-min averaged fields for AQD and M1. See ```time_average_and_rotate_sig1000_RDI_matrix_format.m``` for example.
  - [ ] For Vector, want to use heading from AQD to make ENU velocities. Then generate 15-30 minute spectra. 

