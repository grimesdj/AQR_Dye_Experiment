First Readme:


Scratch Directory: use this to play around 

aquaplot.m: 
	
	takes arguments ('releaseNum', 'version', plots on/off, 'binStart', 'binEnd')
	
	plots timeseries and histograms for user specified release and version. (contains my own QC flag to generate simplified 'L0') 
	plots defaults to on, bin start and bin end default to the entire range

fetch_M1:
	called by aquaplot to fetch M1 data from a specific release based on time

plotRelease_func:
	called by aquaplot. does the actual work of generating and saving plots

rotation:
	called by aquaplot. converts vector XYZ to ENU using AquaDopp heading pitch and roll

statsTable:
	generates a table that displays stats from aquaplot

