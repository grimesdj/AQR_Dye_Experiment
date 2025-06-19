First Readme: Jason Rooker Summer 2025 
Scratch Directory: code development and testing.

Scripts to load and make quick plots:


aquaplot.m: make quick line plots/histograms and output stats.
USAGE: Stats = aquaplot(releaseNum, version, plots, binStart, binEnd)
       releaseNum = [1 2 or 3] dye release designation
       version    = ['raw', 'L0'] options for Aquadopp without/with qcFlag applied,
       		    'M1' to plot the SIG1000 at M1,
		    'Vector' to plot the Vector velocity at release
       plots      = 'on'/'off' plots timeseries and histograms for user specified release and version. (contains my own QC flag to generate simplified 'L0'); defaults to on, bin start and bin end default to the entire range

       Stats      = [avg_north_velocity, ...]
       


fetch_M1:
	called by aquaplot to fetch M1 data from a specific release based on time

plotRelease_func:
	called by aquaplot. does the actual work of generating and saving plots

rotation:
	called by aquaplot. converts vector XYZ to ENU using AquaDopp heading pitch and roll

statsTable:
	generates a table that displays stats from aquaplot

