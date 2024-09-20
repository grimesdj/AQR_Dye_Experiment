% 1) first stab at processing ctdF. I'm still working on the GPS points for Release 1
preprocess_Release1_CTDF_casts
preprocess_Release2_CTDF_casts
preprocess_Release3_CTDF_casts
% 
% 1a) make some preliminary plots of the CTD data.
%     modify releaseNumber & ctdSerialNumber at the top
%
plot_L0_CTDF_casts
%
% 2) process the signature 1000 data into 24hr ensembles
%    also estimate hourly wave stats and 5-min average currents
%    modify adcpID to switch between instruments.
%
main_process_SIG1000
%
% 3) estimate fluorometer calibration fits based on the buckets.
%
fluorometer_dye_calibration
%
% 4) load each moored instrument SN and plot L0 data
%    combine L0 data of each instrument on moorings (M1,...,N,S,...,NW,SW,...)
%    and save as L1, make a prelim plot of all vars for each mooring.
%
process_moored_data
%
% 5) load and pre-process release Aquadopp & Vector data
%
load_AquadoppHR_Release1
load_AquadoppHR_Release2
load_Vector_Release1
load_Vector_Release2
