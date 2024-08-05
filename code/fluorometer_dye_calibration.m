clear all
close all
%
% 1) load the calibration info (serial numbers, times, thermistors)
infoRoot = '~/OneDriveUNCW/KELP-vingilote/info/dye_mixing_cals_and_releases/'
fluoRoot = '~/OneDriveUNCW/KELP-vingilote/info/dye_mixing_cals_and_releases/'
thrmRoot =
%
calibration_times_file = [infoRoot,'dye_bucket_times_vs_serial_number.csv'];
thermistor_info_file   = [infoRoot,'dye_bucket_thermistor_serial_numbers.csv'];
%
% 2) for each serial number: determine the instrument type, load the calibration data,...
% 3) segment into 30 second or 5 minute bucket time-series.
% 4) least squares fit regression (slope/intercept) and save in a fluorometer calibration file.