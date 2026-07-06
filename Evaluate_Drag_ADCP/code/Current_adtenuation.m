%% 

clear all
close all


%% Load
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M1', 'L1', 'ADCP');
fname = '.mat';
M1 = load(fullfile(fpath, fname));