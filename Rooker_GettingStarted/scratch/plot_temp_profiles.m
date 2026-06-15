%% Plot Temperature Profiles

clear all
close all

%% Load Data


M1 = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L1/mooring_M1.mat');
M1_ADCP = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_L0_10min.mat');
M1.Pressure = M1_ADCP.Pressure;


M2 = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M2/L1/mooring_M2.mat');
M3 = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M3/L1/mooring_M3.mat');


% fix M1 Temps
top_T = M1.Temperature(1, :);
M1.Temperature = M1.Temperature([1 2 3 4 5 7 8 9 10], :);
top_mab = M1.Temperature_mab(1);
M1.Temperature_mab = M1.Temperature_mab([1 2 3 4 5 6 7 8 10]);

