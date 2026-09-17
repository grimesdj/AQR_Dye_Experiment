%% main temp anf velocity

clear all
close all

%% Grid Vars

% velocity
fprintf('interpolating velocity to depth grid...\n')
grid_velocity

% temperature
fprintf('interpolating temperature to depth grid...\n')

%% EOFs

% velocity
fprintf('generating velocity EOFs...\n')
velocity_EOF

% temp
fprintf('generating temp EOFs...\n')
temp_EOFs

%% Compare  Vel and Temp
fprintf('Comparing Temp and Velocity...\n')
temp_vs_velocity_spectra

%% Estimate variance
fprintf('Estimating Variance breakdown...\n')
Variance_Bar_Plot
