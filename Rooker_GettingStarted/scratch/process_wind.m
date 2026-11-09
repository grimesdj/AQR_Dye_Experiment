%% Process Wind data

clear all
close all

%% load
Data = load('../../../../Kelp_data/data/wind/2024_wind_data_nearshore.mat');

%% fliter direction
X = Data.WSPD_m_s .* cosd(Data.WDIR_degT);
Y = Data.WSPD_m_s .* sind(Data.WDIR_degT);

%% low pass
dt = seconds(Data.Time(2)-Data.Time(1));
fs = 1/dt;
wpass = 1/(6*3600);

[y, fsd, idx] = hamming_filter(X, wpass, fs);
North = y;

[y, fsd, idx] = hamming_filter(Y, wpass, fs);
East = y;

Time = Data.Time(idx);

figure
plot(Time, North, 'k', 'LineWidth', 2)
yline(0, 'r--')

