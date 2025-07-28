% Let's make some 10 min averages:
clear all
close all

% Define the files to be loaded
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release1/L0/*.mat');

% Load the data
data = load('../../../../Kelp_data/Summer2025/Rooker/Release1/L0/KELP1_AquadoppHR_L0.mat');

% Define the filter
Nf = 600/data.Config.dt;
flt = hamming(Nf);
flt = flt/sum(flt);

% Flag NaN values
nanFlag = ~isnan(data.Velocity_East(:,11));
% Convolve the filter and NaN flag
flag_flt = conv(nanFlag, flt, 'same');
% Set NaNs equal to 0
data.Velocity_East(find(~nanFlag), 11)  = 0;
data.Velocity_North(find(~nanFlag), 11) = 0;
% Convolve the filter and data
LPF_East  = conv(data.Velocity_East(:, 11), flt, 'same');
LPF_North = conv(data.Velocity_North(:, 11), flt, 'same');
% Divide the filtered data by the filtered flag
LPF_East  = LPF_East./flag_flt;
LPF_North = LPF_North./flag_flt;


figure;
ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
plot(ax1, data.Time, LPF_East)
plot(ax2, data.Time, LPF_North)
datetick(ax1, 'x', 'keeplimits')
datetick(ax2, 'x', 'keeplimits')
sgtitle('10 min avg velocities')
ylabel(ax1, 'East Velocity')
ylabel(ax2,'North Velocity')
xlabel(ax2, 'Time')