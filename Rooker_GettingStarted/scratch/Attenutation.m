%% Attenuation Coef

clear all
close all
fprintf("==================================================\nWelcome to Jason's Attenuation Calculator!\n==================================================\n")

%% Load data
filestem = '../../../../Kelp_data/Summer2025/Rooker/Release2/LPF';

% Get all .mat files that aren't PCA
files_all = dir([filestem, '/*.mat']);
names = string({files_all.name});
keep_idx = ~contains(names, 'PCA');
files = files_all(keep_idx);

% distrubute
AQD = load(fullfile(files(1).folder, files(1).name));
M1  = load(fullfile(files(2).folder, files(2).name));
M2  = load(fullfile(files(3).folder, files(3).name));
VEC = load(fullfile(files(4).folder, files(4).name));

%% Process data
fprintf('Just going to look at M1 and M2 for now...\n')

% allocate vars
time    = M1.Time;
dtime   = datetime(time, 'ConvertFrom', 'datenum');
u_out   = M1.Velocity_X;
v_out   = M1.Velocity_Y;
u_in    = M2.Velocity_X;
v_in    = M2.Velocity_Y;

% plot so it makes some sense
fprintf('\nFigure 1 shows the LPF values of M1 and M2 over each other\n')
figure
ax1 = subplot(2, 1, 1);
plot(ax1, dtime, u_out, 'r', 'LineWidth', 1)
hold on
plot(ax1, dtime, u_in, 'b', 'LineWidth', 1)
legend('M1 outside', 'M2 inside')
title('East Velocities', 'FontSize', 16)

ax2 = subplot(2, 1, 2);
plot(ax2, dtime, v_out, 'r', 'LineWidth', 1)
hold on
plot(ax2, dtime, v_in, 'b', 'LineWidth', 1)
legend('M1 outside', 'M2 inside')
title('North Velocities', 'FontSize', 16)

% I downloaded some tide data, maybe I'll take a look at that later
% (consider removing or isolating for tidal current velocity comparisons?)
% find the .csv at
% '../../../../Kelp_data/Summer2025/Rooker/Release2/tide/Tides_MTL_9411340_GMT.csv'
fprintf('You can clearly see the tidal variations...\n')


% detrend
r_out   = detrend(u_out);
k_out   = detrend(v_out);
r_in    = detrend(u_in);
k_in    = detrend(v_in);

fprintf('\nFigure 2 shows the de-trended verison\n')
figure
ax1 = subplot(2, 1, 1);
plot(ax1, dtime, r_out, 'r', 'LineWidth', 1)
hold on
plot(ax1, dtime, r_in, 'b', 'LineWidth', 1)
legend('M1 outside', 'M2 inside')
title('East Velocities', 'FontSize', 16)

ax2 = subplot(2, 1, 2);
plot(ax2, dtime, k_out, 'r', 'LineWidth', 1)
hold on
plot(ax2, dtime, k_in, 'b', 'LineWidth', 1)
legend('M1 outside', 'M2 inside')
title('North Velocities', 'FontSize', 16)

fprintf('Does not look much different, not surprising...\n\n')

figstem = '../../../../Kelp_data/Summer2025/Rooker/figures/attenuation/';
exportgraphics(gcf, [figstem, 'detrend_EN_M1M2.png'])

%% Caclulate Attenuation

% make sigmas
sig_r_out = std(r_out);
sig_r_in  = std(r_in);
sig_k_out = std(k_out);
sig_k_in  = std(k_in);

A_u = sig_r_in/sig_r_out;
A_v = sig_k_in/sig_k_out;

fprintf('East Attenuation Coef = %.2f\n', A_u)
fprintf('North Attenuation Coef = %.2f\n', A_v)

%% Run Correlation Stats
Corr_u = corrcoef(r_in, r_out);
Corr_v = corrcoef(k_in, k_out);

fprintf('\nCorrelations still look pretty bad...\n')

%% Lets take a quick look at magnitude

mag_in  = sqrt(u_in.^2 + v_in.^2);
mag_out = sqrt(u_out.^2 + v_out.^2);

figure
scatter(mag_in/std(mag_in), mag_out/std(mag_out), 'black', 'filled')
axis equal
title('Magnitudes in vs out')

[cxy, f] = mscohere(mag_in, mag_out, [], [], [], 1/600);
T = 1./f;
figure
semilogx(T, cxy, 'LineWidth', 1.5)
xlim([600 86400]);
xline((360/28.984104)*60*60, 'r--', 'Principal lunar semidiurnal constituent', 'LabelVerticalAlignment','bottom')
xline((360/15.041069)*60*60, 'r--', 'Lunar diurnal constituent', 'LabelVerticalAlignment','bottom')
set(gca, 'XDir', 'reverse')
title('M1 vs M2 Magnitudes Coherence', 'FontSize', 20)
xlabel('Period, T [s]', 'FontSize', 16)
ylabel('Correlation Coefficient, r [\phi]', 'FontSize', 16)

exportgraphics(gcf, [figstem, 'magnitude_coherence.png'])

% Looks like theres some correlation on the order of 2 hour periods and 1
% hour periods. That kind of harmonics seems like stokes 2nd order waves?

get(figure(2))
fig2 = gcf;
fig5 = gcf;
get(fig5)
linkaxes([ax1 ax2], 'x')
xlim([datetime('09-Jul-2024 00:00:00') datetime('09-Jul-2024 10:00:00')])

exportgraphics(gcf, [figstem, 'detrend_zoom.png'])
