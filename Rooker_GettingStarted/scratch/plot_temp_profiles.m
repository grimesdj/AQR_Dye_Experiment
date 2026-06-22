clear all
close all

%% Load

prof_fig = figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);
strat_fig = figure;

for mooring_ID = 1
moorings = {'M1', 'M2', 'M3', 'M4'}; 
mooring = moorings{mooring_ID};

% ADCP
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
fname = "mooring_" + mooring + ".mat";
M = load(fullfile(fpath, fname));

%% Smooth Temp Data
fprintf('smoothing Temp data...\n')
x = M.Temperature';
wpass = 1/600;
dt = round(median(diff(M.Time)*86400));
fs = 1/dt;
fig = 0;
ds = 1;
nan_filt = 90;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
M.Temperature = y';

% grab float
top_T = M.Temperature(1, :);
top_mab = M.Temperature_mab(1);

% get bottom pressure
if isempty(M.Pressure)
    tmp = load(fullfile(fpath, '..', 'L0', 'ADCP', "ADCP_" + mooring + "_L0_10min.mat"), 'Pressure' , 'Time');
    M.Pressure = interp1(tmp.Time, tmp.Pressure, M.Time, 'linear');
end

% top sensor is floating
P = M.Pressure;
b = top_mab;
d = P + b;

% grid Temp
[t_grid, z_grid] = meshgrid(M.Time(idx), M.Temperature_mab);
z_grid(1, :) = d(idx);

%% z/h grid
sig = -1 * z_grid + nanmean(M.Pressure(idx));
dz = 1:1:nanmean(sig(1, :));
%dz = min(sig(end-1, :)):1:nanmean(sig(1, :));
sig = fillmissing(sig, "spline", 2);

% interpolate to regular grid
F = scatteredInterpolant(t_grid(:), sig(:), M.Temperature(:), 'linear', 'nearest');
[Tq, Zq] = meshgrid(M.Time, dz);
Vq = F(Tq, Zq);

% make avg profile
mean_profile = nanmean(Vq, 2);
std_profile = nanstd(Vq,[], 2);

figure(prof_fig)

plot(ax1, mean_profile, dz, '-s', 'LineWidth', 2, 'DisplayName', mooring)
hold(ax1, "on")


plot(ax2, std_profile, dz, '-s', 'LineWidth', 2, 'DisplayName', mooring)
hold(ax2, "on")

figure(strat_fig)
plot(diff(mean_profile), dz(1:end-1), '-s', 'LineWidth', 2)
hold on
end

title(ax1, '$\mu$ Profile', 'Interpreter','latex', 'FontSize', 20)
xlabel(ax1, '$^\circ$C', 'Interpreter','latex')
ylabel(ax1, 'Relative Depth z/h')
set(ax1, 'FontSize', 18)
grid(ax1, 'minor')
lgd1 = legend(ax1);
lgd1.Location = 'northwest';
ylim([0 1])

title(ax2, '$\sigma$ Profile', 'Interpreter','latex', 'FontSize', 20)
xlabel(ax2, '$^\circ$C', 'Interpreter','latex')
ylabel(ax2, 'Relative Depth z/h')
set(ax2, 'FontSize', 18)
grid(ax2, "minor")
lgd2 = legend(ax2);
lgd2.Location = 'northeast';
linkaxes([ax1 ax2], 'y')

figure(strat_fig)
axis equal
grid minor
ylabel('Relative Depth z/h')
xlabel('$^\circ$C/m', 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
