clear all
close all

%% Load

prof_fig = figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);
strat_fig = figure;

for mooring_ID = 1:3
moorings = {'M1', 'M2', 'M3'}; % M3 doesnt have ADCP data yet
mooring = moorings{mooring_ID};

% ADCP
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L0', 'ADCP');
fname = "ADCP_" + mooring + "_L0_10min.mat";
M = load(fullfile(fpath,fname));
%M.Config = load(fullfile(fpath, filesep, "ADCP_" + mooring + "_config.mat"));


% Mooring
mpath = fullfile(fpath, '..', '..', 'L1', "mooring_" + mooring + ".mat");
Moor = load(mpath);

% apply qcFlag
fprintf('applying qcFlag...\n')
M.qcFlag(isnan(M.qcFlag)) = 0;
qcFlag = round(M.qcFlag); % omg the qcFlag got smoothed
fields = fieldnames(M);
for i = 1:length(fields)
    field = fields{i};
    dum = M.(field);
    if size(dum) == size(qcFlag) & ~strcmp(field, 'qcFlag')
        dum(~qcFlag) = NaN;
        M.(field) = dum;
    end
end

%% Smooth Temp Data
fprintf('smoothing Temp data...\n')
x = Moor.Temperature';
wpass = 1/600;
dt = round(median(diff(Moor.Time)*86400));
fs = 1/dt;
fig = 0;
ds = 0;
nan_filt = 90;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
Moor.Temperature = y';

ADCP_temp = fillmissing(M.Temperature, 'linear');

% common time
fprintf('interpolating to common grid...\n')
Time = M.Time;
moor_fields = fieldnames(Moor);
for i = 1:length(moor_fields)
    field = moor_fields{i};
    dum = Moor.(field);
    if ~strcmp(field, 'Time') && ~isempty(dum)
        if size(dum, 2) == size(Moor.Time, 2)
            if isvector(dum)
                tmp = interp1(Moor.Time, dum, Time);
            else
                tmp = interp1(Moor.Time, dum', Time)';
            end
        else
            tmp = dum;
        end
        M.(field) = tmp;
    end
end

% add ADCP temp to bottom
dum = M.Temperature;
M.Temperature = [dum ; ADCP_temp];
dum = M.Temperature_mab;
M.Temperature_mab = [dum ; 0];

% grab float
top_T = M.Temperature(1, :);
top_mab = M.Temperature_mab(1);

% top sensor is floating
P = M.Pressure;
b = top_mab;
d = P + b;

% grid Temp
[t_grid, z_grid] = meshgrid(Time, M.Temperature_mab);
z_grid(1, :) = d;

%% z/h grid
sig = z_grid./M.Pressure;
dz = min(sig(end-1, :)):0.1:min(sig(1, :));
sig = fillmissing(sig, "spline", 2);

% interpolate to regular grid
F = scatteredInterpolant(t_grid(:), sig(:), M.Temperature(:));
[Tq, Zq] = meshgrid(Time, dz);
Vq = F(Tq, Zq);

% make avg profile
mean_profile = mean(Vq, 2);
std_profile = std(Vq,[], 2);

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
