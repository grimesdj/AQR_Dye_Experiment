clear all
close all

%% Load

prof_fig = figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);
strat_fig = figure;

for mooring_ID = 1
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
X = Time(1):1/3600:Time(end);
Y = M.bin_mab(1):0.01:M.bin_mab(end)';
moor_fields = fieldnames(Moor);
for i = 1:length(moor_fields)
    field = moor_fields{i};
    dum = Moor.(field);
    if ~strcmp(field, 'Time') && ~isempty(dum)
        if size(dum, 2) == size(Moor.Time, 2)
            if isvector(dum)
                tmp = interp1(Moor.Time, dum, X);
            else
                tmp = interp1(Moor.Time, dum', X)';
            end
        else
            tmp = dum;
        end
        M.(field) = tmp;
    end
end


%% Plot Velocity
[Xq, Yq] = meshgrid(X, Y);
Vq = interp2(Time, M.bin_mab, M.Velocity_North, Xq, Yq);

fprintf('plot Velocity...\n')
vel_fig(mooring_ID) = figure;
a1 = axes('Position',[0.1229 0.1801 0.7327 0.7449]);
img = imagesc(X, Y, Vq);
mask = ~isnan(Vq);
set(img, "AlphaData", mask)
set(gca, 'YDir', 'normal')
colormap(cmocean('balance'))
cmax = max(std(M.Velocity_North, 'omitmissing'));
clim([-0.1 0.1])
cb = colorbar;
ylabel('Meters above bottom')
ylabel(cb, 'North Velocity, [m/s]')
set(gca, 'FontSize', 18)
Rstart = datenum(datetime('03-Jul-2024 18:00:00'));
Rend = datenum(datetime('03-Jul-2024 23:00:00'));
xlim([Rstart Rend]);
ylim([0 12.5])
drawnow

% if mooring_ID == 1
%     % have to manually fix Temp for now
%     print('Removing Duplicate Sensors...')
% 
%     M.Temperature_mab = M.Temperature_mab([1 2 3 4 5 6 7 8 10]);
% elseif mooring_ID ==2
%     print('Switcing mixed signals...')
%     dum = M.Temperature(5, :);
%     M.Temperature(5, :) = M.Temperature(6, :);
%     M.Temperature(6, :) = dum;
% end


% add ADCP temp to bottom
dum = M.Temperature;
tmp = interp1(Time, ADCP_temp, X);
M.Temperature = [dum ; tmp];
dum = M.Temperature_mab;
M.Temperature_mab = [dum ; 0];

% grab float
top_T = M.Temperature(1, :);
top_mab = M.Temperature_mab(1);

% top sensor is floating
P = M.Pressure;
b = top_mab;
d = P + b;
d = interp1(Time, d, X);

% grid Temp
[t_grid, z_grid] = meshgrid(X, M.Temperature_mab);
z_grid(1, :) = d;

%Tq = interp2(t_grid, z_grid, M.Temperature, Xq, Yq);

% add Temp contours
fprintf('add Temp contours...\n')
hold on
% minT = ceil(min(M.Temperature, [], 'all'));
% maxT = floor(max(M.Temperature, [], 'all'));
% Tvec = minT:maxT;
Tstd = round(std(M.Temperature, [], "all"));
Tbar = round(mean(M.Temperature, 'all'));
minT = Tbar - Tstd;
maxT = Tbar + Tstd;
Tvec = minT:0.5:maxT;


cm = cmocean('themral');
clen = size(cm, 1);
step = floor(clen/length(Tvec));
cstep = 1:step:clen;
cm = cm(cstep, :);
for ib = 1:length(Tvec)
    fprintf('adding %d C contour...\n', Tvec(ib))
    contour(t_grid, z_grid, M.Temperature, [Tvec(ib) Tvec(ib)], 'EdgeColor',cm(ib, :), 'LineWidth', 2.5, 'DisplayName', [sprintf('%.1f ', Tvec(ib)) '$^\circ$C'])
    drawnow
end

% c1 = colorbar;
% c1.Location = 'north';
% c1.Color = cmocean('Thermal');

% add free surf
fprintf('adding free surface...\n')
plot(Time, M.Pressure, 'k', 'LineWidth', 2.5, 'DisplayName', 'Free Surface')
lgd = legend;
lgd.Interpreter = 'latex';
lgd.NumColumns = 4;

% uGhhh date ticks
ax = gca;
ax.XTick = ceil(min(Time)):1/24:floor(max(Time));
datetick('x','HH:MM','keeplimits', 'keepticks')
xtickangle(315)


%% Figure
%print(gcf, '../Presentations/tea_06_12_2026/figures/Velocity_contoured', '-dpng', '-r600')
%xlim

return
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



%% PSD









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
axis square













return
%% Export Figs
print(prof_fig, '../../../../Kelp_data/Summer2025/Rooker/figures/mooring_avg_and_std_profiles.png', '-dpng', '-r600')
print(vel_fig(1), '../../../../Kelp_data/Summer2025/Rooker/figures/M1_velocity_and_temp.png', '-dpng', '-r600')
print(vel_fig(2), '../../../../Kelp_data/Summer2025/Rooker/figures/M2_velocity_and_temp.png', '-dpng', '-r600')
print(vel_fig(3), '../../../../Kelp_data/Summer2025/Rooker/figures/M3_velocity_and_temp.png', '-dpng', '-r600')


