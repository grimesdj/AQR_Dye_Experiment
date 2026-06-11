clear all
close all

%% Load

% ADCP
fprintf('loading data...\n')
M1 = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_L0_10min.mat');
M1.Config = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_config.mat');

% Mooring
Moor = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L1/mooring_M1.mat');

% apply qcFlag
fprintf('applying qcFlag...\n')
M1.qcFlag(isnan(M1.qcFlag)) = 0;
qcFlag = round(M1.qcFlag); % omg the qcFlag got smoothed
fields = fieldnames(M1);
for i = 1:length(fields)
    field = fields{i};
    dum = M1.(field);
    if size(dum) == size(qcFlag)
        dum(~qcFlag) = NaN;
        M1.(field) = dum;
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
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt);
Moor.Temperature = y';

% common time
fprintf('interpolating to common grid...\n')
Time = M1.Time;
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
        M1.(field) = tmp;
    end
end

%% Plot Velocity
fprintf('plot Velocity...\n')
figure
img = imagesc(Time, M1.bin_mab, M1.Velocity_North);
mask = ~isnan(M1.Velocity_North);
set(img, "AlphaData", mask)
set(gca, 'YDir', 'normal')
colormap(cmocean('balance'))
cmax = max(std(M1.Velocity_North, 'omitmissing'));
clim([-cmax cmax])
cb = colorbar;
ylabel('Meters above bottom')
ylabel(cb, 'North Velocity, [m/s]')
set(gca, 'FontSize', 18)

% have to manually fix Temp for now
top_T = M1.Temperature(1, :);
M1.Temperature = M1.Temperature([1 2 3 4 5 7 8 9 10], :);
top_mab = M1.Temperature_mab(1);
M1.Temperature_mab = M1.Temperature_mab([1 2 3 4 5 6 7 8 10]);

% top sensor is floating
P = M1.Pressure;
b = top_mab;
d = P + b;

% grid Temp
[t_grid, z_grid] = meshgrid(Time, M1.Temperature_mab);
z_grid(1, :) = d;



% add Temp contours
fprintf('add Temp contours...\n')
hold on
minT = ceil(min(M1.Temperature, [], 'all'));
maxT = floor(max(M1.Temperature, [], 'all'));
Tvec = minT:maxT;

cm = cmocean('thermal');
clen = size(cm, 1);
step = floor(clen/length(Tvec));
cstep = 1:step:clen;
cm = cm(cstep, :);
for ib = 1:length(Tvec)
    contour(t_grid, z_grid, M1.Temperature, [Tvec(ib) Tvec(ib)], 'EdgeColor',cm(ib, :), 'LineWidth', 2.5, 'DisplayName', [sprintf('%d ', Tvec(ib)) '$^\circ$C'])
end

% add free surf
plot(Time, M1.Pressure, 'b', 'LineWidth', 2.5, 'DisplayName', 'Free Surface')
lgd = legend;
lgd.Interpreter = 'latex';
lgd.NumColumns = 4;

% uGhhh date ticks
ax = gca;
ax.XTick = min(Time):1/12:max(Time);
datetick('x','HH:MM mm/dd','keeplimits', 'keepticks')
xtickangle(315)


%% Figure
%print(gcf, '../Presentations/tea_06_12_2026/figures/Velocity_contoured', '-dpng', '-r600')
%xlim

%% z/h grid
% pvar = P - mean(P, "omitmissing");
% zh = z_grid + pvar;
% zh(1,:) = zh(1,:) - 2*pvar;

z = z_grid./M1.Pressure;

figure, plot(Time, z)
% interp
dz = 0:0.1:1;

