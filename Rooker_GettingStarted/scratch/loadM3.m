%% Load M3

clear all
close all

%% get vars

fin = '../../../../Kelp_data/data/FullExperiment/raw/RDI/ADCP_M3_Full.nc';
info = ncinfo(fin);
vars = info.Variables;

% load every var into a struct
for i = 1:length(vars)
    var = vars(i);
    name = var.Name;
    vname = matlab.lang.makeValidName(name);
    fprintf('loading %s...\n', name)
    data.(vname) = ncread(fin, name);
    if length(size(data.(vname))) == 3
        dum = data.(vname);
        E = sprintf('%s_East', vname);
        N = sprintf('%s_North', vname);
        U = sprintf('%s_Up', vname);
        data.(E) = squeeze(dum(:, :, 1));
        data.(N) = squeeze(dum(:, :, 2));
        data.(U) = squeeze(dum(:, :, 3));
    end
end


% make readable time
data.time = datetime(data.time, 'ConvertFrom', 'posixtime');

% % extract vel
% N = data.vel(:, :, 2);
% M3.Velocity_North = squeeze(N);



% figure
% img = imagesc(data.time, data.range, M3.Velocity_North');
% set(img, 'AlphaData', ~isnan(M3.Velocity_North'))
% set(gca, 'YDir', 'normal')
% set(gca, 'FontSize', 18)
% ylabel('mab')
% colorbar
% clim([-0.5 0.5])
% colormap(cmocean('balance'))
% 
% hold on
% plot(data.time, data.depth, 'k', 'LineWidth', 1, 'DisplayName', 'Free Surface')

% clear surface
mask = data.range' < data.depth-4*mean(diff(data.range));  

% build qcFlag
minCorr = min(data.corr, [], 3)';
minAmp  = min(data.amp, [], 3)';
qcFlag = ones(size(data.vel_North));
qcFlag(~mask) = 0;
qcFlag(minCorr < 60) = 0;
qcFlag(minAmp < 30) = 0;

fields = fieldnames(data);

for i = 1:length(fields)
    field = fields{i};
    x = data.(field);
    if length(size(x)) < 3 && length(x) == length(data.time) & ~strcmp(field, 'time')
        if all(size(x) == size(data.vel_North))
            x(~qcFlag) = NaN;
        end

        wpass = 1/600;
        fs = 1;
        fig = 0;
        ds = 1;
        nan_filt = 80;
        hl = 1;
        [y, ~, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
    else
        y = x;
    end
    data.(field) = y;
end
data.time = data.time(idx);

M3.Velocity_North = data.vel_North';
M3.Velocity_East = data.vel_East';
M3.Velocity_Up = data.vel_Up';
M3.qcFlag = qcFlag';
M3.Pressure = data.depth';
M3.bin_mab = data.range;
M3.Temperature = data.temp';
M3.Time = datenum(data.time);
M3.Time = M3.Time';

% trim to deployment
deployTime  =  datenum('03-Jul-2024 00:00:00');
recoverTime  = datenum('25-Jul-2024 00:00:00');
dep = find(datenum(data.time) > deployTime & datenum(data.time) < recoverTime);

fields = fieldnames(M3);
for i = 1:length(fields)
    field = fields{i};
    dum = M3.(field);
    if size(dum, 2) == length(data.time)
        M3.(field) = dum(:, dep);
    end
end

figure
img = imagesc(M3.Time, M3.bin_mab, M3.Velocity_North);
set(img, 'AlphaData', ~isnan(M3.Velocity_North))
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 18)
ylabel('mab')
colorbar
clim([-0.1 0.1])
colormap(cmocean('balance'))
hold on
plot(M3.Time, M3.Pressure, 'k', 'LineWidth', 1.5)

save('../../../../Kelp_data/data/2024_PROCESSED_DATA/M3/L0/ADCP/ADCP_M3_L0_10min.mat', '-struct', "M3")