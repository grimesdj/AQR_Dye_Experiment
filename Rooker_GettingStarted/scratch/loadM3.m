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
end


% make readable time
data.time = datetime(data.time, 'ConvertFrom', 'posixtime');

% extract time
N = data.vel(:, :, 2);
M3.Velocity_North = squeeze(N);

% clear abover surface
mask = data.range' < data.depth-4*mean(diff(data.range));   
M3.Velocity_North(~mask) = NaN;

figure
img = imagesc(data.time, data.range, M3.Velocity_North');
set(img, 'AlphaData', ~isnan(M3.Velocity_North'))
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 18)
ylabel('mab')
colorbar
clim([-0.5 0.5])
colormap(cmocean('balance'))

hold on
plot(data.time, data.depth, 'k', 'LineWidth', 1, 'DisplayName', 'Free Surface')

%% Smooth data
x = M3.Velocity_North;
wpass = 1/600;
fs = 1;
fig = 0;
ds = 0;
nan_filt = 80;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);


figure
img = imagesc(data.time, data.range, y');
set(img, 'AlphaData', ~isnan(y'))
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 18)
ylabel('mab')
colorbar
clim([-0.1 0.1])
colormap(cmocean('balance'))


M3.Velocity_North = y';
M3.qcFlag = ~isnan(M3.Velocity_North);
% finish constucting qcFlag
% trim pre-dep
M3.Pressure = data.depth';
M3.bin_mab = data.range;
M3.Temperature = data.temp';
M3.Time = datenum(data.time);
M3.Time = M3.Time';

save('../../../../Kelp_data/data/2024_PROCESSED_DATA/M3/L0/ADCP/ADCP_M3_L0_10min.mat', '-struct', "M3")