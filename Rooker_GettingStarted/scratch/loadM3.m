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

    
% build qcFlag
minCorr = min(data.corr, [], 3)';
minAmp  = min(data.amp, [], 3)';
qcFlag = ones(size(M3.Velocity_North));
qcFlag(isnan(M3.Velocity_North)) = 0;
qcFlag(minCorr < 60) = 0;
qcFlag(minAmp < 30) = 0;
VelN = M3.Velocity_North;
VelN(~qcFlag) = NaN;

%% Smooth data
x = VelN;
wpass = 1/600;
fs = 1;
fig = 0;
ds = 0;
nan_filt = 80;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);


M3.Velocity_North = y';
M3.qcFlag = qcFlag';
% trim pre-dep
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
        M3.(field) = M3.(field)(:, 1:600:end);
    end
end

figure
img = imagesc(data.time, data.range, y');
set(img, 'AlphaData', ~isnan(y'))
set(gca, 'YDir', 'normal')
set(gca, 'FontSize', 18)
ylabel('mab')
colorbar
clim([-0.1 0.1])
colormap(cmocean('balance'))

save('../../../../Kelp_data/data/2024_PROCESSED_DATA/M3/L0/ADCP/ADCP_M3_L0_10min.mat', '-struct', "M3")