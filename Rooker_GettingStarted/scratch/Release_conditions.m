%% release conditions

clear all
close all

%% load
releaseNum = 1;

load('../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat')
m10 = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_L0_10min.mat');
inputDir = '../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/';
m = fetch_sig1k(inputDir, datetime(release_times{2, 4}), datetime(release_times{2, 5}));

x = m.Velocity_North;
wpass = 1/800;
fs = 1/0.25;
fig = 0;
ds = 0;
nan_filt = 80;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
N = y';

x = R2.Temperature;
wpass = 1/800;
fs = 1;
fig = 0;
ds = 0;
nan_filt = 80;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
T = y';

figure
imagesc(datenum(m.Time), m10.bin_mab, N)
set(gca, 'YDir', 'normal')
colorbar
colormap(cmocean('balance'))
clim([-0.05 0.05]);
hold on
temps = ceil(min(T, [], "all")):floor(max(T, [], 'all'));
cmap = cmocean('haline');
step = floor(size(cmap, 1)./length(temps));
c = cmap(1:step:end, :);
for l = 1:length(temps)
    contour(R2.Time, R2.mab, T, [temps(l) temps(l)], 'Color', c(l, :), 'linewidth', 2.5)
end
xlim([R2.Time(1) R2.Time(end)+5/24])
datetick('x', 'keeplimits')
