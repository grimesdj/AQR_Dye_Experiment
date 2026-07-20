%% release conditions

clear all
close all

%% load
releaseNum = 1;

load('../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat')
load('../../../../Kelp_data/data/Release2/L0/KELP2_Aquadopp_L0.mat')

x = north;
wpass = 1/300;
fs = 1;
fig = 0;
ds = 0;
nan_filt = 80;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
N = y';

x = R2.Temperature;
wpass = 1/300;
fs = 1;
fig = 0;
ds = 0;
nan_filt = 80;
hl = 1;
[y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt, hl);
T = y';

figure
imagesc(datenum(time), dbins, N)
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
    contour(R2.Time, R2