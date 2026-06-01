%% Plot Profiles

clear all
close all

% Load M1
M1.ADCP = load('../../../../Kelp_data/data/Release2/L0/ADCP/M1_ADCP.mat');
M1.LPF = load('../../../../Kelp_data/data/Release2/L1/ADCP/M1_LPF_600sec.mat');
M1.Moor = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L1/mooring_M1.mat');

releasenum = 2;
release = string(releasenum);


depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00') ;
            datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];


%trim times
M1.LPF.Time   = M1.ADCP.Time(M1.LPF.idx); % pls fix this in the LPF script omg
M1.LPF.valid  = M1.LPF.Time >= depTime(releasenum, 1) & M1.LPF.Time <= depTime(releasenum, 2);
M1.ADCP.valid = M1.ADCP.Time >= depTime(releasenum, 1) & M1.ADCP.Time <= depTime(releasenum, 2);
M1.Moor.valid = M1.Moor.Time >= depTime(releasenum, 1) & M1.Moor.Time <= depTime(releasenum, 2);

% turn so it fits imagesc
North = M1.LPF.Velocity_North';
North = North(:,M1.LPF.valid);
Vtime = M1.LPF.Time(M1.LPF.valid)';

% some errors with the loading so lets fix that for now
Temp_cont = M1.Moor.Temperature([2 3 4 5 7 8 9 10], M1.Moor.valid);
%M1.Moor.Temperature_mab(1) = rms(M1.ADCP.maxRNG, 'all')-M1.Moor.Temperature_mab(1);
%Temp_cont(6, :) = nan;
Temp_mab = M1.Moor.Temperature_mab([2 3 4 5 7 8 9 10]);


% smooth temp??
dt = median(diff(M1.Moor.Time)) * 86400;
[Temp, fsd, idx] = hamming_filter(Temp_cont', 1/600, 1/dt, 1, 1);
Temp = Temp';



% plot
figure
img = imagesc(Vtime, M1.ADCP.bin_mab, North);
set(img, 'AlphaData', ~isnan(North))
set(gca, 'YDir', 'normal')
cb = colorbar;
colormap(cmocean('balance'))
clim([-0.1 0.1])
ylabel('$z$ [m]', 'Interpreter','latex')
ylabel(cb, 'Cross-shore Velocity, $u$ [m/s]', 'Interpreter','latex', 'FontSize', 18)
set([gca cb], 'fontsize', 18)
ylim([0 max(M1.ADCP.maxRNG)+std(M1.ADCP.maxRNG)])
xt = min(Vtime):6/24:max(Vtime);
xticks(xt)
datetick('x', 'mm/dd HH:MM', 'keepticks')
xlabel('Time (UTC)')


% add temp contours
hold on
Ttime = M1.Moor.Time(M1.Moor.valid);
Ttime = Ttime(idx);
contour(Ttime,Temp_mab,  Temp, [16 17 17.5 18], 'k', 'LineWidth', 1.5)
yline(1, 'k--', 'LineWidth', 1.5, 'label', 'Dye Release Depth', 'FontSize', 18, 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom', 'FontWeight','bold')

% hold on
% [C,h] = contour( ...
%     M1.Moor.Time(M1.Moor.valid), ...
%     M1.Moor.Temperature_mab(:), ...
%     M1.Moor.Temperature(:,M1.Moor.valid), ...
%     1, ...
%     'k', ...
%     'LineWidth',1);
%figure, imagesc(M1.Moor.Time(M1.Moor.valid), M1.Moor.Temperature_mab, M1.Moor.Temperature(:,M1.Moor.valid))

% interpolate to common grid
bin_mab = M1.ADCP.bin_mab;
[T_grid, Z_grid]        = meshgrid(Ttime, Temp_mab);
[Vtime_grid, bin_grid]  = meshgrid(Vtime, bin_mab);
terp = interp2(T_grid, Z_grid, Temp, Vtime_grid, bin_grid);
T = terp;
N = North;
% remove nans
% 
% mask = ~isnan(terp) & ~isnan(North);
% T = terp(mask);
% N = North(mask);

badRows = all(isnan(terp),2);

T  = terp(~badRows,:);
N = North(~badRows,:);
bin_mab = bin_mab(~badRows);

T = fillmissing(T, 'linear', 2);
N = fillmissing(N, 'linear', 2);

% remove mean
T = T - mean(T);
N = N - mean(N);

% detrend
T = detrend(T);
N = detrend(N);

% remove diurnal thermal exchange
%dt = median(diff(Vtime)) * 86400;
%DTE = hamming_filter(T, 1/(6*60*60), 1/dt, 1, 0);
%ST = T-DTE;


figure
scatter(terp(:), North(:), 'cyan', 'filled')
xlabel('Tempurature $^{\circ}$C', 'Interpreter','latex')
ylabel('North Velocity [m/s]')
set(gca, 'FontSize', 18)

T_norm = T ./ std(T, [], 'all');
N_norm = N ./ std(N, [], 'all');

figure
scatter(T_norm, N_norm, 'black', 'filled')
xlabel('Normalized Temperature')
ylabel('Normalized North Velocity')
set(gca, "FontSize", 18)
grid on
axis equal


dt = median(diff(Vtime)) * 86400;
dt = dt * 4;
T_ds = T(:,1:4:end);
N_ds = N(:,1:4:end);

L = length(Vtime)/4;
M = round(L/8);
window = hann(M);
noverlap = M/2;

K = (L-noverlap)/(M-noverlap);
C95 = 1 - 0.05^(1/(K-1)); 

[Cxy, f] = mscohere(T_ds(10,:), N_ds(10,:), window, noverlap, [], 1/dt);
% add period?
figure
coh = semilogx(f, Cxy, 'k', 'LineWidth',2);
hold on
yline(C95, 'r--', 'LineWidth', 1.5, 'Label', '$C_{95}$', 'Interpreter','latex', 'FontSize', 18)
grid on
%coh.XAxis.

