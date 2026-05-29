%% Plot Profiles

clear all
%close all

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
%[M1.Moor.Temperature_mab, idx] = sort(M1.Moor.Temperature_mab);
%M1.Moor.Temperature = M1.Moor.Temperature(idx, :);

% some errors with the loading so lets fix that for now

Temp_cont = M1.Moor.Temperature(:, M1.Moor.valid);
M1.Moor.Temperature_mab(1) = rms(M1.ADCP.maxRNG, 'all')-M1.Moor.Temperature_mab(1);
Temp_cont(6, :) = nan;



% plot
figure
img = imagesc(M1.ADCP.Time(M1.ADCP.valid), M1.ADCP.bin_mab, North(:,M1.LPF.valid));
set(img, 'AlphaData', ~isnan(North(:, M1.LPF.valid)))
set(gca, 'YDir', 'normal')
cb = colorbar;
colormap(cmocean('balance'))
clim([-0.1 0.1])
ylabel('Meters Above Bottom')
ylabel(cb, 'Cross-shore Velocity, $u$ [m/s]', 'Interpreter','latex', 'FontSize', 18)
set([gca cb], 'fontsize', 18)
ylim([0 max(M1.ADCP.maxRNG)+std(M1.ADCP.maxRNG)])

hold on
%contour(M1.Moor.Time(M1.Moor.valid), M1.Moor.Temperature_mab, M1.Moor.Temperature(:,M1.Moor.valid), [], 'k', 'LineWidth', 1.5)
contour(M1.Moor.Time(M1.Moor.valid),M1.Moor.Temperature_mab,  Temp_cont, [10 14 18 22], 'k', 'LineWidth', 0.5)
yline(1, 'k--', 'LineWidth', 1.5, 'label', 'Dye Release Depth', 'FontSize', 18, 'LabelHorizontalAlignment','left')

% hold on
% [C,h] = contour( ...
%     M1.Moor.Time(M1.Moor.valid), ...
%     M1.Moor.Temperature_mab(:), ...
%     M1.Moor.Temperature(:,M1.Moor.valid), ...
%     1, ...
%     'k', ...
%     'LineWidth',1);


 figure, imagesc(M1.Moor.Time(M1.Moor.valid), M1.Moor.Temperature_mab, M1.Moor.Temperature(:,M1.Moor.valid))