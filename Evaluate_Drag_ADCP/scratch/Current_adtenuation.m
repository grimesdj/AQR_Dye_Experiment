%% 

clear all
close all


%% Load
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M1', 'L1', 'ADCP');
fname = 'M1_10min_gridded_PCA.mat';
M1 = load(fullfile(fpath, fname));

%% alongshore

M1.East = M1.U_grid;
M1.East(M1.U_grid < 0) = NaN;
M1.West = M1.U_grid;
M1.West(M1.U_grid > 0) = NaN;

M1.East_std = std(M1.East, [], 2, 'omitnan');
M1.West_std = std(M1.West, [], 2, 'omitnan');

%% Load
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M3', 'L1', 'ADCP');
fname = 'M3_10min_gridded_PCA.mat';
M3 = load(fullfile(fpath, fname));

%% alongshore

M3.East = M3.U_grid;
M3.East(M3.U_grid < 0) = NaN;
M3.West = M3.U_grid;
M3.West(M3.U_grid > 0) = NaN;

M3.East_std = std(M3.East, [], 2, 'omitnan');
M3.West_std = std(M3.West, [], 2, 'omitnan');


figure
subplot(1, 2, 1)
plot(M1.East_std, M1.dz, 'r-s', 'LineWidth', 2)
hold on
plot(M3.East_std, M3.dz, 'b-s', 'LineWidth', 2)
title('East')
legend('M1','M3')
xline(0, 'k--', 'LineWidth', 1)
grid on


subplot(1, 2, 2)
plot(M1.West_std, M1.dz, 'r-s', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 1)
hold on
plot(M3.West_std, M3.dz, 'b-s', 'LineWidth', 2)
title('West')
grid on