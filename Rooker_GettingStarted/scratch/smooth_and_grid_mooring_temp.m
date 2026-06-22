%% EOF analysis of Temp data

clear all
close all

%% Load Data
mooring_ID = 1;
moorings = {'M1', 'M2', 'M3', 'M4'};
mooring = moorings{mooring_ID};
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_10min_gridded.mat";
load(fullfile(fpath, savestr))

%% EOF
Y = Temp_grid';
[L, EOFs, EC, Error, Skill,lam] = EOF(Y);

% total variance
sig = var(Y(:));

% variance explained
FOV = L/sig;
figure
plot(FOV, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','black')

% first two EOFs
figure
plot(EOFs(:, 1), dz, 'k', 'LineWidth', 2)
hold on
plot(EOFs(:, 2), dz, 'r', 'LineWidth', 2)
axis ij
axis square
legend('1st Mode', '2nd Mode')
