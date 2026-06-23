%% EOF analysis of Temp data

clear all
close all

% initialize
FOV_fig = figure("Position", [2250 150 1000 800]);
t1 = tiledlayout(2, 2);

EOF_fig = figure("Position", [2250 150 1000 800]);
t2 = tiledlayout(2, 2);

%% Load Data
for mooring_ID = 1:4
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
%sig = var(Y(:));

% variance explained
FOV = L/sum(L);
figure(FOV_fig)
nexttile
plot(FOV, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','black')
ylabel('FOV')
xlabel('Mode #')
set(gca, 'FontSize', 18)
axis square
grid minor
title(sprintf('%s', mooring), 'FontSize', 18)
ylim([0 1])

% first two EOFs
for i = 1:size(EOFs, 2)
    if EOFs(end, i) < 0
        EOFs(:, i) = -1* EOFs(:,i);
        EC(:, i) = -1 * EC(:, i);
    end
end


figure(EOF_fig)
nexttile
plot(EOFs(:, 1), dz, 'k', 'LineWidth', 2)
hold on
plot(EOFs(:, 2), dz, 'r', 'LineWidth', 2)
plot(EOFs(:, 3), dz, 'k--', 'LineWidth', 2)
axis ij
axis square
ylabel('Depth [m]')
xlabel('$^\circ\mathrm{C}^2$', 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
title(sprintf('%s', mooring), 'FontSize', 18)


Data(mooring_ID).L = L;
Data(mooring_ID).EOFs = EOFs;
%Data(mooring_ID).sig = sig;
Data(mooring_ID).EC = EC;
Data(mooring_ID).Error = Error;
Data(mooring_ID).Skill = Skill;
Data(mooring_ID).lam = lam;
Data(mooring_ID).FOV = FOV;
Data(mooring_ID).Y = Y;


end
lgd = legend('1st Mode', '2nd Mode', '3rd Mode');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 3;