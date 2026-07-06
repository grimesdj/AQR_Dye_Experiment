%% Variance Bar Plot

clear all
close all



sites = 1:3;

% Example data (replace with yours)
load('../../../../Kelp_data/data/2024_PROCESSED_DATA/Velocity_Variance_North.mat');
North = variance;
clear variance 
load('../../../../Kelp_data/data/2024_PROCESSED_DATA/Velocity_Variance_East.mat');
East = variance;
figure('Position', [50 50, 500, 1000]); 
hold on

NV = sum(North, 2);
EV = sum(East, 2);
V = [EV NV];

subplot(3, 1, 1);
bar(sites, V, 'grouped')
legend('Alongshore Variance', 'Cross-Shore Variance', 'Location','northeastoutside')
xticks(sites)
xticklabels([])
% set(gca, 'Fontsize', 16)
title('Total Variance by Mooring')
grid on
ylabel('$\left[\frac{\mathrm{m}}{\mathrm{s}}\right]^2$', 'Interpreter','latex','Rotation', 0)

subplot(3, 1, 2);
bar(East, 'stacked')
legend('Barotropic Variance', 'Mode 1 BC Variance', 'BC Noise Variance', 'Location','northeastoutside')
xticks(sites)
xticklabels([])
% set(gca, 'Fontsize', 16)
title('Alongshore Variance')
grid on
ylabel('$\left[\frac{\mathrm{m}}{\mathrm{s}}\right]^2$', 'Interpreter','latex','Rotation', 0)

subplot(3, 1, 3);
bar(North, 'stacked')
legend('Barotropic Variance', 'Mode 1 BC Variance', 'BC Noise Variance', 'Location','northeastoutside')
xticks(sites)
xticklabels({'M1', 'M2', 'M3'})
% set(gca, 'Fontsize', 16)
title('Cross-Shore Variance')
grid on
ylabel('$\left[\frac{\mathrm{m}}{\mathrm{s}}\right]^2$', 'Interpreter','latex','Rotation', 0)

return


w = 0.35; % bar width offset

% Left stacked bar at each site
N = bar(sites - w/2, North, 'stacked', 'BarWidth', w);
C = vertcat(N.FaceColor);

% Right stacked bar at each site
E = bar(sites + w/2, East, 'stacked', 'BarWidth', w);
E(1).FaceColor = N(1).FaceColor;
E(2).FaceColor = N(2).FaceColor;
E(3).FaceColor = N(3).FaceColor;

xticks(sites)
xticklabels({'M1', 'M2', 'M3'})

legend('Barotropic Variance', 'Mode 1 Variance', 'Variance Due to Noise')
grid on



