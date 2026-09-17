%% Variance Bar Plot

clear all
close all

sites = 1:3;


load('../../../../Kelp_data/data/2024_PROCESSED_DATA/Velocity_Variance.mat');
North = variance_V;
East = variance_U;
figure('Position', [50 50, 800, 1200]); 
%figure
%hold on

NV = sum(North, 2);
EV = sum(East, 2);
V = [EV NV];

subplot(3, 1, 1);
b = bar(sites, V, 'grouped');
b(2).FaceColor = [0 0 0];
b(1).FaceColor = [0.9 0.1 0];
legend('Alongshore Variance', 'Cross-Shore Variance', 'Location','northeastoutside')
xticks(sites)
xticklabels([])
set(gca, 'Fontsize', 18)
title('Total Variance', 'FontSize', 18)
grid on
axis square
ylabel('$\left[\frac{\mathrm{m}}{\mathrm{s}}\right]^2$', 'Interpreter','latex','Rotation', 0)
%yl = ylim;
xlim([0 4])

subplot(3, 1, 2);
b = bar(East, 'stacked');
cmap = cmocean('haline');
step = floor(size(cmap, 1)/3);
c = cmap(1:step:end, :);
b(1).FaceColor = c(1, :);
b(2).FaceColor = c(2, :);
b(3).FaceColor = c(3, :);
legend('Barotropic Variance', 'Mode 1 BC Variance', 'BC Noise Variance', 'Location','northeastoutside')
xticks(sites)
xticklabels([])
set(gca, 'Fontsize', 18)
title('Alongshore Variance', 'FontSize', 18)
grid on
axis square
ylabel('$\left[\frac{\mathrm{m}}{\mathrm{s}}\right]^2$', 'Interpreter','latex','Rotation', 0)
%ylim(yl)

subplot(3, 1, 3);
b = bar(North, 'stacked');
b(1).FaceColor = c(1, :);
b(2).FaceColor = c(2, :);
b(3).FaceColor = c(3, :);
legend('Barotropic Variance', 'Mode 1 BC Variance', 'BC Noise Variance', 'Location','northeastoutside')
xticks(sites)
xticklabels({'M1', 'M2', 'M3'})
set(gca, 'Fontsize', 18)
title('Cross-Shore Variance', 'FontSize', 18)
grid on
axis square
ylabel('$\left[\frac{\mathrm{m}}{\mathrm{s}}\right]^2$', 'Interpreter','latex','Rotation', 0)
%ylim(yl)
yticks([0:0.005:0.01])

%%export for poster
fpath = fullfile('..', '..', '..', '..', 'Documents', 'YCSECA', '2026', 'figures');
print(gcf, fullfile(fpath, 'variance_breakdown.png'), '-dpng', '-r600')
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



