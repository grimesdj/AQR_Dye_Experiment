%% 10 min Stats

clear all
close all

AQD = load('../../../../Kelp_data/Summer2025/Rooker/Release2/LPF/AQD9788.mat');
M1 = load('../../../../Kelp_data/Summer2025/Rooker/Release2/LPF/M1SIG1K.mat');
M2 = load('../../../../Kelp_data/Summer2025/Rooker/Release2/LPF/M2SIG1K.mat');
VEC = load('../../../../Kelp_data/Summer2025/Rooker/Release2/LPF/VEC17861.mat');

%% Are we even using ENU???????


% Display Data
figure
ax1 = subplot(2, 2, 1);
scatter(ax1, AQD.Velocity_X, AQD.Velocity_Y, 'b.')
axis equal
title(AQD.label, 'FontSize', 18)

ax2 = subplot(2, 2, 2);
scatter(ax2, M1.Velocity_X, M1.Velocity_Y, 'r.')
axis equal
title(M1.label, 'FontSize', 18)

ax3 = subplot(2, 2, 3);
scatter(ax3, M2.Velocity_X, M2.Velocity_Y, 'm.')
axis equal
title(M2.label, 'FontSize', 18)

ax4 = subplot(2, 2, 4);
scatter(ax4, VEC.Velocity_X, VEC.Velocity_Y, 'g.')
axis equal
title(VEC.label, 'FontSize', 18)

% Calculate standard deviation
AQD.std = [std(AQD.Velocity_X) std(AQD.Velocity_Y)];
M1.std = [std(M1.Velocity_X) std(M1.Velocity_Y)];
M2.std = [std(M2.Velocity_X) std(M2.Velocity_Y)];
VEC.std = [std(VEC.Velocity_X) std(VEC.Velocity_Y)];

xlabel(ax1, {'\sigma = ' AQD.std(1)}, 'FontSize', 16)
ylabel(ax1, {'\sigma = ' AQD.std(2)}, 'FontSize', 16)

xlabel(ax2, {'\sigma = ' M1.std(1)}, 'FontSize', 16)
ylabel(ax2, {'\sigma = ' M1.std(2)}, 'FontSize', 16)

xlabel(ax3, {'\sigma = ' M2.std(1)}, 'FontSize', 16)
ylabel(ax3, {'\sigma = ' M2.std(2)}, 'FontSize', 16)

xlabel(ax4, {'\sigma = ' VEC.std(1)}, 'FontSize', 16)
ylabel(ax4, {'\sigma = ' VEC.std(2)}, 'FontSize', 16)

figure
hold on
scatter(AQD.Velocity_X, AQD.Velocity_Y, 'b.', 'LineWidth', 1)
scatter(M1.Velocity_X, M1.Velocity_Y, 'r.', 'LineWidth', 1)
scatter(M2.Velocity_X, M2.Velocity_Y, 'm.', 'LineWidth', 1)
scatter(VEC.Velocity_X, VEC.Velocity_Y, 'g.', 'LineWidth', 1)
axis equal
title('Data and \sigma of All Instruments', "FontSize", 20)

xline(AQD.std(1), 'b', 'Label', '\sigma AQD')
xline(-1*AQD.std(1), 'b', 'Label', '\sigma AQD')
yline(AQD.std(2), 'b', 'Label', '\sigma AQD')
yline(-1*AQD.std(2), 'b', 'Label', '\sigma AQD')

xline(M1.std(1), 'r', 'Label','\sigma M1')
xline(-1*M1.std(1), 'r', 'Label','\sigma M1')
yline(M1.std(2), 'r', 'Label','\sigma M1')
yline(-1*M1.std(2), 'r', 'Label','\sigma M1')

xline(M2.std(1), 'm', 'Label','\sigma M2')
xline(-1*M2.std(1), 'm', 'Label','\sigma M2')
yline(M2.std(2), 'm', 'Label','\sigma M2')
yline(-1*M2.std(2), 'm', 'Label','\sigma M2')

xline(VEC.std(1), 'g', 'Label','\sigma VEC')
xline(-1*VEC.std(1), 'g', 'Label','\sigma VEC')
yline(VEC.std(2), 'g', 'Label','\sigma VEC')
yline(-1*VEC.std(2), 'g', 'Label','\sigma VEC')

exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Compare_All_std.pdf')