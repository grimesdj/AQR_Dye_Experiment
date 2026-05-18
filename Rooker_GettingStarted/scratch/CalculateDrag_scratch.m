%% Calculate Drag
clear all
close all

%% Load data
filestem = '../../../../Kelp_data/Summer2025/Rooker/Release1/LPF';
files = dir([filestem, '/*_PCA.mat']);


AQD1    = load([files(1).folder, '/', files(1).name]);
M11     = load([files(1).folder, '/', files(2).name]);
M21     = load([files(1).folder, '/', files(3).name]);
VEC1    = load([files(1).folder, '/', files(4).name]);

filestem = '../../../../Kelp_data/Summer2025/Rooker/Release2/LPF';
files = dir([filestem, '/*_PCA.mat']);

AQD2    = load([files(1).folder, '/', files(1).name]);
M12     = load([files(1).folder, '/', files(2).name]);
M22     = load([files(1).folder, '/', files(3).name]);
VEC2    = load([files(1).folder, '/', files(4).name]);

% for standarazing puropses
outV    = [M21.Velocity_X M21.Velocity_Y; M22.Velocity_X M22.Velocity_Y];
fprintf('Using East and North instead of PCA\n')

% Averaging M2 and VEC
fprintf('Seeing what this looks like without the avg\n')
meanx   = [VEC1.Velocity_X; VEC2.Velocity_X]; %mean([M2.PCA_X, VEC.PCA_X(1:length(M2.PCA_X))], 2);
meany   = [VEC1.Velocity_Y; VEC2.Velocity_Y];%mean([M2.PCA_Y, VEC.PCA_Y(1:length(M2.PCA_Y))], 2);
inV     = [meanx meany];

%% Find Representative Stats

% RMS
rmsout  = rms(outV, 1);
rmsin   = rms(inV, 1);
rms_redux = (1-(rmsin./rmsout))*100;
fprintf('RMS redux: %.3f%%\n', rms_redux)

% Absolute Average
absout  = mean(abs(outV), 1);
absin   = mean(abs(inV), 1);
abs_redux = (1-(absin./absout))*100;
fprintf('Abs redux: %.3f%%\n', abs_redux)


% Calculate Energy Lost
E_out   = (1/2).*1000.*(rmsout.^2);
E_in    = (1/2).*1000.*(rmsin.^2);
E_loss  = E_out-E_in;

% Get an ellipse together
Cout    = cov(outV);
muout   = mean(outV);
Cin     = cov(inV);
muin    = mean(inV);

% plot it!
figure
plot(outV(:,1), outV(:,2), 'r.', 'MarkerSize', 10)

hold on
plot(inV(:,1), inV(:,2), 'b.', 'MarkerSize', 10)

% plot error ellipse
t = 0:.01:2*pi;
ellipseout = [rmsout(:, 1) * cos(t); rmsout(:, 2) * sin(t)]';
ellipsein = [rmsin(:, 1) * cos(t); rmsin(:, 2) * sin(t)]';

plot(ellipseout(:,1), ellipseout(:,2), 'r--', 'LineWidth', 2)
plot(ellipsein(:,1), ellipsein(:,2), 'b--', 'LineWidth', 2)
legend('Outside Velocity', 'Inside Velocity', 'Outside RMS', 'Inside RMS')
title('Velocities Outside vs. Inside Kelp Forest', 'FontSize', 25)
xlabel('Velocity, East (m/s)', 'FontSize', 18)
ylabel('Velocity, N (m/s)', 'FontSize', 18)

%exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Release1_RMS.png')

% make scatters in vs out

% Allocate Vars
oU = outV(:,1);
iU = inV(:,1);

oV = outV(:,2);
iV = inV(:,2); 

% u scatter
ufig = figure;
scatter(oU(1:100:end), iU(1:100:end), 25,  'b', 'filled')
ylabel('Inside Velocity, $E_{in}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
xlabel('Outside Velocity, $E_{out}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
sgtitle('Alongshore Inside vs Outside Velocities', 'FontSize', 30)


% u best fit
xmin = min(oU) - 0.05;
xmax = max(oU) + 0.05;

x_ext = linspace(xmin, xmax, 100);

Pu = polyfit(oU, iU, 1); 
uline = polyval(Pu, x_ext);
hold on
[xs, idx] = sort(oU);
plot(x_ext, uline, 'k--', 'LineWidth', 2)
legend('$u$ Velocity Data', sprintf('Best Fit Slope = %.2f', Pu(1)), ...
    'Interpreter', 'latex', 'FontSize', 14, 'location', 'southeast');

grid on
axis equal
xlim(gca, [-0.25, 0.25])

Ucoef = corrcoef(oU, iU);
text(mean(oU) - 4*std(oU), mean(iU) + 4*std(iU), sprintf('r = %.3f', Ucoef(2)), 'FontSize', 20, 'EdgeColor', 'r')


%exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Release1_East_rdux.png')

% v scatter
vfig = figure;
scatter(oV(1:100:end), iV(1:100:end), 25, 'r', 'filled')
ylabel('Inside Velocity, $U_{in}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
xlabel('Outside Velocity, $U_{out}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
sgtitle('Cross-Shore Inside vs Outside Velocities', 'FontSize', 30)

% v best fit

% define extended x-range (e.g., 10% beyond data)
xmin = min(oV) - 0.05;
xmax = max(oV) + 0.05;

x_ext = linspace(xmin, xmax, 100);

Pv = polyfit(oV, iV, 1);
vline = polyval(Pv, x_ext);
hold on
[xs, idx] = sort(oV);
plot(x_ext, vline, 'k--', 'LineWidth', 2)
legend('$v$ Velocity Data', sprintf('Best Fit Slope = %.2f', Pv(1)), ...
    'Interpreter', 'latex', 'FontSize', 14, 'location', 'southeast');

grid on
axis equal
xlim(gca, [-0.25, 0.25])

% Slope Redux
m = [Pu(1) Pv(1)];
m_redux = (1-m)*100;
fprintf('Slope redux: %.3f%%\n', m_redux)

Vcoef = corrcoef(oV, iV);
text(mean(oV) - 4*std(oV), mean(iV) + 4*std(iV), sprintf('r = %.3f', Vcoef(2)), 'FontSize', 20, 'EdgeColor', 'r')

%exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Release1_North_rdux.png')

% Normalize Vars
norm_oU = (oU-mean(oU)) ./ std(oU);
norm_iU = (iU-mean(iU)) ./ std(iU);

norm_oV = (oV-mean(oV)) ./ std(oV);
norm_iV = (iV-mean(iV)) ./ std(iV);

% u scatter
unfig = figure;
scatter(norm_oU, norm_iU, 25,  'b', 'filled')
ylabel('Inside Velocity, $u_{in}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
xlabel('Outside Velocity, $u_{out}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
sgtitle('Normalized East Inside vs Outside Velocities', 'FontSize', 30)


% u best fit
hold on % why do about orth regressions return 1??????????
[coef, ~, lat] = pca([norm_oU norm_iU]);
a_orth = coef(2,1) / coef(1,1);
xx = linspace(-3, 3, 100);
yy = a_orth * xx;
plot(xx, yy, 'k--', 'LineWidth', 2)
legend('$u$ Velocity Data', sprintf('Best Fit Slope = %.2f', a_orth), ...
    'Interpreter', 'latex', 'FontSize', 14, 'location', 'southeast');
% u orthoganal fit



grid on
axis equal
xlim(gca, [-10, 10])


%exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Release1_norm_East_rdux.png')

% v scatter
vnfig = figure;
scatter(norm_oV, norm_iV, 25, 'r', 'filled')
ylabel('Inside Velocity, $v_{in}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
xlabel('Outside Velocity, $v_{out}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
sgtitle('Normalized North Inside vs Outside Velocities', 'FontSize', 30)

% v best fit
hold on
[coef, ~, lat] = pca([norm_oV norm_iV]);
a_orth = coef(2,1) / coef(1,1);
xx = linspace(-3, 3, 100);
yy = a_orth * xx;
plot(xx, yy, 'k--', 'LineWidth', 2)
legend('$v$ Velocity Data', sprintf('Best Fit Slope = %.2f', a_orth), ...
    'Interpreter', 'latex', 'FontSize', 14, 'location', 'southeast');
% v orthoganal fit


grid on
axis equal
xlim(gca, [-10, 10])


%exportgraphics(gcf, ['../../../../Kelp_data/Summer2025/Rooker/figures/Release1_norm_North_rdux.png'])

% Slope Redux
m = [Pu(1) Pv(1)];
m_redux = (1-m)*100;
fprintf('Slope redux: %.3f%%\n', m_redux)


% % u scatter
% figure;
% ufig = subplot(1,2,1);
% scatter(ufig, outV(:,1), inV(:,1), 25,  'b', 'filled')
% ylabel('Inside Velocity, $u_{in}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
% xlabel('Outside Velocity, $u_{out}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
% title('Inside vs Outside Velocities', 'FontSize', 30)
% 
% 
% % u best fit
% Pu = polyfit(outV(:,1), inV(:,1), 1); 
% uline = polyval(Pu, outV(:,1));
% hold on
% [xs, idx] = sort(outV(:,1));
% plot(ufig, xs, uline(idx), 'r--', 'LineWidth', 1)
% legend('$u$ Velocity Data', sprintf('Best Fit Slope = %.2f', Pu(1)), ...
%     'Interpreter', 'latex', 'FontSize', 14, 'location', 'southeast');
% axis equal
% 
% % v scatter
% vfig = subplot(1, 2, 2);
% scatter(vfig, outV(:,2), inV(:,2), 25, 'r', 'filled')
% ylabel('Inside Velocity, $v_{in}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
% xlabel('Outside Velocity, $v_{out}$ (m/s)', 'FontSize', 20, 'Interpreter','latex')
% 
% % v best fit
% Pv = polyfit(outV(:,2), inV(:,2), 1);
% vline = polyval(Pv, outV(:,2));
% hold on
% [xs, idx] = sort(outV(:,2));
% plot(vfig, xs, vline(idx), 'b--', 'LineWidth', 1)
% legend('$v$ Velocity Data', sprintf('Best Fit Slope = %.2f', Pv(1)), ...
%     'Interpreter', 'latex', 'FontSize', 14, 'location', 'southeast');
% axis equal
% 
% % Slope Redux
% m = [Pu(1) Pv(1)];
% m_redux = (1-m)*100;
% fprintf('Slope redux: %.3f%%\n', m_redux)
% 
% 
% exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Redux_Scatter.png')
