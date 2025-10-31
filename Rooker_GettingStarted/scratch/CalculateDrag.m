%% Calculate Drag
clear all
close all

% Load data
filestem = '../../../../Kelp_data/Summer2025/Rooker/Release2/LPF';
files = dir([filestem, '/*PCA.mat']);

AQD     = load([files(1).folder, '/', files(1).name]);
M1      = load([files(1).folder, '/', files(2).name]);
M2      = load([files(1).folder, '/', files(3).name]);
VEC     = load([files(1).folder, '/', files(4).name]);


% for standarazing puropses
outV    = [M1.PCA_X M1.PCA_Y];

% Averaging M2 and VEC
meanx   = mean([M2.PCA_X, VEC.PCA_X(1:length(M2.PCA_X))], 2);
meany   = mean([M2.PCA_Y, VEC.PCA_Y(1:length(M2.PCA_Y))], 2);
inV     = [meanx meany];

% Find RMS
rmsout  = rms(outV, 1);
rmsin   = rms(inV, 1);

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
plot(outV(:,1), outV(:,2), 'r.', 'LineWidth', 2)

hold on
plot(inV(:,1), inV(:,2), 'b.', 'LineWidth', 2)

% plot error ellipse
t = 0:.01:2*pi;
ellipseout = [rmsout(:, 1) * cos(t); rmsout(:, 2) * sin(t)]';
ellipsein = [rmsin(:, 1) * cos(t); rmsin(:, 2) * sin(t)]';

plot(ellipseout(:,1), ellipseout(:,2), 'r--', 'LineWidth', 2)
plot(ellipsein(:,1), ellipsein(:,2), 'b--', 'LineWidth', 2)
legend('Outside Velocity', 'Inside Velocity', 'Outside RMS', 'Inside RMS')
title('Velocities Outside vs. Inside Kelp Forest', 'FontSize', 25)
xlabel('Velocity, u (m/s)', 'FontSize', 18)
ylabel('Velocity, v (m/s)', 'FontSize', 18)

exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/RMS.png')

