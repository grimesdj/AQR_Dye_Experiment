clear all; 
close all;

%% File Setup
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/*.mat');
colors = {[0, 0, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0]};

for i = 1:length(files)
    dataCell{i} = load(fullfile(files(i).folder, files(i).name));
    labels{i} = dataCell{i}.Config.SN;
    
    % Determine appropriate bin
    if size(dataCell{i}.Velocity_East, 2) == 1
        dataCell{i}.bin = 1;
    else
        dataCell{i}.bin = round((1.237 - dataCell{i}.Config.blank) / dataCell{i}.Config.binSize);
    
    end
end


%% Averaged Velocities
addpath('../code')
[T_all, U_all, V_all] = LPF(dataCell, labels, colors);

% normalizing vector size
[~, sizeLimit] = size(U_all{2}); % gonna make this a max function later

AQD = [U_all{1}(1, 1:sizeLimit); V_all{1}(1, 1:sizeLimit)];
M1  = [U_all{2}(1, 1:sizeLimit); V_all{2}(1, 1:sizeLimit)];
M2  = [U_all{3}(1, 1:sizeLimit); V_all{3}(1, 1:sizeLimit)];
VEC = [U_all{4}(1, 1:sizeLimit); V_all{4}(1, 1:sizeLimit)];


% Compare inside instruments to M1 (Inst - M1)
AQDvsM1 = AQD-M1;
M2vsM1  = M2-M1;
VECvsM1 = VEC-M1;




%% Plot difference scatter in-out
figure
scatter(AQDvsM1(1,:), AQDvsM1(2,:), 'b.')
hold on
scatter(M2vsM1(1,:),M2vsM1(2,:) , 'm.')
scatter(VECvsM1(1,:),VECvsM1(2,:) , 'g.')

% Generate and plot trendlines
AQDtrend = polyfit(AQDvsM1(1,:), AQDvsM1(2,:), 1);
yfit = polyval(AQDtrend, AQDvsM1(1,:));
plot(AQDvsM1(1,:), yfit, 'b-', 'LineWidth', 2)

M2trend = polyfit(M2vsM1(1,:), M2vsM1(2,:), 1);
yfit = polyval(M2trend, M2vsM1(1,:));
plot(M2vsM1(1,:), yfit, 'm-', 'LineWidth', 2)

VECtrend = polyfit(VECvsM1(1,:), VECvsM1(2,:), 1);
yfit = polyval(VECtrend, VECvsM1(1,:));
plot(VECvsM1(1,:), yfit, 'g-', 'LineWidth', 2)

slope = mean([AQDtrend M2trend VECtrend]);
theta = atan2d(slope, 1);
sprintf('Angle of max difference variability is %f degrees', theta)


% Formatting and such
title("Inside Instruments vs M1 ADCP ( X - M1 )", "FontSize", 18)
legend(labels{[1 3 4]})
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)

%% Compare AQD - VEC
AQDvsVEC = AQD-VEC;
figure, scatter(AQDvsVEC(1,:), AQDvsVEC(2,:), '.')

% mean and covariance
mu = [mean(AQDvsVEC(1,:)), mean(AQDvsVEC(2,:))];
Sigma = cov(AQDvsVEC(1,:), AQDvsVEC(2,:));

% eigen decomposition
[eigvec,eigval] = eig(Sigma);

% parametric ellipse
theta = linspace(0,2*pi,200);
ellipse = [cos(theta); sin(theta)]';

% scale by sqrt of eigenvalues (std dev)
ellipse = ellipse * sqrt(eigval) * eigvec';

% shift to mean
ellipse = ellipse + mu;
hold on
plot(ellipse(:,1), ellipse(:,2), 'r','LineWidth',2)

% Formatting
title([labels{1} ' vs ' labels{4}], "FontSize", 18)
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)
legend('AQD vs  VEC', '1 \sigma')

%% Compare AQD - M2
AQDvsM2 = AQD-M2;
figure
scatter(AQDvsM2(1,:), AQDvsM2(2,:), '.')

% mean and covariance
mu = [mean(AQDvsM2(1,:)), mean(AQDvsM2(2,:))];
Sigma = cov(AQDvsM2(1,:), AQDvsM2(2,:));

% eigen decomposition
[eigvec,eigval] = eig(Sigma);

% parametric ellipse
theta = linspace(0,2*pi,200);
ellipse = [cos(theta); sin(theta)]';

% scale by sqrt of eigenvalues (std dev)
ellipse = ellipse * sqrt(eigval) * eigvec';

% shift to mean
ellipse = ellipse + mu;
hold on
plot(ellipse(:,1), ellipse(:,2), 'r','LineWidth',2)

% Formatting
title([labels{1} ' vs ' labels{3}], "FontSize", 18)
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)
legend('AQD vs M2', '1 \sigma')

%% Compare VEC - M2
VECvsM2 = VEC-M2;
figure
scatter(VECvsM2(1,:), VECvsM2(2,:), '.')

% mean and covariance
mu = [mean(VECvsM2(1,:)), mean(VECvsM2(2,:))];
Sigma = cov(VECvsM2(1,:), VECvsM2(2,:));

% eigen decomposition
[eigvec,eigval] = eig(Sigma);

% parametric ellipse
theta = linspace(0,2*pi,200);
ellipse = [cos(theta); sin(theta)]';

% scale by sqrt of eigenvalues (std dev)
ellipse = ellipse * sqrt(eigval) * eigvec';

% shift to mean
ellipse = ellipse + mu;
hold on
plot(ellipse(:,1), ellipse(:,2), 'r','LineWidth',2)

% Formatting
title([labels{4} ' vs ' labels{3}], "FontSize", 18)
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)
legend('VEC vs M2', '1 \sigma')

% figure(3)
% exportgraphics(figure(1), '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/timeseries/10_min_avg.png')
exportgraphics(figure(2), '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/QuiverPlot.png')
% exportgraphics(figure(3), '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/DifferenceScatter.png')
% exportgraphics(figure(4), '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/AQDvsVEC.png')
% exportgraphics(figure(5), '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/AQDvsM2.png')
% exportgraphics(figure(6), '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/VECvsM2.png')


filestem = '../../../../Kelp_data/Summer2025/Rooker/Release2/LPF';

labels = replace(labels, ' ', '');

for i = 1:length(labels)
    
    Velocity_X = U_all{i}';
    Velocity_Y = V_all{i}';
    Time = datenum(T_all{i}');
    label = labels{i};
    sprintf('Saving %s...', labels{i})
    save([filestem '/' labels{i} '.mat'], ...
         'Velocity_X', 'Velocity_Y', 'Time', 'label');
    pca_function(filestem, labels{i})
    title(labels{i})
end

