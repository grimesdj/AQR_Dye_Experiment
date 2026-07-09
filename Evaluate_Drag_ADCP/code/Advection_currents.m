%% Advection Spectum

clear all
close all

%% Load Data

% files = dir('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/*.mat');
% for i = 1:length(files)
%     if ~contains(files(i).name, 'min') && ~contains(files(i).name, 'config')
%         fprintf('loading %s...\n', files(i).name)
%         dum = load(fullfile(files(i).folder, files(i).name));
%         if length(size(dum)) < 3
%             if i == 1
%                 M1 = dum;
%             else
%                 fields = fieldnames(dum);
%                 for num = 1:length(fields)
%                     field = fields{num};
% 
%                     M1.(field) = [M1.(field) dum.(field)];
%                 end
%             end
%         end
%     end
% end
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M1', 'L1', 'ADCP');
fname = 'M1_10min_gridded_PCA.mat';
M1 = load(fullfile(fpath, fname));

fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M2', 'L1', 'ADCP');
fname = 'M2_10min_gridded_PCA.mat';
M2 = load(fullfile(fpath, fname));

fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M3', 'L1', 'ADCP');
fname = 'M3_10min_gridded_PCA.mat';
M3 = load(fullfile(fpath, fname));

% 
% deployTime  =  datenum('03-Jul-2024 00:00:00');
% recoverTime  = datenum('25-Jul-2024 00:00:00');
% dep = find(M1.Time > deployTime & M1.Time < recoverTime);
% 



dt = round(mode(diff(M1.Tq(1, :))*86400));

% North Velocity Excursion
fprintf('Generating PSD\n')
x = mean(M1.V_grid, 1);
w = round(( 86400 * 10 ) / dt);
window = hamming(w);
noverlap = w/2;
nfft = [];
fs = 1/dt;

[pxx, f, pxxc] = pwelch(x, window, noverlap, nfft, fs, 'ConfidenceLevel', 0.95);

figure
subplot(2, 2, 1)
loglog(f, pxx, 'k-', 'LineWidth', 1.5)
%hold on
%loglog(f, pxxc, 'm--', 'LineWidth', 0.4)

xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18, 'LineWidth', 1, 'LabelVerticalAlignment','bottom')
xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1, 'LabelVerticalAlignment','bottom')
xlabel('Hz')
ylabel('$\left[\frac{\mathrm{m}^2}{\mathrm{s}}\right]$', 'Interpreter','latex', 'Rotation', 0)
grid on
set(gca, 'FontSize', 18)
title('Cross-Shore PSD')

% advection 
fprintf('Calculating Advection...\n')
df = mean(diff(f));

u_rms = sqrt(pxx * df);
x_rms = u_rms ./ (2*pi*f);

% Avoid the DC bin
x_rms(f==0) = NaN;

subplot(2, 2, 3)
loglog(f, x_rms, 'k-', 'LineWidth', 2)

%xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18,  'LineWidth', 1)
%xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1)
yline(25, 'g--', 'Label', 'Significant Excursion', 'FontSize', 18, 'LineWidth', 1)
xlabel('Hz')
ylabel('$L_x$ [m]', 'Interpreter','latex')
grid on
set(gca, 'FontSize', 18)

title('Cross-Shore Excursion Length Spectum')












%% alongshore

% 
% deployTime  =  datenum('03-Jul-2024 00:00:00');
% recoverTime  = datenum('25-Jul-2024 00:00:00');
% dep = find(M1.Time > deployTime & M1.Time < recoverTime);
% 



dt = round(mode(diff(M1.Tq(1, :))*86400));

% North Velocity Excursion
fprintf('Generating PSD\n')
x = mean(M1.U_grid, 1);
w = round(( 86400 * 10 ) / dt);
window = hamming(w);
noverlap = w/2;
nfft = [];
fs = 1/dt;

[pxx, f, pxxc] = pwelch(x, window, noverlap, nfft, fs, 'ConfidenceLevel', 0.95);

subplot(2, 2, 2)
loglog(f, pxx, 'k-', 'LineWidth', 1.5)
%hold on
%loglog(f, pxxc, 'm--', 'LineWidth', 0.4)

xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18, 'LineWidth', 1, 'LabelVerticalAlignment','bottom')
xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1, 'LabelVerticalAlignment','bottom')
xlabel('Hz')
ylabel('$\left[\frac{\mathrm{m}^2}{\mathrm{s}}\right]$', 'Interpreter','latex', 'Rotation', 0)
grid on
set(gca, 'FontSize', 18)
title('Alongshore PSD')

% advection 
fprintf('Calculating Advection...\n')
df = mean(diff(f));

u_rms = sqrt(pxx * df);
x_rms = u_rms ./ (2*pi*f);

% Avoid the DC bin
x_rms(f==0) = NaN;

subplot(2, 2, 4)
loglog(f, x_rms, 'k-', 'LineWidth', 2)

%xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18,  'LineWidth', 1)
%xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1)
yline(50, 'g--', 'Label', 'Significant Excursion', 'FontSize', 18, 'LineWidth', 1)
xlabel('Hz')
ylabel('$x$ [m]', 'Interpreter','latex')
grid on
set(gca, 'FontSize', 18)

title('Alongshore Excursion Length Spectum')


%% Atten
figure

x1 = mean(M1.U_grid(:, :), 1);
x3 = mean(M3.U_grid(:, 1:end-1), 1);

plot(x3, x1, 'k.', 'MarkerSize', 10)
xlabel('$\bar{U}_{\mathrm{M}3}$ [m/s]', 'Interpreter','latex')
ylabel('$\bar{U}_{\mathrm{M}1}$ [m/s]', 'Interpreter','latex')

X = [x3' x1'];
coeff = pca(X);
v = coeff(:, 1);
m = v(2)/v(1);
t = linspace(-0.2, 0.2);
y = m .* t;

hold on
plot(t, y, 'r--', 'LineWidth', 1.5)
axis square
ylim([-.3 .3])
xlim([-.3 .3])
grid minor
set(gca, 'FontSize', 18)
text(mean(t) + std(t), mean(y)-2*std(y), sprintf('m = %.2f', m), 'FontSize', 18, 'Color', 'red', 'EdgeColor', 'black', 'LineWidth', 1)

% try again but seperate east and west
west = x1 < 0 & x3 < 0;
WX = [x3(west)', x1(west)'];

east = x1 > 0 & x3 > 0;
EX = [x3(east)', x1(east)'];

figure
plot(x3, x1, 'k.', 'MarkerSize', 10)
xlabel('$\bar{U}_{\mathrm{M}3}$ [m/s]', 'Interpreter','latex')
ylabel('$\bar{U}_{\mathrm{M}1}$ [m/s]', 'Interpreter','latex')

coeff = pca(EX);
v = coeff(:, 1);
Em = v(2)/v(1);
Et = linspace(0, 0.2);
Ey = Et' .*v';


coeff = pca(WX);
v = coeff(:, 1);
Wm = v(2)/v(1);
Wt = linspace(-0.2, 0);
Wy = Wt' .*v';

hold on
plot(Ey(:,1), Ey(:,2), 'r--', 'LineWidth', 2)
plot(Wy(:, 1), Wy(:, 2), 'r--', 'LineWidth', 2)

axis square
ylim([-.3 .3])
xlim([-.3 .3])
grid minor
set(gca, 'FontSize', 18)
text(mean(t) + std(t), mean(y)-2*std(y), sprintf('Em = %.2f', Em), 'FontSize', 18, 'Color', 'red', 'EdgeColor', 'black', 'LineWidth', 1)
text(mean(t) + std(t), mean(y)-3*std(y), sprintf('Wm = %.2f', Wm), 'FontSize', 18, 'Color', 'red', 'EdgeColor', 'black', 'LineWidth', 1)





%% cross-shore
figure

x1 = M1.V_grid(end-1, :);
x2 = M2.V_grid(end-1, :);

plot(x1, x2, 'k.', 'MarkerSize', 10)
xlabel('$\bar{V}_{\mathrm{M}1}$ [m/s]', 'Interpreter','latex')
ylabel('$\bar{V}_{\mathrm{M}2}$ [m/s]', 'Interpreter','latex')

X = [x1' x2'];
coeff = pca(X);
v = coeff(:, 1);
m = v(2)/v(1);
t = linspace(-0.2, 0.2);
y = m .* t;

hold on
plot(t, y, 'r--', 'LineWidth', 1.5)
axis square
ylim([-.3 .3])
xlim([-.3 .3])
grid minor
set(gca, 'FontSize', 18)
text(mean(t) + std(t), mean(y)-2*std(y), sprintf('m = %.2f', m), 'FontSize', 18, 'Color', 'red', 'EdgeColor', 'black', 'LineWidth', 1)