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
hold on
loglog(f, pxxc, 'm--', 'LineWidth', 0.4)

xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18, 'LineWidth', 1)
xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1)
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

xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18,  'LineWidth', 1)
xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1)
yline(25, 'g--', 'Label', 'Significant Excursion', 'FontSize', 18, 'LineWidth', 1)
xlabel('Hz')
ylabel('$x$ [m]', 'Interpreter','latex')
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
hold on
loglog(f, pxxc, 'm--', 'LineWidth', 0.4)

xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18, 'LineWidth', 1)
xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1)
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

xline(1/86400, 'b--', 'Label', 'Diurnal', 'FontSize', 18,  'LineWidth', 1)
xline(2/86400, 'b--', 'Label', 'Semi-Diurnal', 'FontSize', 18, 'LineWidth', 1)
yline(50, 'g--', 'Label', 'Significant Excursion', 'FontSize', 18, 'LineWidth', 1)
xlabel('Hz')
ylabel('$x$ [m]', 'Interpreter','latex')
grid on
set(gca, 'FontSize', 18)

title('Alongshore Excursion Length Spectum')


