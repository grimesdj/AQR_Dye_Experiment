






%% DOESNT WORK -- delete













%% rosman spectra

clear all
close all

% initalize psd
ds = 30;

w = 1024;
window = hamming(w);
noverlap = w/2;
nfft = [];
fs = 1/600;

%% Load
for mooring_ID = 1:4
moorings = {'M1', 'M2', 'M3', 'M4'}; 
mooring = moorings{mooring_ID};

fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
fname = mooring + "_EOF.mat";
M.(mooring) = load(fullfile(fpath, fname));
end

% prof = std(M.M3.Temperature');
% idx = find(prof==max(prof));
idx = ceil(length(M.M3.dz)./2);
id = [idx-1 idx idx+1];

% x = M.M3.Temperature(id, 1:ds:end);
% x = mean(x);
x = M.M3.EC(:, 1);
wpass = 1/ds;
%[y, ~, ~] = hamming_filter(x, wpass, fs);
T1 = x;
[pxx, f] = pwelch(T1, window, noverlap, nfft, fs);
loglog(f, movmean(pxx, 5), 'k', 'LineWidth', 2)
%grid on
axis square
% 
% prof = std(M.M2.Temperature');
% idx = find(prof==max(prof));
idx = ceil(length(M.M2.dz)./2);
id = [idx-1 idx idx+1];


hold on
% x = M.M2.Temperature(id, 1:ds:end);
% x = mean(x);
x = M.M2.EC(:, 1);
wpass = 1/ds;
%[y, ~, ~] = hamming_filter(x, wpass, fs);
T2 = x;
[pxx, f] = pwelch(T2, window, noverlap, nfft, fs);
loglog(f, movmean(pxx, 5), 'r', 'LineWidth', 2)
grid on
axis square
ylabel('$\left[\frac{^{\circ}\mathrm{C}}{\mathrm{Hz}}\right]$', 'Rotation', 0, 'Interpreter','latex')
xlabel('$f$ [Hz]', 'Interpreter','latex')

% conf
N = length(T1);
S = w * (1-noverlap/w);
K = 1 + (N-w)/S;
dof = 1.5 * K;
dof = dof * 2.5;

Xp = chi2inv(0.9, dof);
Xm = chi2inv(0.1, dof);
err = (pxx(end-1) * [dof/Xp dof/Xm]) * 10000;

loglog(f(3)*[1 1], err, 's-k', 'LineWidth', 2)

% bouyancy freq
%prof = mean(M.M2.Temperature, 2);
%strat = diff(prof)./diff(M.M2.dz0.3160
alpha = 2 * 10^(-4);
g = 9.81;

N2 = g * alpha * 0.3160; % max strat
N = sqrt(N2);
f = N/(2*pi);
xline(f, '--', 'LineWidth', 1)

% coriolis freq
% omega = 7.2921e-5;
% O = omega/(2*pi);
% Cf = 2 * O * sind(34.5); % approximate lat
% xline(Cf, '-.', 'LineWidth', 1)

lgd = legend('M3, Outside Kelp', 'M2, Inside Kelp', '90% Confidence', 'Buoyancy Freqency', 'Coriolis Frequency');
lgd.Location = 'bestoutside';
set(gca, 'FontSize', 18)
ylim([10e-3 10e5])

fpath = fullfile('..', '..', '..', '..', 'Documents', 'YCSECA', '2026', 'figures');
%print(gcf, fullfile(fpath, '30s_Spectra.png'), '-dpng', '-r600')