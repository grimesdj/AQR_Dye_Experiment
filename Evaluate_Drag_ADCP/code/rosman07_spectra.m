%% rosman spectra

clear all
close all

% initalize psd
ds = 30;

w = 720 * 600/ds;
window = hamming(w);
noverlap = w/2;
nfft = [];
fs = 1/ds;

%% Load
for mooring_ID = 1:4
moorings = {'M1', 'M2', 'M3', 'M4'}; 
mooring = moorings{mooring_ID};

fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
fname = "mooring_" + mooring + ".mat";
M.(mooring) = load(fullfile(fpath, fname));
end

idx = round(length(M.M1.Temperature_mab)/2);

x = M.M1.Temperature(idx, 1:ds:end);
wpass = 1/30;
[y, ~, ~] = hamming_filter(x', wpass, fs);
T1 = y';
[pxx, f] = pwelch(T1, window, noverlap, nfft, fs);
loglog(f, movmean(pxx, 5), 'k', 'LineWidth', 2)
%grid on
axis square
% 
% prof = mean(M.M2.Temperature, 2);
% strat = diff(prof)./diff(M.M2.Temperature_mab);
% idx = find(strat == max(strat));
idx = round(length(M.M2.Temperature_mab)/2);

hold on
x = M.M2.Temperature(idx, 1:ds:end);
wpass = 1/30;
[y, ~, ~] = hamming_filter(x', wpass, fs);
T2 = y';
[pxx, f] = pwelch(T2, window, noverlap, nfft, fs);
loglog(f, movmean(pxx, 5), 'r', 'LineWidth', 2)
%grid on
axis square

% conf
N = length(T1);
S = w * (1-noverlap/w);
K = 1 + (N-w)/S;
dof = 2 * K;

Xp = chi2inv(0.9, dof);
Xm = chi2inv(0.1, dof);
err = (pxx(end) * [dof/Xp dof/Xm]);

loglog(f(2)*[1 1], err, 's-k', 'LineWidth', 2)

% bouyancy freq

 prof = mean(M.M2.Temperature, 2);
 strat = diff(prof)./diff(M.M2.Temperature_mab);
alpha = 2 * 10^(-4);
g = 9.81;

N2 = g * alpha * strat(idx);
N = sqrt(N2);
f = N/(2*pi);
xline(f, '--', 'LineWidth', 1)

legend('M3, Outside Kelp', 'M2, Inside Kelp', '90% Confidence', 'Buoyancy Freqency')
set(gca, 'FontSize', 18)
ylim([10e-3 10e5])

fpath = fullfile('..', '..', '..', '..', 'Documents', 'YCSECA', '2026', 'figures');
print(gcf, fullfile(fpath, '30s_Spectra.png'), '-dpng', '-r600')