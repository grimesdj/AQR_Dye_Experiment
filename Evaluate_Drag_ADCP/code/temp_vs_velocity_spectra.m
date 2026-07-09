% Temp vs Velocity Spectra

clear all
close all

%% Load Data

moorings = {'M1', 'M2', 'M3'};
for mooring_ID = 1:length(moorings)

% load temp
mooring = moorings{mooring_ID};
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_10min_gridded.mat";
Temp(mooring_ID) = load(fullfile(fpath, savestr));

% load temp EOF
savestr = mooring + "_EOF.mat";
TEOF = load(fullfile(fpath, savestr));

% ADCP
if mooring_ID == 3;
    % buh
else
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L0', 'ADCP');
fname = "ADCP_" + mooring + "_L0_10min.mat";
if exist(fullfile(fpath, fname), 'file')
    Vel(mooring_ID) = load(fullfile(fpath,fname));
else
   % Vel(mooring_ID) = Vel(mooring_ID-1);
end
end


% load vel EOF
savestr = mooring + "_EOF_depth_coords_U.mat";
path = fullfile(fpath, '..', '..', 'L1', 'ADCP');
    if exist(fullfile(path, savestr), 'file')
        VEOF = load(fullfile(path, savestr));
    end


%% Spectra

% M1

% Temp middle
A = TEOF.EC(:, 1);
phi = TEOF.EOFs(:, 1);
TBC = A * phi';
%BT = mean(Temp(mooring_ID).Temp_grid, 1);
num = length(phi);
idx = round(num/2);
TBC = TBC(:, idx);
%TBC = TBC - BT';
TBC = detrend(TBC);


psd_fig = figure;
subplot(2, 1, 1)
w = 720;
window = hamming(w);
noverlap = round(w*(2/3));
nfft = [];
fs = 1/600;
[pxx,f, pxxc] = pwelch(TBC,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
loglog(f(3:end), pxx(3:end), 'LineWidth', 1)
hold on
loglog(f(3:end), pxxc(3:end, :), 'k--', 'LineWidth', 0.5)
grid on
title(sprintf('Mode 1 Temperature %s', mooring))
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
set(gca, 'FontSize', 14)

% Velocity Bottom
A = VEOF.EC(:, 1);
phi = VEOF.EOFs(:,1);
VBC = A * phi';

VBC = VBC(:, end-1);

subplot(2, 1, 2)
[pxx,f, pxxc] = pwelch(VBC,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
loglog(f(3:end), pxx(3:end), 'LineWidth', 1)
hold on
loglog(f(3:end), pxxc(3:end, :), 'k--', 'LineWidth', 0.5)
grid on
title(sprintf('Mode 1 Velocity %s', mooring))
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
set(gca, 'FontSize', 14)

% Coherence

% match sizes i guess
id = 1:length(VBC);

x = VBC;
y = TBC(id);

[Cxy,F] = mscohere(x,y,window,noverlap,nfft,fs);

% coherence null hyp threshold
N = id(end);

step = w - noverlap;
numW = floor((N-noverlap)/step);
L = numW * (noverlap/w);
alpha = 0.05;
Ccrit = 1 - alpha^(1/(L-1));

% cospec

[Pxy,F] = cpsd(x,y,window,noverlap,nfft,fs);

% plot coherence
figure
ax1 = subplot(3,1,1);
semilogx(ax1, F,Cxy, 'LineWidth', 1)
title(sprintf('Mode 1 Magnitude-Squared Coherence %s', mooring))

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
yline(Ccrit, 'r--', 'Label', 'Significance Threshold')
set(gca, 'FontSize', 14)
ylim([0 1])

% phase
phase = angle(Pxy);
%phase(Cxy < Ccrit) = NaN;

ax2 = subplot(3,1,2);
semilogx(ax2, F, phase, 'LineWidth', 1)

xlabel('Hz')
ylabel('\Theta(f)')
title('Cross Spectrum Phase')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
yline([pi, pi/2, 0, -pi/2, -pi], 'k--', 'LineWidth', 0.4)
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'$-\pi$', '$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$','$\pi$'})
ax2.TickLabelInterpreter = "latex";
%ax2.YTickLabelRotation = 90;
set(gca, 'FontSize', 14)

ax3 = subplot(3,1,3);
loglog(ax3, F,abs(Pxy), 'LineWidth', 1)

xlabel('Hz')
title('Cross Spectrum magnitude')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
set(gca, 'FontSize', 14)


linkaxes([ax1 ax2 ax3], 'x')

%xlim([1/(40 * 3600) 1/(8 * 3600)])

%% Spectra of raw
% 
% Vel = Vel(1).Velocity_North;
% Temp = Temp(1).Temp_grid;
% 
% Vel = Vel(1, :); % bottom vel
% Temp = Temp(:, idx);

%% mode 2??

% M1

% Temp bottom
A = TEOF.EC(:, 2);
phi = TEOF.EOFs(:,2);
TBC = A * phi';
%BT = mean(Temp(1).Temp_grid, 1);
TBC = TBC(:, end-1);
%TBC = TBC - BT';
TBC = detrend(TBC);


psd_fig = figure;
subplot(2, 1, 1)
w = 720;
window = hamming(w);
noverlap = round(w*(2/3));
nfft = [];
fs = 1/600;
[pxx,f, pxxc] = pwelch(TBC,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
loglog(f(3:end), pxx(3:end), 'LineWidth', 1)
hold on
loglog(f(3:end), pxxc(3:end, :), 'k--', 'LineWidth', 0.5)
grid on
title(sprintf('Mode 2 Temperature %s', mooring))
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
set(gca, 'FontSize', 14)

% Velocity middle
A = VEOF.EC(:, 2);
phi = VEOF.EOFs(:,2);
VBC = A * phi';
num = length(phi);
idx = round(num/2);
VBC = VBC(:, idx);

subplot(2, 1, 2)
[pxx,f, pxxc] = pwelch(VBC,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
loglog(f(3:end), pxx(3:end), 'LineWidth', 1)
hold on
loglog(f(3:end), pxxc(3:end, :), 'k--', 'LineWidth', 0.5)
grid on
title(sprintf('Mode 2 Velocity %s', mooring))
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
set(gca, 'FontSize', 14)
% Coherence

% match sizes i guess
id = 1:length(VBC);

x = VBC;
y = TBC(id);

[Cxy,F] = mscohere(x,y,window,noverlap,nfft,fs);

% coherence null hyp threshold
N = id(end);

step = w - noverlap;
numW = floor((N-noverlap)/step);
L = numW * (noverlap/w);
alpha = 0.05;
Ccrit = 1 - alpha^(1/(L-1));

% cospec

[Pxy,F] = cpsd(x,y,window,noverlap,nfft,fs);

% plot coherence
figure
ax1 = subplot(3,1,1);
semilogx(ax1, F,Cxy, 'LineWidth', 1)
title(sprintf('Mode 2 Magnitude-Squared Coherence %s', mooring))

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
yline(Ccrit, 'r--', 'Label', 'Significance Threshold')
set(gca, 'FontSize', 14)
ylim([0 1])

% phase
phase = angle(Pxy);
%phase(Cxy < Ccrit) = NaN;

ax2 = subplot(3,1,2);
semilogx(ax2, F, phase, 'LineWidth', 1)

xlabel('Hz')
ylabel('\Theta(f)')
title('Cross Spectrum Phase')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
yline([pi, pi/2, 0, -pi/2, -pi], 'k--', 'LineWidth', 0.4)
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'$-\pi$', '$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$','$\pi$'})
ax2.TickLabelInterpreter = "latex";
set(gca, 'FontSize', 14)

ax3 = subplot(3,1,3);
loglog(ax3, F,abs(Pxy), 'LineWidth', 1)

xlabel('Hz')
title('Cross Spectrum magnitude')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
set(gca, 'FontSize', 14)


linkaxes([ax1 ax2 ax3], 'x')
%xlim([1/(40 * 3600) 1/(8 * 3600)])

end

%% Spectra of raw

Vel = Vel(1).Velocity_North;
Temp = Temp(1).Temp_grid;

Vel = Vel(1, :); % bottom vel
Temp = Temp(:, idx);