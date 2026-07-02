% Temp vs Velocity Spectra

clear all
close all

%% Load Data

moorings = {'M1', 'M2', 'M3', 'M4'};
for mooring_ID = 1:length(moorings)

% load temp
mooring = moorings{mooring_ID};
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_10min_gridded.mat";
Temp(mooring_ID) = load(fullfile(fpath, savestr));

% load temp EOF
savestr = mooring + "_EOF.mat";
TEOF(mooring_ID) = load(fullfile(fpath, savestr));

% ADCP
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L0', 'ADCP');
fname = "ADCP_" + mooring + "_L0_10min.mat";
% if exist(fullfile(fpath, fname), 'file')
%     Vel(mooring_ID) = load(fullfile(fpath,fname));
% else
%     Vel(mooring_ID) = Vel(mooring_ID-1);
% end

% load vel EOF
savestr = mooring + "_EOF_depth_coords.mat";
path = fullfile(fpath, '..', '..', 'L1', 'ADCP');
    if exist(fullfile(path, savestr), 'file')
        VEOF(mooring_ID) = load(fullfile(path, savestr));
    end

end

%% Spectra

% M1

% Temp middle
A = TEOF(1).EC(:, 1);
phi = TEOF(1).EOFs(:,1);
TBC = A * phi';
BT = mean(Temp(1).Temp_grid, 1);
num = length(phi);
idx = round(num/2);
TBC = TBC(:, idx);
TBC = TBC - BT';
TBC = detrend(TBC);


figure
w = 720;
window = hamming(w);
noverlap = round(w*(2/3));
nfft = [];
fs = 1/600;
[pxx,f, pxxc] = pwelch(TBC,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
semilogx(f(3:end), pxx(3:end), 'LineWidth', 1)
hold on
semilogx(f(3:end), pxxc(3:end, :), 'k--', 'LineWidth', 0.5)
grid on
title('Baroclinic Temperature M1')
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')

% Velocity Bottom
A = VEOF(1).EC(:, 1);
phi = VEOF(1).EOFs(:,1);
VBC = A * phi';

VBC = VBC(:, end-1);

figure
[pxx,f, pxxc] = pwelch(VBC,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
semilogx(f(3:end), pxx(3:end), 'LineWidth', 1)
hold on
semilogx(f(3:end), pxxc(3:end, :), 'k--', 'LineWidth', 0.5)
grid on
title('Baroclinic Velocity M1')
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')

% Coherence

% match sizes i guess
idx = 1:length(VBC);

x = VBC;
y = TBC(idx);

[Cxy,F] = mscohere(x,y,[],[],[],fs);

[Pxy,F] = cpsd(x,y,[],[],[],fs);

figure
subplot(3,1,1)
semilogx(F,Cxy)
title('Magnitude-Squared Coherence')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')

subplot(3,1,2)
semilogx(F,angle(Pxy))

xlabel('Hz')
ylabel('\Theta(f)')
title('Cross Spectrum Phase')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
yline([0, -pi/2, -pi], 'r--', 'LineWidth', 0.4)
yticks([-pi -pi/2 0])

subplot(3,1,3)
semilogx(F,abs(Pxy))

xlabel('Hz')
ylabel('\Theta(f)')
title('Cross Spectrum magnitude')

xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')



% 
% [pxy, f] = cpsd(TBC(idx), VBC, window, noverlap, nfft, fs);
% figure
% semilogx(f, pxy, 'k-', 'LineWidth', 2)
% grid on
% %ylim([0 1])
% xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
% xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
% xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
