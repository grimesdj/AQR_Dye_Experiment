%% bandpass Temp

clear all
close all

% intiallize psd
w = 720;
window = hamming(w);
noverlap = w/2;
nfft = [];
fs = 1/600;

spectra_fig = figure;
t1 = tiledlayout(2, 2);

std_fig = figure("Position", [2250 150 1000 800]);
t2 = tiledlayout(2, 2);

%% Load Data
moorings = {'M1', 'M2', 'M3', 'M4'};
for i = 1:length(moorings)
    mooring = moorings{i};
    
    fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
    savestr = mooring + "_EOF.mat";
    M.(mooring) = load(fullfile(fpath, savestr));


    %% Recontruct Mode 1
    
    phi = M.(mooring).EC(:, 1);
    A   = M.(mooring).EOFs(:, 1);
    S = phi * A';
    sig = std(S, [], 1);
    
    M.(mooring).U = S;
    M.(mooring).U_std = sig;

    figure(spectra_fig)
    nexttile
    [pxx, f] = pwelch(M.(mooring).U(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'k-', 'LineWidth', 2)
    hold on
    grid on
    title(sprintf('%s', mooring))

    %% Bandpass
    % Subtidal
    x = S;
    fpass = 1/(32 * 3600);
    fs = 1/600;
    M.(mooring).UST =  lowpass(x,fpass,fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    [pxx, f] = pwelch(M.(mooring).UST(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
    yl = ylim;
    xl = xlim;
    p1 = patch([xl(1) fpass fpass xl(1)], [ yl(1) yl(1) yl(2) yl(2)], [0.1 0.1 0.1], 'FaceAlpha', 0.2, 'edgecolor', 'none');
    
    
    % diurnal
    x = M.(mooring).U - M.(mooring).UST;
    fpass = [1/(32 * 3600) 1/(16 * 3600)];
    M.(mooring).UDU = bandpass(x,fpass,fs);
    [pxx, f] = pwelch(M.(mooring).UDU(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
    p2 = patch([fpass(1) fpass(2) fpass(2) fpass(1)], [yl(1) yl(1) yl(2) yl(2)], [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'edgecolor', 'none');
    x1 = xline(fpass, 'LineWidth', 1);
    tl = min(pxx);
    text(fpass(1)/1.5, tl, 'ST', 'FontSize', 20, 'HorizontalAlignment','center')
    text(mean(fpass)*0.95, tl, 'DU', 'FontSize', 20, 'HorizontalAlignment','center')

    % semi-diurnal
    x = M.(mooring).U - M.(mooring).UST - M.(mooring).UDU;
    fpass = [1/(16 * 3600) 1/(8 * 3600)];
    M.(mooring).USD = bandpass(x,fpass,fs);
    [pxx, f] = pwelch(M.(mooring).USD(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
    p3 = patch([fpass(1) fpass(2) fpass(2) fpass(1)], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'edgecolor', 'none');
    x2 = xline(fpass(2), 'LineWidth', 1);
    text(mean(fpass)*0.95, tl, 'SD', 'FontSize', 20, 'HorizontalAlignment','center')

    % Mid-high
    x = M.(mooring).U - M.(mooring).UST - M.(mooring).UDU;
    fpass = [1/(8 * 3600) 1/(1 * 3600)];
    M.(mooring).UMH = bandpass(x,fpass,fs);
    [pxx, f] = pwelch(M.(mooring).UMH(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
    p4 = patch([fpass(1) fpass(2) fpass(2) fpass(1)], [yl(1) yl(1) yl(2) yl(2)], [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'edgecolor', 'none');
    x3 = xline(fpass(2), 'LineWidth', 1);
    text(mean(fpass)*0.95, tl, 'MH', 'FontSize', 20, 'HorizontalAlignment','center')

    % high freq
    M.(mooring).UHF = M.(mooring).U - M.(mooring).UST - M.(mooring).UDU - M.(mooring).USD - M.(mooring).UMH;
    [pxx, f] = pwelch(M.(mooring).UHF(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
    text(1.5*fpass(2), tl, 'HF', 'FontSize', 20, 'HorizontalAlignment','center')
    set(gca, 'FontSize', 18)
    
    %% std profiles

    M.(mooring).UST_std = std(M.(mooring).UST, [], 1);
    M.(mooring).UDU_std = std(M.(mooring).UDU, [], 1);
    M.(mooring).USD_std = std(M.(mooring).USD, [], 1);
    M.(mooring).UMH_std = std(M.(mooring).UMH, [], 1);
    M.(mooring).UHF_std = std(M.(mooring).UHF, [], 1);

    %% compare std profiles at each band
    figure(std_fig);
    dum = M.(mooring);
    nexttile
    plot(dum.U_std, dum.dz, 'k-s', 'LineWidth', 2)
    hold on
    plot(dum.UST_std, dum.dz, '-s', 'LineWidth', 2)
    plot(dum.UDU_std, dum.dz, '-s', 'LineWidth', 2)
    plot(dum.USD_std, dum.dz, '-s', 'LineWidth', 2)
    plot(dum.UMH_std, dum.dz, '-s', 'LineWidth', 2)
    plot(dum.UHF_std, dum.dz, '-s', 'LineWidth', 2)
    axis square
    grid minor
    axis ij
    ylabel('$h$ [m]', 'Interpreter','latex')
    xlabel('$^{\circ}$C', 'Interpreter','latex')
    set(gca, 'FontSize', 18)
    title(sprintf('%s', mooring))
    
end
figure(spectra_fig)
uistack(p1, 'top')
uistack(p2, 'top')
uistack(p3, 'top')
uistack(p4, 'top')
uistack(x2, 'top')
uistack(x1, 'top')
uistack(x3, 'top')
lgd = legend('original', 'subtidal', 'diurnal', 'semi-diurnal', 'Mid-High', 'high-frequency');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 5;
lgd.AutoUpdate = "off";
uistack(p1, 'bottom')
uistack(p2, 'bottom')
uistack(p3, 'bottom')
uistack(p4, 'bottom')
uistack(x1, 'bottom')
uistack(x2, 'bottom')
uistack(x3, 'bottom')


figure(std_fig)
lgd = legend('original', 'subtidal', 'diurnal', 'semi-diurnal', 'Mid-High', 'high-frequency');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 5;
lgd.AutoUpdate = "off";

%% Original
figure
plot(M.M1.U_std, M.M1.dz, '-s', 'LineWidth', 2)
hold on
plot(M.M2.U_std, M.M2.dz, '-s', 'LineWidth', 2)
plot(M.M3.U_std, M.M3.dz, '-s', 'LineWidth', 2)
plot(M.M4.U_std, M.M4.dz, '-s', 'LineWidth', 2)
axis square
grid minor
axis ij
set(gca, 'FontSize', 18)
legend('M1', 'M2', 'M3', 'M4', 'Location', 'eastoutside')
title('original')

%% Subtidal
figure
plot(M.M1.UST_std, M.M1.dz, '-s', 'LineWidth', 2)
hold on
plot(M.M2.UST_std, M.M2.dz, '-s', 'LineWidth', 2)
plot(M.M3.UST_std, M.M3.dz, '-s', 'LineWidth', 2)
plot(M.M4.UST_std, M.M4.dz, '-s', 'LineWidth', 2)
axis square
grid minor
axis ij
set(gca, 'FontSize', 18)
legend('M1', 'M2', 'M3', 'M4', 'Location', 'eastoutside')
title('Subtidal')


%% Diurnal
figure
plot(M.M1.UDU_std, M.M1.dz, '-s', 'LineWidth', 2)
hold on
plot(M.M2.UDU_std, M.M2.dz, '-s', 'LineWidth', 2)
plot(M.M3.UDU_std, M.M3.dz, '-s', 'LineWidth', 2)
plot(M.M4.UDU_std, M.M4.dz, '-s', 'LineWidth', 2)
axis square
grid minor
axis ij
set(gca, 'FontSize', 18)
legend('M1', 'M2', 'M3', 'M4', 'Location', 'eastoutside')
title('diurnal')

%% Semidiurnal
figure
plot(M.M1.USD_std, M.M1.dz, '-s', 'LineWidth', 2)
hold on
plot(M.M2.USD_std, M.M2.dz, '-s', 'LineWidth', 2)
plot(M.M3.USD_std, M.M3.dz, '-s', 'LineWidth', 2)
plot(M.M4.USD_std, M.M4.dz, '-s', 'LineWidth', 2)
axis square
grid minor
axis ij
set(gca, 'FontSize', 18)
legend('M1', 'M2', 'M3', 'M4', 'Location', 'eastoutside')
title('semi-diurnal')

%% Mid-High
figure
plot(M.M1.UMH_std, M.M1.dz, '-s', 'LineWidth', 2)
hold on
plot(M.M2.UMH_std, M.M2.dz, '-s', 'LineWidth', 2)
plot(M.M3.UMH_std, M.M3.dz, '-s', 'LineWidth', 2)
plot(M.M4.UMH_std, M.M4.dz, '-s', 'LineWidth', 2)
axis square
grid minor
axis ij
set(gca, 'FontSize', 18)
legend('M1', 'M2', 'M3', 'M4', 'Location', 'eastoutside')
title('Mid-High')


%% HF
figure
plot(M.M1.UHF_std, M.M1.dz, '-s', 'LineWidth', 2)
hold on
plot(M.M2.UHF_std, M.M2.dz, '-s', 'LineWidth', 2)
plot(M.M3.UHF_std, M.M3.dz, '-s', 'LineWidth', 2)
plot(M.M4.UHF_std, M.M4.dz, '-s', 'LineWidth', 2)
axis square
grid minor
axis ij
set(gca, 'FontSize', 18)
legend('M1', 'M2', 'M3', 'M4', 'Location', 'eastoutside')
title('High-Frequency')



%% Compare Percent changes

% M1 is shortest
id = 1:length(M.M1.dz);

%% Original

% M4 --> M1
M4M1 = -(1 - M.M4.U_std(id) ./ M.M1.U_std) * 100;

% M4 --> M3
M4M3 = -(1 - M.M4.U_std(id) ./ M.M3.U_std(id)) *  100;

% M4 --> M2
M4M2 = -(1 - M.M4.U_std(id) ./ M.M2.U_std(id)) * 100;
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M3, dz, '-s', 'LineWidth', 2)
plot(M4M2, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
legend('M4 --> M1', 'M4 --> M3', 'M4 --> M2', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('% Change')
ylabel('$h$ [m]', 'Interpreter','latex')
title('original')

%% Subtidal

% M4 --> M1
M4M1 = -(1 - M.M4.UST_std(id) ./ M.M1.UST_std) * 100;

% M4 --> M3
M4M3 = -(1 - M.M4.UST_std(id) ./ M.M3.UST_std(id)) *  100;

% M4 --> M2
M4M2 = -(1 - M.M4.UST_std(id) ./ M.M2.UST_std(id)) * 100;
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M3, dz, '-s', 'LineWidth', 2)
plot(M4M2, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
legend('M4 --> M1', 'M4 --> M3', 'M4 --> M2', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('% Change')
ylabel('$h$ [m]', 'Interpreter','latex')
title('Subtidal')

%% Diurnal

% M4 --> M1
M4M1 = -(1 - M.M4.UDU_std(id) ./ M.M1.UDU_std) * 100;

% M4 --> M3
M4M3 = -(1 - M.M4.UDU_std(id) ./ M.M3.UDU_std(id)) *  100;

% M4 --> M2
M4M2 = -(1 - M.M4.UDU_std(id) ./ M.M2.UDU_std(id)) * 100;
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M3, dz, '-s', 'LineWidth', 2)
plot(M4M2, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
legend('M4 --> M1', 'M4 --> M3', 'M4 --> M2', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('% Change')
ylabel('$h$ [m]', 'Interpreter','latex')
title('Diurnal')

%% SD

% M4 --> M1
M4M1 = -(1 - M.M4.USD_std(id) ./ M.M1.USD_std) * 100;

% M4 --> M3
M4M3 = -(1 - M.M4.USD_std(id) ./ M.M3.USD_std(id)) *  100;

% M4 --> M2
M4M2 = -(1 - M.M4.USD_std(id) ./ M.M2.USD_std(id)) * 100;
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M3, dz, '-s', 'LineWidth', 2)
plot(M4M2, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
legend('M4 --> M1', 'M4 --> M3', 'M4 --> M2', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('% Change')
ylabel('$h$ [m]', 'Interpreter','latex')
title('semidiurnal')

%% MH

% M4 --> M1
M4M1 = -(1 - M.M4.UMH_std(id) ./ M.M1.UMH_std) * 100;

% M4 --> M3
M4M3 = -(1 - M.M4.UMH_std(id) ./ M.M3.UMH_std(id)) *  100;

% M4 --> M2
M4M2 = -(1 - M.M4.UMH_std(id) ./ M.M2.UMH_std(id)) * 100;
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M3, dz, '-s', 'LineWidth', 2)
plot(M4M2, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
legend('M4 --> M1', 'M4 --> M3', 'M4 --> M2', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('% Change')
ylabel('$h$ [m]', 'Interpreter','latex')
title('Mid-High')

%% HF

% M4 --> M1
M4M1 = -(1 - M.M4.UHF_std(id) ./ M.M1.UHF_std) * 100;

% M4 --> M3
M4M3 = -(1 - M.M4.UHF_std(id) ./ M.M3.UHF_std(id)) *  100;

% M4 --> M2
M4M2 = -(1 - M.M4.UHF_std(id) ./ M.M2.UHF_std(id)) * 100;
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M3, dz, '-s', 'LineWidth', 2)
plot(M4M2, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
legend('M4 --> M1', 'M4 --> M3', 'M4 --> M2', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('% Change')
ylabel('$h$ [m]', 'Interpreter','latex')
title('High-Frequency')