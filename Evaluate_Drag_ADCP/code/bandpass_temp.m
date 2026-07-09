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

T = [48 18 8 2];
fp = 1./(T * 3600);

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
    fpass = fp(1);
    fs = 1/600;
    M.(mooring).UST =  lowpass(x,fpass,fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    [pxx, f] = pwelch(M.(mooring).UST(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
   
    
    % diurnal
    x = M.(mooring).U - M.(mooring).UST;
    fpass = fp(2);
    M.(mooring).UDU = lowpass(x,fpass,fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    [pxx, f] = pwelch(M.(mooring).UDU(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)

    % semi-diurnal
    x = M.(mooring).U - M.(mooring).UST - M.(mooring).UDU;
    fpass = fp(3);
    M.(mooring).USD = lowpass(x,fpass,fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    [pxx, f] = pwelch(M.(mooring).USD(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)

    % Mid-high
    x = M.(mooring).U - M.(mooring).UST - M.(mooring).UDU - M.(mooring).USD;
    fpass = fp(4);
    M.(mooring).UMH = lowpass(x,fpass,fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    [pxx, f] = pwelch(M.(mooring).UMH(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)

    % high freq
    M.(mooring).UHF = M.(mooring).U - M.(mooring).UST - M.(mooring).UDU - M.(mooring).USD - M.(mooring).UMH;
    [pxx, f] = pwelch(M.(mooring).UHF(:, 1), window, noverlap, nfft, fs);
    loglog(f, pxx, 'LineWidth', 1)
    set(gca, 'FontSize', 18)
    


    yl = ylim;
    yl(1) = yl(1)/100;
    yl(2) = yl(2)*10;
    xl = xlim;
    p1 = patch([xl(1) fp(1) fp(1) xl(1)], [ yl(1) yl(1) yl(2) yl(2)], [0 0 0], 'FaceAlpha', 0.4, 'edgecolor', 'none');
    p2 = patch([fp(1) fp(2) fp(2) fp(1)], [ yl(1) yl(1) yl(2) yl(2)], [0 0 0], 'FaceAlpha', 0.3, 'edgecolor', 'none');
    p3 = patch([fp(2) fp(3) fp(3) fp(2)], [ yl(1) yl(1) yl(2) yl(2)], [0 0 0], 'FaceAlpha', 0.2, 'edgecolor', 'none');
    p4 = patch([fp(3) fp(4) fp(4) fp(3)], [ yl(1) yl(1) yl(2) yl(2)], [0 0 0], 'FaceAlpha', 0.1, 'edgecolor', 'none');
    text(sqrt(f(2) * fp(1)), yl(1), 'ST', 'FontSize', 20, 'HorizontalAlignment','center', 'VerticalAlignment', 'bottom')
    text(sqrt(fp(1)*fp(2)), yl(1), 'DU', 'FontSize', 20, 'HorizontalAlignment','center', 'VerticalAlignment', 'bottom')
    text(sqrt(fp(2)*fp(3)), yl(1), 'SD', 'FontSize', 20, 'HorizontalAlignment','center', 'VerticalAlignment', 'bottom')
    text(sqrt(fp(3)*fp(4)), yl(1), 'MH', 'FontSize', 20, 'HorizontalAlignment','center', 'VerticalAlignment', 'bottom')
    text(sqrt(fp(4) * f(end)), yl(1), 'HF', 'FontSize', 20, 'HorizontalAlignment','center', 'VerticalAlignment', 'bottom')
    x1 = xline(fp, 'LineWidth', 1);
    ylim(yl)
    xlim(xl)

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
uistack(x1, 'top')

lgd = legend('original', 'subtidal', 'diurnal', 'semi-diurnal', 'Mid-High', 'high-frequency');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 5;
lgd.AutoUpdate = "off";
p = findall(gcf, 'Type', 'patch');
l = findall(gcf, 'Type', 'ConstantLine');
for i = 1:length(p)
    uistack(l(i), 'bottom')
    uistack(p(i), 'bottom')
end
ax = findall(gcf, 'Type', 'axes');
ylabel(ax, '$\left[\frac{^{\circ}\mathrm{C}^2}{\mathrm{s}}\right]$', 'Interpreter','latex', 'Rotation', 0)
xlabel(ax, '$f$ [Hz]', 'Interpreter','latex')


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
id = 0:length(M.M1.dz)-1;

%% Original

% M4 --> M1
M4M1 = M.M1.U_std(end-id) - M.M4.U_std(end-id);

% M4 --> M3
M4M3 = M.M3.U_std(end-id) - M.M4.U_std(end-id);

% M4 --> M2
M4M2 = M.M2.U_std(end-id) - M.M4.U_std(end-id);
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M2, dz, '-s', 'LineWidth', 2)
plot(M4M3, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
xline(0, 'k--')
xl = xlim;
xm = max(abs(xl));
xlim([-xm xm])
legend('M4 --> M1', 'M4 --> M2', 'M4 --> M3', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('$\Delta\sigma_T$', 'Interpreter', 'latex')
ylabel('$h$ [m]', 'Interpreter','latex')
title('original')

%% Subtidal

% M4 --> M1
M4M1 = M.M1.UST_std(end-id) - M.M4.UST_std(end-id);

% M4 --> M3
M4M3 = M.M3.UST_std(end-id) - M.M4.UST_std(end-id);

% M4 --> M2
M4M2 = M.M2.UST_std(end-id) - M.M4.UST_std(end-id);
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M2, dz, '-s', 'LineWidth', 2)
plot(M4M3, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
xline(0, 'k--')
xl = xlim;
xm = max(abs(xl));
xlim([-xm xm])
legend('M4 --> M1', 'M4 --> M2', 'M4 --> M3', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('$\Delta\sigma_T$', 'Interpreter', 'latex')
ylabel('$h$ [m]', 'Interpreter','latex')
title('Subtidal')

%% Diurnal

% M4 --> M1
M4M1 = M.M1.UDU_std(end-id) - M.M4.UDU_std(end-id);

% M4 --> M3
M4M3 = M.M3.UDU_std(end-id) - M.M4.UDU_std(end-id);

% M4 --> M2
M4M2 = M.M2.UDU_std(end-id) - M.M4.UDU_std(end-id);
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M2, dz, '-s', 'LineWidth', 2)
plot(M4M3, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
xline(0, 'k--')
xl = xlim;
xm = max(abs(xl));
xlim([-xm xm])
legend('M4 --> M1', 'M4 --> M2', 'M4 --> M3', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('$\Delta\sigma_T$', 'Interpreter', 'latex')
ylabel('$h$ [m]', 'Interpreter','latex')
title('Diurnal')

%% SD

% M4 --> M1
M4M1 = M.M1.USD_std(end-id) - M.M4.USD_std(end-id);

% M4 --> M3
M4M3 = M.M3.USD_std(end-id) - M.M4.USD_std(end-id);

% M4 --> M2
M4M2 = M.M2.USD_std(end-id) - M.M4.USD_std(end-id);
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M2, dz, '-s', 'LineWidth', 2)
plot(M4M3, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
xline(0, 'k--')
xl = xlim;
xm = max(abs(xl));
xlim([-xm xm])
legend('M4 --> M1', 'M4 --> M2', 'M4 --> M3', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('$\Delta\sigma_T$', 'Interpreter', 'latex')
ylabel('$h$ [m]', 'Interpreter','latex')
title('semidiurnal')

%% MH

% M4 --> M1
M4M1 = M.M1.UMH_std(end-id) - M.M4.UMH_std(end-id);

% M4 --> M3
M4M3 = M.M3.UMH_std(end-id) - M.M4.UMH_std(end-id);

% M4 --> M2
M4M2 = M.M2.UMH_std(end-id) - M.M4.UMH_std(end-id);
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M2, dz, '-s', 'LineWidth', 2)
plot(M4M3, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
xline(0, 'k--')
xl = xlim;
xm = max(abs(xl));
xlim([-xm xm])
legend('M4 --> M1', 'M4 --> M2', 'M4 --> M3', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('$\Delta\sigma_T$', 'Interpreter', 'latex')
ylabel('$h$ [m]', 'Interpreter','latex')
title('Mid-High')

%% HF

% M4 --> M1
M4M1 = M.M1.UHF_std(end-id) - M.M4.UHF_std(end-id);

% M4 --> M3
M4M3 = M.M3.UHF_std(end-id) - M.M4.UHF_std(end-id);

% M4 --> M2
M4M2 = M.M2.UHF_std(end-id) - M.M4.UHF_std(end-id);
dz = M.M1.dz;

figure
plot(M4M1, dz, '-s', 'LineWidth', 2)
hold on
plot(M4M2, dz, '-s', 'LineWidth', 2)
plot(M4M3, dz, '-s', 'LineWidth', 2)
axis ij
axis square
grid minor
xline(0, 'k--')
xl = xlim;
xm = max(abs(xl));
xlim([-xm xm])
legend('M4 --> M1', 'M4 --> M2', 'M4 --> M3', 'Location', 'eastoutside')
set(gca, "FontSize", 18)
xlabel('$\Delta\sigma_T$', 'Interpreter', 'latex')
ylabel('$h$ [m]', 'Interpreter','latex')
title('High-Frequency')

%% Expansion Coefficents

