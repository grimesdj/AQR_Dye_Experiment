%% EOF analysis of Temp data

clear all
close all

% initialize
FOV_fig = figure("Position", [2250 150 1000 800]);
t1 = tiledlayout(2, 2);

EOF_fig = figure("Position", [2250 150 1000 800]);
t2 = tiledlayout(2, 2);

Spectra_fig(1) = figure;
t3 = tiledlayout(2, 2);
Spectra_fig(2) = figure;
t4 = tiledlayout(2, 2);


%% Load Data
moorings = {'M1', 'M2', 'M3', 'M4'};
for mooring_ID = 1:4
mooring = moorings{mooring_ID};
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_10min_gridded.mat";
load(fullfile(fpath, savestr))




%% EOF
Y = Temp_grid';

[L, EOFs, EC, Error, Skill,lam] = EOF(Y);

% Ybar_fig(mooring_ID) = figure;
% 
% w = 1024;
% window = hamming(w);
% noverlap = w/2;
% nfft = [];
% fs = 1/600;
% 
% [pxx, f, pxxc] = pwelch(Ybar, window, noverlap, nfft, fs, 'ConfidenceLevel', 0.95);
% loglog(f, movmean(pxx, 3), 'LineWidth', 1)
% hold on
% loglog(f, movmean(pxxc, 3), 'k--', 'LineWidth', 0.5)
% grid on
% total variance
%sig = var(Y(:));

% variance explained
FOV = L/sum(L);
figure(FOV_fig)
nexttile
plot(FOV, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','black')
ylabel('FOV')
xlabel('Mode #')
set(gca, 'FontSize', 18)
axis square
grid minor
title(sprintf('%s', mooring), 'FontSize', 18)
ylim([0 1])

% monte-carlo noise floor
N = size(Y,1);
M = size(Y,2);
Ni = 1000;

F = zeros(Ni,M);

for ii = 1:Ni
    X = randn(N,M);

    [~,S,~] = svd(X,'econ');
    lam = diag(S).^2;

    F(ii,:) = lam'/sum(lam);
end

noise95 = prctile(F,95);
hold on
plot(noise95,'r--','LineWidth',2)
%text(mean(1:length(noise95)), mean(noise95) + 0.25*std(FOV), 'White Noise Floor', 'Color', 'r', 'EdgeColor', 'k')


% consistent sign convention
for i = 1:size(EOFs, 2)
    if EOFs(end, i) > 0
        EOFs(:, i) = -1* EOFs(:,i);
        EC(:, i) = -1 * EC(:, i);
    end
end

% first three EOFs
figure(EOF_fig)
nexttile
plot(EOFs(:, 1), dz, 'k', 'LineWidth', 2)
hold on
plot(EOFs(:, 2), dz, 'r--', 'LineWidth', 2)
%plot(EOFs(:, 3), dz, 'k--', 'LineWidth', 2)
axis ij
axis square
ylabel('Depth [m]')
xlabel('$^\circ\mathrm{C}^2$', 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
title(sprintf('%s', mooring), 'FontSize', 18)

% EC Spectra
for modenum = 1:2
    figure(Spectra_fig(modenum))
    x = EC(:, modenum);
    w = 720;
    window = hamming(w);
    noverlap = round(w*(2/3));
    nfft = [];
    fs = 1/600;
    [pxx,f, pxxc] = pwelch(x,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
    
    nexttile
    semilogx(f(3:end), pxx(3:end), 'k', 'LineWidth', 1, 'DisplayName', sprintf('Mode %d', modenum))
    hold on
    semilogx(f(3:end), pxxc(3:end, :), 'm--', 'LineWidth', 0.25, 'DisplayName','95% Confidence')
    xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
    xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
    %legend(sprintf('Mode %d', modenum), '95% Confidence')
    title(sprintf('%s EC Spectrum', mooring), 'FontSize', 18)
    set(gca, 'FontSize', 18)
    grid on
end

% store in struct
Data(mooring_ID).L = L;
Data(mooring_ID).EOFs = EOFs;
%Data(mooring_ID).sig = sig;
Data(mooring_ID).EC = EC;
Data(mooring_ID).Error = Error;
Data(mooring_ID).Skill = Skill;
Data(mooring_ID).lam = lam;
Data(mooring_ID).FOV = FOV;
Data(mooring_ID).Y = Y;
Data(mooring_ID).dz = dz;
Data(mooring_ID).Time = Time;

end
figure(FOV_fig)
lgd = legend('Fraction of Variance', 'MC White Noise Floor');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 2;

figure(EOF_fig)
lgd = legend('1st Mode', '2nd Mode');%, '3rd Mode');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 2;

figure(Spectra_fig(1))
lgd = legend('Mode 1 Spectra', '95% Confidence');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 2;

figure(Spectra_fig(2))
lgd = legend('Mode 2 Spectra', '95% Confidence');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 2;



for j = 1:2
%% Compare modes at all moorings
modes(j) = figure;
for i = 1:length(Data)
    EOFs    = Data(i).EOFs;
    dz      = Data(i).dz;
    plot(EOFs(:, j), dz, '-s', 'LineWidth', 2, 'DisplayName', moorings{i})
    hold on
end

axis ij
axis square
ylabel('Depth [m]')
xlabel('$^\circ\mathrm{C}^2$', 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
xlim([-0.5 0.5])
title(sprintf('Mode %d EOF', j))
legend
end

%% Figures
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'Summer2025','Rooker', 'figures');
print(FOV_fig, fullfile(fpath, 'FOV_fig.png'), '-dpng' ,'-r600')
print(EOF_fig, fullfile(fpath, 'EOF_fig.png'), '-dpng' ,'-r600')
print(Spectra_fig(1), fullfile(fpath, 'Spectra_fig_mode_1.png'), '-dpng' ,'-r600')
print(Spectra_fig(2), fullfile(fpath, 'Spectra_fig_mode_2.png'), '-dpng' ,'-r600')

print(modes(1), fullfile(fpath, 'modes_fig_1.png'), '-dpng', '-r600')
print(modes(2), fullfile(fpath, 'modes_fig_2.png'), '-dpng', '-r600')