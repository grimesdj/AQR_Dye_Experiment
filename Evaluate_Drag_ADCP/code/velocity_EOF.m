clear all
close all



% bandpass range IN HOURS (BP = 0 means bandpass off, BP = 1 means bandpass on)
BP = 0; 
lower_bound = 48;
upper_bound = 10;

% remove barotropic? (0 for no, 1 for yes)
rmBT = 1;

% initialize
FOV_fig(1) = figure("Position", [2250 150 1000 800]);
t1 = tiledlayout(2, 2);

FOV_fig(2) = figure("Position", [2250 150 1000 800]);
t2 = tiledlayout(2, 2);


EOF_fig(1) = figure("Position", [2250 150 1000 800]);
t3 = tiledlayout(2, 2);

EOF_fig(2) = figure("Position", [2250 150 1000 800]);
t4 = tiledlayout(2, 2);

Spectra_fig(1, 1) = figure;
t5 = tiledlayout(2, 2);

Spectra_fig(2, 1) = figure;
t6 = tiledlayout(2, 2);

Spectra_fig(1, 2) = figure;
t7 = tiledlayout(2, 2);

Spectra_fig(2, 2) = figure;
t8 = tiledlayout(2, 2);

BT_fig = figure;


moorings = {'M1', 'M2', 'M3'};
for mooring_ID = 1:length(moorings)
%% Load
mooring = moorings{mooring_ID};

% ADCP
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1', 'ADCP');
fname = mooring + "_10min_gridded_PCA.mat";
Moor = load(fullfile(fpath,fname));
dz = Moor.Zq(:, 1);

%% EOF
card = {'U', 'V'};
for d = 1:length(card)
    direction = card{d};
    field = direction + "_grid";

Vel_grid = fillmissing(Moor.(field), 'linear', 2, 'EndValues','nearest');
Y = Vel_grid';
if BP
    Y = bandpass(Y, [1/(lower_bound * 3600) 1/(upper_bound * 3600)], 1/600);
end

[L, EOFs, EC, Error, Skill,lam, Barotropic] = EOF(Y, [], rmBT);

% FOV BT vs BC

% Velocities
U  = Y;
UBT = Barotropic;
UBC = U - UBT;

% Variance
V = sum(var(U, 0, 1));
VBT = sum(var(UBT, 0, 1));
VBC = sum(var(UBC, 0, 1));

% FOV
FBT = VBT/V;
FBC = VBC/V;

% barotropic spectra
figure(BT_fig)
w = 720;
window = hamming(w);
noverlap = round(w*(2/3));
nfft = [];
fs = 1/600;
[pxx,f, pxxc] = pwelch(Barotropic,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
%subplot(2, 1, d)
loglog(f, pxx,'-', 'LineWidth', 2)
hold on
%loglog(f, pxxc, 'k--', 'LineWidth', 0.5)
grid on
title(sprintf('%s %s Barotropic Spectra', mooring, direction))
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')
ylim([0 10^2])


% FOV
FOV = L/sum(L);
figure(FOV_fig(d))
nexttile
plot(FOV, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','black')
ylabel('FOV')
xlabel('Mode #')
set(gca, 'FontSize', 18)
axis square
grid minor
ylim([0 1])
title(sprintf('%s', mooring))
sgtitle(sprintf('%s', direction))


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


%% plot modes

% consistent sign convention
for i = 1:size(EOFs, 2)
    if EOFs(end, i) > 0
        EOFs(:, i) = -1* EOFs(:,i);
        EC(:, i) = -1 * EC(:, i);
    end
end

figure(EOF_fig(d))
nexttile
plot(EOFs(:, 1), dz, 'k', 'LineWidth', 2)
hold on
plot(EOFs(:, 2), dz, 'r', 'LineWidth', 2)
%plot(EOFs(:, 3), dz, 'k--', 'LineWidth', 2)
axis square
ylabel('$h$ [m]', 'Interpreter','latex')
xlabel('$\mathrm{m}^2/\mathrm{s}^2$', 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
grid minor
axis ij
title(sprintf('%s', mooring))
sgtitle(sprintf('%s', direction))

%% save
savestr = mooring + "_EOF_depth_coords_" + direction + ".mat";
fpath = fullfile(fpath, '..', '..', 'L1', 'ADCP');
%save(fullfile(fpath, savestr), 'L', 'EOFs', 'EC', 'Error', 'Skill','lam','Barotropic', 'dz')


for modenum = 1:2
    figure(Spectra_fig(d, modenum))
    x = EC(:, modenum);
    w = 720;
    window = hamming(w);
    noverlap = round(w*(2/3));
    nfft = [];
    fs = 1/600;
    [pxx,f, pxxc] = pwelch(x,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
    
    nexttile
    loglog(f(3:end), pxx(3:end), 'k', 'LineWidth', 1, 'DisplayName', sprintf('Mode %d', modenum))
    hold on
    loglog(f(3:end), pxxc(3:end, :), 'm--', 'LineWidth', 0.25, 'DisplayName','95% Confidence')
    xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
    xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
    %xline(1/(21.2*3600), 'g--', 'Label', 'Interial Tide')
    title(sprintf('%s Velocity EC Spectrum %s', direction, mooring))
    grid minor
    set(gca, 'FontSize', 16)
end





sEOF1 = std(EC(:, 1)) * EOFs(:, 1);
sEOF2 = std(EC(:, 2)) * EOFs(:, 2);
figure, plot(Moor.mean_profile_V, dz, 'k', 'LineWidth', 1)
axis square
axis ij
hold on, plot(Moor.mean_profile_V + sEOF1, dz, 'LineWidth', 1)
hold on, plot(Moor.mean_profile_V - sEOF1, dz, 'LineWidth', 1)
hold on, plot(Moor.mean_profile_V + sEOF1 + sEOF2, dz, 'LineWidth', 1)
hold on, plot(Moor.mean_profile_V + sEOF1 - sEOF2, dz, 'LineWidth', 1)
hold on, plot(Moor.mean_profile_V - sEOF1 - sEOF2, dz, 'LineWidth', 1)
hold on, plot(Moor.mean_profile_V - sEOF1 + sEOF2, dz, 'LineWidth', 1)

%% Variance Bar Plot
FM1 = FBC * FOV(1);
Fnoise = FBC - FM1;

variance.(direction)(mooring_ID, :) = [FBT FM1 Fnoise] .* V;


end

figure(FOV_fig(1))
lgd = legend('Fraction of Variance', 'MC White Noise Floor');
lgd.Layout.Tile = 4;
lgd.NumColumns = 1;

figure(FOV_fig(2))
lgd = legend('Fraction of Variance', 'MC White Noise Floor');
lgd.Layout.Tile = 4;
lgd.NumColumns = 1;

figure(EOF_fig(1))
lgd = legend('1st Mode', '2nd Mode');%, '3rd Mode');
lgd.Layout.Tile = 4;
lgd.NumColumns = 1;

figure(EOF_fig(2))
lgd = legend('1st Mode', '2nd Mode');%, '3rd Mode');
lgd.Layout.Tile = 4;
lgd.NumColumns = 1;

figure(Spectra_fig(1, 1))
ax  = findall(gcf, 'Type', 'axes');
ylims = vertcat(ax.YLim); 
ymax = max(ylims(:,2));
set(ax, 'YLim', [0 ymax])
lgd = legend('Mode 1 Spectra', '95% Confidence');
lgd.Layout.Tile = 4;

figure(Spectra_fig(1, 2))
ax  = findall(gcf, 'Type', 'axes');
ylims = vertcat(ax.YLim); 
ymax = max(ylims(:,2));
set(ax, 'YLim', [0 ymax])
lgd = legend('Mode 2 Spectra', '95% Confidence');
lgd.Layout.Tile = 4;

figure(Spectra_fig(2, 1))
ax  = findall(gcf, 'Type', 'axes');
ylims = vertcat(ax.YLim); 
ymax = max(ylims(:,2));
set(ax, 'YLim', [0 ymax])
lgd = legend('Mode 1 Spectra', '95% Confidence');
lgd.Layout.Tile = 4;

figure(Spectra_fig(2, 2))
ax  = findall(gcf, 'Type', 'axes');
ylims = vertcat(ax.YLim); 
ymax = max(ylims(:,2));
set(ax, 'YLim', [0 ymax])
lgd = legend('Mode 2 Spectra', '95% Confidence');
lgd.Layout.Tile = 4;




%% save variance data
variance_U = variance.U;
variance_V = variance.V;

save('../../../../Kelp_data/data/2024_PROCESSED_DATA/Velocity_Variance', 'variance_U', 'variance_V');

end




%return
cmap = cmocean('haline');
cl = size(cmap, 1);
id = round(cl/4);
c = cmap(id, :);

figure(EOF_fig(1))
ax  = findall(gcf, 'Type', 'axes');

figure
m1 = EOFs(:, 1);
plot(m1, dz, 'r', 'LineWidth', 8)
Bt = mean(Barotropic, 1)*1e18;
hold on
plot(Bt, dz, 'k', 'LineWidth', 8)
Tp = m1 + Bt;
plot(Tp, dz, '-', 'LineWidth', 8, 'Color', c)
axis square
axis ij
xline(0, '--k', 'LineWidth', 3, 'label', '$u = 0$', 'Interpreter','latex', 'LabelOrientation','horizontal', 'FontSize', 16)
xticklabels([])
yticklabels([])
legend('Baroclinic Flow', 'Barotropic Flow', 'Actual Flow', 'Location', 'southeast', 'fontSize', 20)



fpath = fullfile('..', '..', '..', '..', 'Documents', 'YCSECA', '2026', 'figures');
print(gcf, fullfile(fpath, 'EOF_illustration.png'), '-dpng', '-r600')