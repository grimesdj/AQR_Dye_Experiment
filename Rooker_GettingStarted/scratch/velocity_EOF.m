%% Velocity_EOF

clear all
close all

%% User Input data

% bandpass range IN HOURS (BP = 0 means bandpass off, BP = 1 means bandpass on)
BP = 1; 
lower_bound = 48;
upper_bound = 8;

% remove barotropic? (0 for no, 1 for yes)
rmBT = 1;

% initialize
FOV_fig = figure("Position", [2250 150 1000 800]);
t1 = tiledlayout(2, 2);

EOF_fig = figure("Position", [2250 150 1000 800]);
t2 = tiledlayout(2, 2);


moorings = {'M1', 'M2', 'M3'};
for mooring_ID = 1:length(moorings)
%% Load
mooring = moorings{mooring_ID};

% ADCP
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L0', 'ADCP');
fname = "ADCP_" + mooring + "_L0_10min.mat";
M = load(fullfile(fpath,fname));
%M.Config = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_config.mat');


% apply qcFlag
fprintf('applying qcFlag...\n')
M.qcFlag(isnan(M.qcFlag)) = 0;
qcFlag = round(M.qcFlag); % omg the qcFlag got smoothed
fields = fieldnames(M);
for i = 1:length(fields)
    field = fields{i};
    dum = M.(field);
    if size(dum) == size(qcFlag) & ~strcmp(field, 'qcFlag')
        dum(~qcFlag) = NaN;
        M.(field) = dum;
    end
end

%% Grid data

% original reference
hbar = nanmean(M.Pressure); % mean depth
bin_dep = hbar - M.bin_mab; % each bin at mean depth
bin_dep = [hbar; bin_dep]; % add bottom 0 bin
o = M.Velocity_East(1, :);%zeros(1, length(M.Time));
U = [o;M.Velocity_East]; % add 0 velocity at seafloor

% reference grid
[t_grid, z_grid] = meshgrid(M.Time, bin_dep);

% clear nans
mask = ~isnan(U);
U = U(mask);
t_grid = t_grid(mask);
z_grid = z_grid(mask);

% interpolate to a surface
F = scatteredInterpolant(t_grid(:), z_grid(:), U(:), 'linear', 'nearest');

% query grid
z_max = max(bin_dep);
z_min = hbar - min(M.Pressure);
step = 0.5;
dz = z_min+2:step:z_max;
[Tq, Zq] = meshgrid(M.Time, dz);

% Eval at query grid
Vel_grid = F(Tq, Zq);


% mean prfiles
mean_profile = nanmean(Vel_grid, 2);
std_profile = nanstd(Vel_grid,[], 2);

figure
subplot(1, 2, 1);
plot(mean_profile, dz,'-s', 'LineWidth', 2)
xlabel('$u$ [m/s]', 'Interpreter','latex')
ylabel('$h$ [m]', 'Interpreter','latex')
set(gca, 'FontSize', 18)
title('Mean Profile')
grid minor
axis ij

subplot(1, 2, 2);
plot(std_profile, dz, '-s', 'LineWidth', 2)
xlabel('$u$ [m/s]', 'Interpreter','latex')
ylabel('$h$ [m]', 'Interpreter','latex')
set(gca, 'FontSize', 18)
title('$\sigma$ Profile', 'Interpreter','latex')
grid minor
axis ij
%% EOF

Vel_grid = fillmissing(Vel_grid, 'linear', 2, 'EndValues','nearest');
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
figure
w = 720;
window = hamming(w);
noverlap = round(w*(2/3));
nfft = [];
fs = 1/600;
[pxx,f, pxxc] = pwelch(Barotropic,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
semilogx(f, pxx, 'LineWidth', 1)
hold on
semilogx(f, pxxc, 'k--', 'LineWidth', 0.5)
grid on
title(sprintf('%s Barotropic Spectra', mooring))
xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
xline(1/(21.2 * 3600), 'g--', 'Label', 'intertial tide')



% FOV
FOV = L/sum(L);
figure(FOV_fig)
nexttile
plot(FOV, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','black')
ylabel('FOV')
xlabel('Mode #')
set(gca, 'FontSize', 18)
axis square
grid minor
ylim([0 1])
title(sprintf('%s', mooring))


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

figure(EOF_fig)
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

%% save
savestr = mooring + "_EOF_depth_coords.mat";
fpath = fullfile(fpath, '..', '..', 'L1', 'ADCP');
save(fullfile(fpath, savestr), 'L', 'EOFs', 'EC', 'Error', 'Skill','lam','Barotropic', 'dz')
save(fullfile(fpath, mooring + "_10min_gridded_East.mat"), 'Tq', 'Zq', 'Vel_grid')


figure
for modenum = 1:2
    x = EC(:, modenum);
    w = 720;
    window = hamming(w);
    noverlap = round(w*(2/3));
    nfft = [];
    fs = 1/600;
    [pxx,f, pxxc] = pwelch(x,window,noverlap,nfft,fs, 'ConfidenceLevel', 0.95);
    
    subplot(2, 1, modenum)
    semilogx(f(3:end), pxx(3:end), 'k', 'LineWidth', 1, 'DisplayName', sprintf('Mode %d', modenum))
    hold on
    semilogx(f(3:end), pxxc(3:end, :), 'm--', 'LineWidth', 0.25, 'DisplayName','95% Confidence')
    xline(1/86400, 'b--', 'label', 'Diurnal', 'LineWidth', 1)
    xline(2/86400, 'b--', 'label', 'Semi-Diurnal', 'LineWidth', 1)
    xline(1/(21.2*3600), 'g--', 'Label', 'Interial Tide')
    title('Velocity EC Spectrum')
    grid minor
    legend
end



sEOF1 = std(EC(:, 1)) * EOFs(:, 1);
sEOF2 = std(EC(:, 2)) * EOFs(:, 2);
figure, plot(mean_profile, dz, 'k', 'LineWidth', 1)
axis square
axis ij
hold on, plot(mean_profile + sEOF1, dz, 'LineWidth', 1)
hold on, plot(mean_profile - sEOF1, dz, 'LineWidth', 1)
hold on, plot(mean_profile + sEOF1 + sEOF2, dz, 'LineWidth', 1)
hold on, plot(mean_profile + sEOF1 - sEOF2, dz, 'LineWidth', 1)
hold on, plot(mean_profile - sEOF1 - sEOF2, dz, 'LineWidth', 1)
hold on, plot(mean_profile - sEOF1 + sEOF2, dz, 'LineWidth', 1)

%% Variance Bar Plot
FM1 = FBC * FOV(1);
Fnoise = FBC - FM1;

variance(mooring_ID, :) = [FBT FM1 Fnoise] .* V;


end

figure(FOV_fig)
lgd = legend('Fraction of Variance', 'MC White Noise Floor');
lgd.Layout.Tile = 4;
lgd.NumColumns = 1;

figure(EOF_fig)
lgd = legend('1st Mode', '2nd Mode');%, '3rd Mode');
lgd.Layout.Tile = 4;
lgd.NumColumns = 1;




%% save variance data
save('../../../../Kelp_data/data/2024_PROCESSED_DATA/Velocity_Variance_East', 'variance');





