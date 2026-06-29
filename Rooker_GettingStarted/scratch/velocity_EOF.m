%% Velocity_EOF

clear all
close all

%% Load
mooring_ID = 1;
moorings = {'M1', 'M2', 'M3'}; % M3 doesnt have ADCP data yet
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
[t_grid, z_grid] = meshgrid(M.Time, M.bin_mab);

h = M.Pressure;
h = fillmissing(h, 'spline');
sig = z_grid./h;

dz = max(sig(1, :)):0.1:0.8;

% interpolate
F = scatteredInterpolant(t_grid(:), sig(:), M.Velocity_North(:), 'linear', 'nearest');
[Tq, Zq] = meshgrid(M.Time, dz);
Vel_grid = F(Tq, Zq);

% plot mean Profiles

mean_profile = nanmean(Vel_grid, 2);
std_profile = nanstd(Vel_grid,[], 2);

figure
subplot(1, 2, 1);
plot(mean_profile, dz,'-s', 'LineWidth', 2)
xlabel('$u$ [m/s]', 'Interpreter','latex')
ylabel('$z/h$', 'Interpreter','latex')
set(gca, 'FontSize', 18)
title('Mean Profile')
grid minor

subplot(1, 2, 2);
plot(std_profile, dz, '-s', 'LineWidth', 2)
xlabel('$u$ [m/s]', 'Interpreter','latex')
ylabel('$z/h$', 'Interpreter','latex')
set(gca, 'FontSize', 18)
title('Std Profile')
grid minor
%% EOF

Vel_grid = fillmissing(Vel_grid, 'linear', 2, 'EndValues','nearest');
Y = Vel_grid';
[L, EOFs, EC, Error, Skill,lam, Barotropic] = EOF(Y, [], 1);

% FOV
FOV = L/sum(L);
figure
plot(FOV, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','black')
ylabel('FOV')
xlabel('Mode #')
set(gca, 'FontSize', 18)
axis square
grid minor
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
legend('FOV', 'MC Noise Floor')


%% plot modes

% consistent sign convention
for i = 1:size(EOFs, 2)
    if EOFs(end, i) > 0
        EOFs(:, i) = -1* EOFs(:,i);
        EC(:, i) = -1 * EC(:, i);
    end
end

figure
plot(EOFs(:, 1), dz, 'k', 'LineWidth', 2)
hold on
plot(EOFs(:, 2), dz, 'r', 'LineWidth', 2)
%plot(EOFs(:, 3), dz, 'k--', 'LineWidth', 2)
axis square
ylabel('$z/h$', 'Interpreter','latex')
xlabel('$\mathrm{m}^2/\mathrm{s}^2$', 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
legend('1st mode', '2nd mode', 'Location','eastoutside')
grid minor

%% save
savestr = mooring + "_EOF.mat";
fpath = fullfile(fpath, '..', '..', 'L1', 'ADCP');
save(fullfile(fpath, savestr), 'L', 'EOFs', 'EC', 'Error', 'Skill','lam','Barotropic')



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
    title('Velocity EC Spectrum')
    grid minor
    legend
end