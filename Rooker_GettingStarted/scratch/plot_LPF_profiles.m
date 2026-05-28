%% Plot Profiles

% Load M1
M1.ADCP = load('../../../../Kelp_data/data/Release2/L0/ADCP/M1_ADCP.mat');
M1.LPF = load('../../../../Kelp_data/data/Release2/L1/ADCP/M1_LPF_600sec.mat');
M1.Moor = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L1/mooring_M1.mat');

releasenum = 2;
release = string(releasenum);


depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00') ;
            datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];

M1.ADCP.Time = M1.ADCP.Time(M1.ADCP.Time >= depTime(1, releasenum) & M1.ADCP.Time <= depTime(2, releasenum));
M1.Moor.Time = M1.Moor.Time(M1.Moor.Time >= depTime(1, releasenum) & M1.Moor.Time <= depTime(2, releasenum));
Time = interp1(M1.ADCP.Moor, M1.Moor.Time);


% plot
figure
img = imagesc(datetime(M1.ADCP.Time, 'convertfrom', 'datenum'), M1.ADCP.bin_mab, M1.LPF.Velocity_North');
set(img, 'AlphaData', ~isnan(M1.LPF.Velocity_North'))
set(gca, 'YDir', 'normal')
cb = colorbar;
colormap(cmocean('balance'))
clim([-0.1 0.1])
ylabel('Meters Above Bottom')
ylabel(cb, 'Cross-shore Velocity, $u$ [m/s]', 'Interpreter','latex', 'FontSize', 18)
set([gca cb], 'fontsize', 18)
ylim([0 max(M1.ADCP.maxRNG)+std(M1.ADCP.maxRNG)])





