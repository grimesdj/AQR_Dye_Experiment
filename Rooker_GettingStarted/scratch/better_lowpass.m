% a better lowpass method

clear all
close all

%% Load Data

% grab files
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/*.mat');

% take a look inside (0_o)
for i = 1:length(files)
    fname = [files(i).folder,filesep,files(i).name];
    obj = matfile(fname);
    CF = obj.Config;
    SN = CF.SN;
    SN = matlab.lang.makeValidName(SN);
    tag = SN(1:3);

    % cant use data(i) because structs not exactly the same \fp
    data.(tag) = load(fname);

    bins = size(data.(tag).Velocity_East, 2);
    data.(tag).Config.binnum = 1:bins;
    
    % make smoothie
    fprintf('Smoothing %s\n', tag)
    [LPF(i).Velocity_East, LPF(i).fsd, LPF(i).idx] = hamming_filter(data.(tag).Velocity_East, 1/600, 1/data.(tag).Config.dt, 1);
    LPF(i).Velocity_North = hamming_filter(data.(tag).Velocity_North, 1/600, 1/data.(tag).Config.dt);


    fprintf('done!\n')

    %% East
    
    % O. G.
    figure
    ax1 = subplot(2, 1, 1);
    img = imagesc(data.(tag).Time, data.(tag).Config.binnum, data.(tag).Velocity_East');
    set(img, 'Alphadata', ~isnan(data.(tag).Velocity_East'))
    set(gca, 'YDir','normal')
    cb = colorbar;
    ylabel(cb, 'Velocity, [m/s]', 'FontSize', 18)
    colormap(cmocean('balance'))
    clim([-0.5 0.5])
    datetick('x', 'keeplimits')
    ylabel('Original')
    set(gca, 'FontSize', 18)
    
    % smoothed
    ax2 = subplot(2, 1, 2);
    img = imagesc(data.(tag).Time(LPF(i).idx), data.(tag).Config.binnum, LPF(i).Velocity_East');
    set(img, 'Alphadata', ~isnan(LPF(i).Velocity_East'))
    set(gca, 'YDir','normal')
    cb = colorbar;
    ylabel(cb, 'Velocity, [m/s]', 'FontSize', 18)
    colormap(cmocean('balance'))
    clim([-0.05 0.05])
    datetick('x', 'keeplimits')
    ylabel('Smoothed')
    set(gca, 'FontSize', 18)
    
    % eh
    linkaxes([ax1 ax2], 'x')
    sgtitle([tag ' East'], 'fontsize', 22)
    
    
    %% North
    
    % O.G.
    figure
    ax1 = subplot(2, 1, 1);
    img = imagesc(data.(tag).Time, data.(tag).Config.binnum, data.(tag).Velocity_North');
    set(img, 'Alphadata', ~isnan(data.(tag).Velocity_North'))
    set(gca, 'YDir','normal')
    cb = colorbar;
    ylabel(cb, 'Velocity, [m/s]', 'FontSize', 18)
    colormap(cmocean('balance'))
    clim([-0.5 0.5])
    datetick('x', 'keeplimits')
    ylabel('Original')
    set(gca, 'FontSize', 18)
    
    % smoothed
    ax2 = subplot(2, 1, 2);
    img = imagesc(data.(tag).Time(LPF(i).idx), data.(tag).Config.binnum, LPF(i).Velocity_North');
    set(img, 'Alphadata', ~isnan(LPF(i).Velocity_North'))
    set(gca, 'YDir','normal')
    cb = colorbar;
    ylabel(cb, 'Velocity, [m/s]', 'FontSize', 18)
    colormap(cmocean('balance'))
    clim([-0.05 0.05])
    datetick('x', 'keeplimits')
    ylabel('Smoothed')
    set(gca, 'FontSize', 18)
    
    % ehh
    linkaxes([ax1 ax2], 'x')
    sgtitle([tag, ' North'], 'fontsize', 22)

end


