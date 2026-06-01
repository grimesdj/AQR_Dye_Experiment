%% Produce Low-Pass data

clear all
close all

%% User Input Data

releasenum = 2; % Enter Release number here
releasenum = string(releasenum);

fc = 600; % Enter window length ( ~min period pass ) [ seconds ]


%% Load Data

% grab files
filepath = ['../../../../Kelp_data/data/Release' + releasenum + '/L0/ADCP/*.mat'];
files = dir(filepath);

% take a look inside (0_o)
for i = 1:length(files)
    fname = [files(i).folder,filesep,files(i).name];
    obj = matfile(fname);
    CF = obj.Config;
    SN = files(i).name(1:end-4);
    SN = matlab.lang.makeValidName(SN);
    tag = SN(1:3);
    SN = strip(SN, '_');
    tag = strip(tag, '_');
    fprintf('loading %s...\n', tag)

    % cant use data(i) because structs not exactly the same \fp
    data.(tag) = load(fname);

    bins = size(data.(tag).Velocity_East, 2);
    data.(tag).Config.binnum = 1:bins;
    
    %correct time
    if strcmp(tag, 'AQD')
        thresh = 100; % keeping the AQD nans just because theres so many data gaps and the correlations are good
    else
        thresh = 80;
    end

    % not sure how the old configs for the M1/M2 came to be
    if ~isfield(data.(tag).Config, "dt")
    data.(tag).Config.dt = median(diff(data.(tag).Time)) * 86400;
    end

    [~, data.(tag).Velocity_East]  = correct_burst(data.(tag).Time, data.(tag).Velocity_East, 1/data.(tag).Config.dt);
    [data.(tag).Time, data.(tag).Velocity_North] = correct_burst(data.(tag).Time, data.(tag).Velocity_North, 1/data.(tag).Config.dt);

    % make smoothie
    fprintf('Smoothing %s\n', tag)
    [LPF(i).Velocity_East, LPF(i).fsd, LPF(i).idx] = hamming_filter(data.(tag).Velocity_East, 1/fc, 1/data.(tag).Config.dt, 1, 1, thresh);
    LPF(i).Velocity_North = hamming_filter(data.(tag).Velocity_North, 1/fc, 1/data.(tag).Config.dt, 1, 1, thresh);
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
    drawnow

    % save or something
    savename = string(sprintf('%s_LPF_%dsec.mat', tag, fc));
    savepath = fullfile(files(i).folder, '..', '..', 'L1', 'ADCP', savename);
    fprintf('saving %s\n', savename)
    S = LPF(i);
    save(savepath, '-struct', "S")


end



%% Functions

function [y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt)
% 
% hamming_filter
%     applies a low-pass filter to data using a hamming filter
% 
% USAGE: [y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds, nan_filt)
%        [y, fsd, idx] = hamming_filter(x, wpass, fs);
% 
%     INPUTS:
%             x         = vector data to be low-passed
%             wpass     = frequency cutoff (Hz)
%             fs        = original sampling rate (Hz) (defaults to 1 Hz)
%             fig       = (optional) figure flag: 1 for figures, 0 for none (defaults to 0)
%             ds        = (optional) downsample flag: 1 for downsampling, 0 for none (defaults to 1)
%             nan_filts = (optional) NaN any windows that are over a certain
%                          percent of nans
% 
%     OUTPUTS:
%             y     = filtered data vector
%             fsd   = new sample rate (Hz)
%             idx     = Downsample indicies
%             
% 
%             written by Jason Rooker, May 2026


                
% Defaults
if nargin < 6 || isempty(nan_filt)
    nan_filt = 0;
    if nargin < 5 || isempty(ds)
        ds = 1;
    
        if nargin < 4 || isempty(fig)
            fig = 0;
    
            if nargin < 3 || isempty(fs)
                fs = 1;
            end
        end
    end
end


% window length
N = round((1/wpass) * fs);

w = hamming(N);
w = w/sum(w);

mask = ~isnan(x);

x0 = x;
x0(~mask) = 0;

% normalization factor
norm = conv2(double(mask),w,'same');
y = conv2(x0,w,'same') ./ norm;

% NaN-dense filter to protect statistical integrety for burst data
thresh = 1-nan_filt/100;
fprintf('Any points constucted from more than %d%% NaNs set equal to NaN\n', nan_filt)
y(norm < thresh) = NaN;

% downsample
if ds
    D = max(1, round(fs / (4 * wpass)));
else
    D = 1;
end

idx = 1:D:length(x);
y = y(idx, :);
fsd = fs/D;

% (optional) figs to make sure it worked
if fig

    xfig = 1:length(x);
    figure
    ax1 = subplot(2, 1, 1);
    plot(xfig, x, '.', 'LineWidth', 2)
    ylabel({'Original', 'Data'})
    set(gca, 'FontSize', 18)
    grid minor
    
    ax2 = subplot(2, 1, 2);
    plot(xfig(1:D:end), y, '.', 'LineWidth', 2)
    ylabel({'Filtered', 'Data'})
    set(gca, 'FontSize', 18)
    ylim(ax1.YLim)
    linkaxes([ax1 ax2], 'x')
    grid minor

end

end
% EOF

function [t_fix, x_fix] = correct_burst(t, x, fs)
 % 
 % correct_burst - if data has bursts of data or time gaps, the missing timesteps can
 %                     distort figures and filters! this function fills gaps
 % 
 %    USAGE:  
 %            [t_fix, x_fix] = correct_burst(t, x, fs)
 %            [t_fix, x_fix] = correct_burst(t, x)
 % 
 %    INPUTS:
 %            t = time vector (datenum or datetime)
 %            x = data array
 %            fs = (optional) intended sample frequency [Hz]
 % 
 %    OUTPUTS:
 %            t_fix = time with gaps interpolated (NOTE: length(t_fix) > length(t))
 %            x_fix = data array with NaN columns in time gaps
 % 
 %            written by Jason Rooker May 2026


% force datnum
if strcmp(class(t), datetime)
    t = datenum(t);
end

% default with optional fs
if nargin < 3 || isempty(fs)
    dtds = mode(diff(t))/86400;
else
    dtds = 1/(fs * 86400);
end

% collect segments
[startINDs, endINDs] = Segment(t, 2);

% make a time vector
t_fix = [];
x_fix = [];

% Correct
for i = 1:length(startINDs)-1
    t_fix = [t_fix; t(startINDs(i):endINDs(i)); (t(endINDs(i)):dtds:t(startINDs(i+1)))'];
    x_fix = [x_fix; x(startINDs(i):endINDs(i), :); nan(length(t(endINDs(i)):dtds:t(startINDs(i+1))), size(x, 2))];
end
t_fix = [t_fix; t(startINDs(end):endINDs(end))];
x_fix = [x_fix; x(startINDs(end):endINDs(end), :)];


if i
    if i == 1
        fprintf("%d time gap corrected\n", i)
    else
        fprintf("%d time gaps corrected\n", i)
    end
end

% plot fixed data
figure
subplot(2, 1, 1);
plot(x, '.')
ylabel('Original raw data')

subplot(2, 1, 2);
plot(x_fix, '.')
ylabel({'Time gaps', 'filled by NaNs'})
xlabel('indcies (they should be different)')

% plot fixed time
figure
plot(t, '.')
hold on
plot(t_fix, '.')
legend('Original Time', 'Corrected Time')






%% Outside functions

function [startINDs, endINDs] = Segment(time,cuttoff);
% 
% Usage: [startINDs, endINDs] = Segment(time,cuttoff);
%
% This function will log the start and end indices for nearly
% monotonic sections amidst large steps. time is the vector of time
% stamps to be split up and cuttof is the max number of timesteps to
% allow before generating new segment. 
%
% written by Dr. Derek Grimes

%Clip the timeseries into monotonic segments
[r,c] = size(time);

dt = diff(time);% in seconds

nominal_sample_rate = nanmedian(dt);

seg_ind = find(dt >= ones(size(dt))*cuttoff*nominal_sample_rate);

startINDs = [0; seg_ind]+1;
endINDs = [seg_ind; r];
end
end

% EOF