%% Parse Temps

clear all
close all

load("../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat")


% setup for freq transform
N = length(R23.Time);
timerange = seconds(days(R1.Time(end) - R1.Time(1)));
dt = round(N/timerange);
fs = 1/dt;

time = datetime(R23.Time, 'ConvertFrom', 'datenum');

% original
figure
plot(time, R23.Temperature, 'LineWidth', 1)
title('original')

% de-trend
temp = detrend(R23.Temperature);
figure
%temp = temp - movmean(temp, fs*2000);
plot(time, temp, 'LineWidth', 1)
title('detrend')

% Freq
w = floor(N/2);
noverlap = w/2;
NFFT = floor(w*2);

[PSD, f] = pwelch(temp, w, noverlap ,NFFT, fs);
figure
loglog(f, PSD, 'LineWidth',2)
legend(sprintf('%.3f MAB', R1.mab(1)), ...
    sprintf('%.3f MAB', R23.mab(2)), ...
    sprintf('%.3f MAB', R23.mab(3)), ...
    sprintf('%.3f MAB', R23.mab(4)), ...
    sprintf('%.3f MAB', R23.mab(5)))

title('PSD spec')

%% low-pass?
data.T = temp(:, 2);
data.Time = time;
data.dt = dt;


[T_all, Temp_all] = LPF({data});

figure
plot(T_all{1}, Temp_all{1})

%save() % save PSD in its own folder in data and use coherence.m to analyze

function [t_all, T_all] = LPF(dataCell, labels, colors)
%   
% USAGE: [T_all, U_all, V_all] = LPF(dataCell, labels, colors)
% 
%   takes time series data and 
%  returns 3 cell arrays, containing 10-min averaged data, u, v, and time
% 
%   dataCell = cell array (# of files x 1) containing data structures
%   labels = {optional} Cell array of labels for graphs
%   colors = {optional} Cell array of RGB vals for graphs
%    
%
%   T_all = interpolated time
%   U_all = East Velocities
%   V_all = North Velocities
%   

%% Defaults
if nargin < 3
    colors = num2cell(lines(length(dataCell)), 2);
    if nargin < 2
        labels = num2cell(1:length(dataCell));
    end
end



for i = 1:length(dataCell)

    data = dataCell{i};

   
    % Extract and clean data
    T = data.T;
    t = data.Time;
    dt = data.dt;
    
    % Lowpass filter
    Nf = 600 / dt;
    flt = hamming(Nf); flt = flt / sum(flt);
    nanFlag = ~isnan(T);
    flag_flt = conv(nanFlag, flt, 'same');
    T(~nanFlag) = 0; v(~nanFlag) = 0;
    T_filt = conv(T, flt, 'same') ./ flag_flt;
    % Plot time series
    figure;
    plot(t, T_filt,'LineWidth', 1.5);
  
    % % Downsample to 10-minute intervals
    % [~, unique_idx] = unique(minutes(t - t(1)));
    % t = t(unique_idx);
    % T_filt = T_filt(unique_idx);
    % tq = t(1):minutes(10):t(end);
    % T_ds = interp1(t, T_filt, tq, 'linear');
    
    % Save Global Vars
    T_all{i} = T_filt;
    
    t_all{i} = t;
    
end
xlabel('Time', 'FontSize', 18)
sgtitle('Temperature at 1.2 MAB', 'FontSize', 25)


% Add legends
%legend(gca, h1, labels, 'Location', 'best')

datetick(gca);

% Set colors
for i = 1:length(dataCell)
    h1(i).Color = colors{i};

end
yline(0); 

end

% EOF