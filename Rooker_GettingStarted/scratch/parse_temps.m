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
    sprintf('%.3f MAB', R1.mab(2)), ...
    sprintf('%.3f MAB', R1.mab(3)), ...
    sprintf('%.3f MAB', R1.mab(4)), ...
    sprintf('%.3f MAB', R1.mab(5)))

title('PSD spec')

%% low-pass?
y = lowpass(temp, 1/600, fs);

figure
plot(y)
