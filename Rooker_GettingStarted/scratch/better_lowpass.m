% a better lowpass method

clear all
close all

%% Load Data

% grab files
load('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/KELP2_AquadoppHR_L0.mat');

LPF = hamming_filter(Velocity_North(:,10), 1/600, Config.dt);








function y = hamming_filter(x, wpass, fs)


N = (1/wpass) * fs;

w = hamming(N);
w = w/sum(w);

mask = ~isnan(x);

x0 = x;
x0(~mask) = 0;

y = conv(x0,w,'same') ./ conv(double(mask),w,'same');

y(conv(double(mask),w,'same') < 0.5) = NaN;

figure
subplot(2, 1, 1)
plot(x, 'LineWidth', 2)
ylabel({'Original', 'Data'})
set(gca, 'FontSize', 18)

subplot(2, 1, 2)
plot(y, 'LineWidth', 2)
ylabel({'Original', 'filtered'})
set(gca, 'FontSize', 18)


figure
scatter(x,y)

end


