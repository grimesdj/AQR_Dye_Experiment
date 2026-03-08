%% Coherence
% Analzying the Coherence of several parameters

clear all
close all

%% Load data
filestem = '../../../../Kelp_data/Summer2025/Rooker/Release2/LPF';

% Get all .mat files that aren't PCA
files_all = dir([filestem, '/*.mat']);
names = string({files_all.name});
keep_idx = ~contains(names, 'PCA');
files = files_all(keep_idx);

% downsample factor
d = 300;

% distrubute
AQD = load(fullfile(files(1).folder, files(1).name));
M1  = load(fullfile(files(2).folder, files(2).name));
M2  = load(fullfile(files(3).folder, files(3).name));
VEC = load(fullfile(files(4).folder, files(4).name));

% trim
AQDx = AQD.Velocity_X(1:d:length(M2.Velocity_X));
AQDy = AQD.Velocity_Y(1:d:length(M2.Velocity_X));
VECx = VEC.Velocity_X(1:d:length(M2.Velocity_X));
VECy = VEC.Velocity_Y(1:d:length(M2.Velocity_X));


%% Coherence between M1 and M2 North and East

% allocate vars
time    = M1.Time(1:d:end);
dtime   = datetime(time, 'ConvertFrom', 'datenum');
u_out   = M1.Velocity_X(1:d:end);
v_out   = M1.Velocity_Y(1:d:end);
u_in    = M2.Velocity_X(1:d:end);
v_in    = M2.Velocity_Y(1:d:end);

% assign cohere parameters
N = length(u_in);
M = floor(N/8);
noverlap = floor(M/2); % sticking with matlab defaults but explicitly defining
NFFT = M * 8; 

% run it!
[Ecxy,f] = mscohere(u_in,u_out,hamming(M),noverlap,NFFT,1/d);
figure

% 95% confidence
alpha = 0.05;
L = (2*N)/M;
Ccrit = 1 - alpha^(1/(L-1));
T = 1./(f * 60 * 60);



% plot!
plot(T, Ecxy, 'k', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
xlim([0.15 25])
title('East Coherence M1 vs M2', 'FontSize', 25)


%% Coherence between M1 and M2 North


[Ncxy,f] = mscohere(v_in,v_out,hamming(M),noverlap,NFFT,1/d);


figure
% plot!
plot(T, Ncxy, 'k', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
xlim([0.15 25])
title('North Coherence M1 vs M2', 'FontSize', 25)


%% Compare Both
figure
plot(T, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(T, Ncxy, 'g', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
xlim([0.15 25])
legend('East', 'North', 'Location', 'best')
title('M1 vs M2')

%% Check inside instruments

% M2 and VEC
[Ecxy,f] = mscohere(VECx,u_in,hamming(M),noverlap,NFFT,1/d);
[Ncxy,f] = mscohere(VECy,v_in,hamming(M),noverlap,NFFT,1/d);

figure
tiledlayout(2, 2)
nexttile
plot(T, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(T, Ncxy, 'g', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
xlim([0.15 25])
title('M2 vs VEC')

nexttile

% M2 and AQD
[Ecxy,f] = mscohere(AQDx,u_in,hamming(M),noverlap,NFFT,1/d);
[Ncxy,f] = mscohere(AQDy,v_in,hamming(M),noverlap,NFFT,1/d);

plot(T, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(T, Ncxy, 'g', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
xlim([0.15 25])
title('M2 vs AQD')

nexttile

% AQD and VEC
[Ecxy,f] = mscohere(AQDx,VECx,hamming(M),noverlap,NFFT,1/d);
[Ncxy,f] = mscohere(AQDy,VECy,hamming(M),noverlap,NFFT,1/d);

plot(T, Ecxy, 'b', 'LineWidth', 2.5);
hold on
plot(T, Ncxy, 'g', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
xlim([0.15 25])
title('VEC vs AQD')

lgd = legend(gca, 'East', 'North');
lgd.Layout.Tile = 4;














% Transfer function
Mt = floor(N/8);
[H, f] = tfestimate(u_in, u_out, hamming(M), noverlap, NFFT, 1/d);

% --- Disc start ---

coh_idx = Ncxy > Ccrit;

% FFT of velocity
U = fft(u_in);
N = length(U);


% Only multiply in coherent frequency bins
H_full = zeros(size(U));  
H_full(coh_idx) = H(coh_idx);

% Apply transfer function
OUT_coh_fft = H_full .* U;

% Inverse FFT to get coherent temperature time series
OUT_coherent = ifft(OUT_coh_fft, 'symmetric');


t = (0:length(u_in)-1)/(1/d);

figure;
plot(time, u_out, 'k', 'DisplayName','Raw OUT'); hold on;
plot(time, OUT_coherent, 'r', 'DisplayName','Coherent Outside North Velocities',  'LineWidth', 1.5);
legend;
xlabel('Time [s]');
ylabel('vel');
title('Coherent Outside Velocity Driven by Inside Velocity');
datetick(gca,'keeplimits')

% --- Disc end ---
