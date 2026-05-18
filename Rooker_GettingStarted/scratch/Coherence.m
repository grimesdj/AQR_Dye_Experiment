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
[ALCcxy,f] = mscohere(u_in,u_out,hamming(M),noverlap,NFFT,1/d);
ALC = figure;

% 95% confidence
alpha = 0.05;
T = 1./(f * 60 * 60);
K = floor((N - noverlap) / (M - noverlap));
nu = 1.5 * 2 * K; 
Ccrit = 1 - alpha^(1/(nu - 1));

% plot!
plot(T, ALCcxy, 'k', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 18)
set(gca, "FontSize", 20)
ylim([0 1])
xlim([0.15 25])
title('Alongshore Coherence', 'FontSize', 20)




%% Coherence between M1 and M2 North


[CSCcxy,f] = mscohere(v_in,v_out,hamming(M),noverlap,NFFT,1/d);


CSC = figure;
% plot!
plot(T, CSCcxy, 'k', 'LineWidth', 2.5)
set(gca, 'XDir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
xlabel('Period [hours]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 18)
set(gca, "FontSize", 20)
ylim([0 1])
xlim([0.15 25])
title('Cross-Shore Coherence ', 'FontSize', 20)



%% Compare Both
% %figure
% %plot(T, Ecxy, 'b', 'LineWidth', 2.5)
% %hold on
% plot(T, Ncxy, 'g', 'LineWidth', 2.5)
% set(gca, 'XDir', 'reverse')
% set(gca, 'XScale', 'log')
% set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
% xlabel('Period [hours]')
% ylabel('Coherence')
% yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
% set(gca, "FontSize", 18)
% ylim([0 1])
% xlim([0.15 25])
% legend('East', 'North', 'Location', 'best')
% title('M1 vs M2')

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



%% PSD for wave reduction

% East Outside
[NOpxx,f] = pwelch(v_out,hamming(M),noverlap,NFFT,1/d);

% East inside
[NIpxx,f] = pwelch(v_in,hamming(M),noverlap,NFFT,1/d);

% relative difference
df = NOpxx - NIpxx;
rdf = 1 - sqrt(NIpxx ./ NOpxx);
rdfm = movmean(rdf, length(NIpxx)/8);

T = (1./f)./60^2;

figure;
plot(T, NOpxx, 'b', 'LineWidth', 2)
hold on
plot(T, NIpxx, 'g', 'LineWidth', 2)
%plot(f, df, 'r--', 'LineWidth', 1.5)
set(gca, 'XScale', 'log')
set(gca, 'XDir', 'reverse')
set(gca, 'YScale', 'log')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
%xlim([1/(d*N) 1/(2*d)])
ylim([10^-3 10^2])
ylabel('Power Spectral Density, [m^2/s]')
xlabel('Period, T [hours]')
grid on
set(gca, 'FontSize', 18)
legend('Outside', 'Inside')


figure
plot(T, rdf*100, 'b', 'LineWidth', 1.5)
hold on
%plot(T, rdfm*100, 'r--', 'LineWidth', 2)
set(gca, 'XScale', 'log')
set(gca, 'XDir', 'reverse')
set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
ylabel('Relative Flow Reduction, [%]')
xlabel('Period, T [hours]')
grid on
set(gca, 'FontSize', 20)

threshold = 0.8 * max(CSCcxy); 
idxBand = CSCcxy >= threshold;
T_band = T(idxBand);
T_min = min(T_band)-0.5;
T_max = max(T_band);

xline([T_min T_max], 'm-.', 'LineWidth', 1.5)
text(sqrt(T_min*T_max), 0, {'High Cross-Shore','Coherence'}, 'Color','m', 'HorizontalAlignment','center', 'fontsize', 18);

lgd = legend('Relative Flow Reduction (RFR)', 'Location','southwest');



exportgraphics(gcf, '../../../../Documents/CSURF/2026/RelRed.png')

figure(ALC)
xline([T_min T_max], 'm-.', 'LineWidth', 1.5)
grid on

figure(CSC)
xline([T_min T_max], 'm-.', 'LineWidth', 1.5)
text(sqrt(T_min*T_max), 0.8, {'High Cross-Shore','Coherence'}, 'Color','m', 'HorizontalAlignment','center', 'fontsize', 18);
grid on

exportgraphics(ALC, '../../../../Documents/CSURF/2026/ECohere.png')
exportgraphics(CSC, '../../../../Documents/CSURF/2026/NCohere.png')





% 
% % compute cross-spectral density
% [Sxy,f] = cpsd(u_out, u_in, hamming(M), noverlap, NFFT, 1/d);
% 
% % compute PSD of driver
% [Sxx,~] = pwelch(u_out, hamming(M), noverlap, NFFT, 1/d);
% 
% % compute gain function
% G = Sxy ./ Sxx;
% 
% % predicted inside PSD from outside
% S_pred = abs(G).^2 .* Sxx;
% 
% % reduction PSD
% Reduction = Sxx - S_pred;
% RelRed = (Sxx - S_pred) ./ Sxx * 100;
% 
% figure
% plot(T, RelRed, 'k', 'LineWidth', 2)
% set(gca, 'XScale', 'log')
% set(gca, 'XDir', 'reverse')
% set(gca, 'XTick', [0.25 0.5 1 3 6 12 24])
% ylabel('Relative Flow Reduction, [%]')
% xlabel('Period, T [hours]')
% grid on
% set(gca, 'FontSize', 18)
% 

