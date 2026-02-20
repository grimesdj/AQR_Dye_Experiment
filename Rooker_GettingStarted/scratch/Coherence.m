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

% distrubute
AQD = load(fullfile(files(1).folder, files(1).name));
M1  = load(fullfile(files(2).folder, files(2).name));
M2  = load(fullfile(files(3).folder, files(3).name));
VEC = load(fullfile(files(4).folder, files(4).name));

% trim
AQDx = AQD.Velocity_X(1:length(M2.Velocity_X));
AQDy = AQD.Velocity_Y(1:length(M2.Velocity_X));
VECx = VEC.Velocity_X(1:length(M2.Velocity_X));
VECy = VEC.Velocity_Y(1:length(M2.Velocity_X));


%% Coherence between M1 and M2 North and East

% allocate vars
time    = M1.Time;
dtime   = datetime(time, 'ConvertFrom', 'datenum');
u_out   = M1.Velocity_X;
v_out   = M1.Velocity_Y;
u_in    = M2.Velocity_X;
v_in    = M2.Velocity_Y;

% assign cohere parameters
N = length(u_in);
M = floor(N/8);
noverlap = floor(M/2); % sticking with matlab defaults but explicitly defining
NFFT = M*8;

% run it!
[Ecxy,f] = mscohere(u_in,u_out,hamming(M),noverlap,NFFT,1/600);
figure

% 95% confidence
alpha = 0.05;
L = (2*N)/M;
Ccrit = 1 - alpha^(1/(L-1));
%T = 1./(f * 60 * 60);

% plot!
plot(f, Ecxy, 'k', 'LineWidth', 2.5)
%set(gca, 'XDir', 'reverse')
% set(gca, 'XScale', 'log')
% ax = gca;
% ax.XAxis.Exponent = 0;
% ax.XAxis.TickLabelFormat = '%f';
xlabel('Frequency [Hz]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
title('East Coherence M1 vs M2', 'FontSize', 25)


%% Coherence between M1 and M2 North


[Ncxy,f] = mscohere(v_in,v_out,hamming(M),noverlap,NFFT,1/600);

figure
% plot!
plot(f, Ncxy, 'k', 'LineWidth', 2.5)
xlabel('Frequency [Hz]')
ylabel('Coherence')
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
set(gca, "FontSize", 18)
ylim([0 1])
title('North Coherence M1 vs M2', 'FontSize', 25)


%% Compare Both
figure
plot(f, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(f, Ncxy, 'g', 'LineWidth', 2.5)
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
legend('East', 'North')
set(gca, "FontSize", 18)
ylim([0 1])
title('M1 vs M2')



%% Check inside instruments

% M2 and VEC
[Ecxy,f] = mscohere(VECx,M2.Velocity_X,hamming(M),noverlap,NFFT,1/600);
[Ncxy,f] = mscohere(VECy,M2.Velocity_Y,hamming(M),noverlap,NFFT,1/600);

figure
plot(f, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(f, Ncxy, 'g', 'LineWidth', 2.5)
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
legend('East', 'North')
set(gca, "FontSize", 18)
ylim([0 1])
title('M2 vs VEC')

% M2 and AQD
[Ecxy,f] = mscohere(AQDx,M2.Velocity_X,hamming(M),noverlap,NFFT,1/600);
[Ncxy,f] = mscohere(AQDy,M2.Velocity_Y,hamming(M),noverlap,NFFT,1/600);

figure
plot(f, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(f, Ncxy, 'g', 'LineWidth', 2.5)
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
legend('East', 'North')
set(gca, "FontSize", 18)
ylim([0 1])
title('M2 vs AQD')

% M2 and AQD
[Ecxy,f] = mscohere(AQDx,VECx,hamming(M),noverlap,NFFT,1/600);
[Ncxy,f] = mscohere(AQDy,VECy,hamming(M),noverlap,NFFT,1/600);

figure
plot(f, Ecxy, 'b', 'LineWidth', 2.5)
hold on
plot(f, Ncxy, 'g', 'LineWidth', 2.5)
yline(Ccrit, 'r--', 'Label', '95% Confidence', 'LineWidth', 1.5)
legend('East', 'North')
set(gca, "FontSize", 18)
ylim([0 1])
title('VEC vs AQD')