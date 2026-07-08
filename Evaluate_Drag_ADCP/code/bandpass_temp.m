%% bandpass Temp

clear all
close all

%% Load Data
mooring = 'M1';
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_EOF.mat";
M1 = load(fullfile(fpath, savestr));

mooring = 'M2';
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_EOF.mat";
M2 = load(fullfile(fpath, savestr));

mooring = 'M3';
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_EOF.mat";
M3 = load(fullfile(fpath, savestr));

mooring = 'M4';
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_EOF.mat";
M4 = load(fullfile(fpath, savestr));

%% Recontruct Mode 1

% M1
phi = M1.EC(:, 1);
A   = M1.EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);

M1.S = S;
M1.sig = sig;

% M2
phi = M2.EC(:, 1);
A   = M2.EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);

M2.S = S;
M2.sig = sig;

% M3
phi = M3.EC(:, 1);
A   = M3.EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);

M3.S = S;
M3.sig = sig;

% M4

phi = M4.EC(:, 1);
A   = M4.EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);

M4.S = S;
M4.sig = sig;



%% Bandpass
M(1).S = M1.S;
M(2).S = M2.S;
M(3).S = M3.S;
M(4).S = M4.S;

for i = 1:length(M)
    U = M(i).S;
    figure
    % Subtidal
    x = U;
    fpass = [1/(length(x) * 600) 1/(32 * 3600)];
    fs = 1/600;
    M(i).UST = bandpass(x,fpass,fs);
    [pxx, f] = pwelch(UST(:, 1));
    loglog(f, pxx)
    hold on
    
    % diurnal
    x = U - UST;
    fpass = [1/(32 * 3600) 1/(16 * 3600)];
    M(i).UDU = bandpass(x,fpass,fs);
    [pxx, f] = pwelch(UDU(:, 1));
    loglog(f, pxx)
    
    % semi-diurnal
    x = U - UST - UDU;
    fpass = [1/(16 * 3600) 1/(8 * 3600)];
    M(i).USD = bandpass(x,fpass,fs);
    [pxx, f] = pwelch(USD(:, 1));
    loglog(f, pxx)
    
    % high freq
    M(i).UHF = U - UST - UDU - USD;
    [pxx, f] = pwelch(UHF(:, 1));
    loglog(f, pxx)
    
end

% DU std
M1.UDU_std = std(M(1).UDU, [], 1);
M2.UDU_std = std(M(2).UDU, [], 1);
M3.UDU_std = std(M(3).UDU, [], 1);
M4.UDU_std = std(M(4).UDU, [], 1);

figure
% compare std profiles................