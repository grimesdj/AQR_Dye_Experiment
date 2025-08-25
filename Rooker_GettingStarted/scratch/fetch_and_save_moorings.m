% fetch and save moorings
clear all
close all


%% M1 R1
disp('M1 R1')
Data = fetch_M1(1,1);
load('../../../../Kelp_data/Summer2025/Rooker/M1/L0/ADCP/ADCP_M1_config.mat')
Config.dt = 1/double(Config.Burst_SamplingRate);
Config.blank = Config.Burst_BlankingDistance;
Config.binSize = Config.Burst_CellSize;
Config.SN = 'M1 SIG1K';
fn = fieldnames(Data);
for i = 1:numel(fn)
    eval([fn{i} ' = Data.(fn{i});']);
end
save('../../../../Kelp_data/Summer2025/Rooker/Release1/L0/KELP1_M1ADCP_L0.mat', fn{:}, 'Config')
clear

%% M1 R2
disp('M1 R2')
Data = fetch_M1(2,1);
load('../../../../Kelp_data/Summer2025/Rooker/M1/L0/ADCP/ADCP_M1_config.mat')
Config.dt = 1/double(Config.Burst_SamplingRate);
Config.blank = Config.Burst_BlankingDistance;
Config.binSize = Config.Burst_CellSize;
Config.SN = 'M1 SIG1K';
fn = fieldnames(Data);
for i = 1:numel(fn)
    eval([fn{i} ' = Data.(fn{i});']);
end
save('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/KELP2_M1ADCP_L0.mat', fn{:}, 'Config')
clear

%% M2 R1
disp('M2 R1')
Data = fetch_M1(1,2);
load('../../../../Kelp_data/Summer2025/Rooker/M2/L0/ADCP/ADCP_M2_config.mat')
Config.dt = 1/double(Config.Burst_SamplingRate);
Config.blank = Config.Burst_BlankingDistance;
Config.binSize = Config.Burst_CellSize;
Config.SN = 'M2 SIG1K';
fn = fieldnames(Data);
for i = 1:numel(fn)
    eval([fn{i} ' = Data.(fn{i});']);
end
save('../../../../Kelp_data/Summer2025/Rooker/Release1/L0/KELP1_M2ADCP_L0.mat', fn{:}, 'Config')
clear

%% M2 R2
disp('M2 R2')
Data = fetch_M1(2,2);
load('../../../../Kelp_data/Summer2025/Rooker/M2/L0/ADCP/ADCP_M2_config.mat')
Config.dt = 1/double(Config.Burst_SamplingRate);
Config.blank = Config.Burst_BlankingDistance;
Config.binSize = Config.Burst_CellSize;
Config.SN = 'M2 SIG1K';
fn = fieldnames(Data);
for i = 1:numel(fn)
    eval([fn{i} ' = Data.(fn{i});']);
end
save('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/KELP2_M2ADCP_L0.mat', fn{:}, 'Config')
clear
