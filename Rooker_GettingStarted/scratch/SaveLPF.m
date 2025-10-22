%% Lets Save this LPF Data
clear all; 
close all;

%% File Setup
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/*.mat');
colors = {[0, 0, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0]};

for i = 1:length(files)
    dataCell{i} = load(fullfile(files(i).folder, files(i).name));
    labels{i} = dataCell{i}.Config.SN;
    
    % Determine appropriate bin
    if size(dataCell{i}.Velocity_East, 2) == 1
        dataCell{i}.bin = 1;
    else
        dataCell{i}.bin = round((1.237 - dataCell{i}.Config.blank) / dataCell{i}.Config.binSize);
    
    end
end


%% Averaged Velocities
addpath('../code')
[T_all, U_all, V_all] = LPF(dataCell, labels, colors);

% normalizing vector size
[~, sizeLimit] = size(U_all{2}); % gonna make this a max function later

AQD = [U_all{1}(1, 1:sizeLimit); V_all{1}(1, 1:sizeLimit)];
M1  = [U_all{2}(1, 1:sizeLimit); V_all{2}(1, 1:sizeLimit)];
M2  = [U_all{3}(1, 1:sizeLimit); V_all{3}(1, 1:sizeLimit)];
VEC = [U_all{4}(1, 1:sizeLimit); V_all{4}(1, 1:sizeLimit)];

save('../../../../Kelp_data/Summer2025/Rooker/Release2//*.mat')