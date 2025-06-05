
% making new L0 data from the completed raw.mat

Data = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP2_AquadoppHR_raw.mat");

find(datenum('03-Jul-2024 18:38:00'))

return

minAmp = min(Data.a1, min(Data.a2, Data.a3));
minCor = min(Data.c1, min(Data.c2, Data.c3));
flagA = find(minAmp <= 30);
flagC = find(minCor <= 70);
flagind = unique([flagA' flagC']);
flagind; %iterate through bins smh
    
Data.a1(flagind) = NaN;
Data.a2(flagind) = NaN;
Data.a3(flagind) = NaN;
Data.b1(flagind) = NaN;
Data.b2(flagind) = NaN;
Data.b3(flagind) = NaN;
Data.c1(flagind) = NaN;
Data.c2(flagind) = NaN;
Data.c3(flagind) = NaN;
Data.east(flagind) = NaN;
Data.north(flagind) = NaN;
Data.up(flagind) = NaN;

plotRelease(0, Data)
