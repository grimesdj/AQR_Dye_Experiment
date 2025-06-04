
% making new L0 data from the completed raw.mat

Data = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat");


minAmp = min(Data.a1, min(Data.a2, Data.a3));
minCor = min(c1, min(c2, c3));
flagA = find(minAmp <= 30);
flagC = find(minCor <= 70);
flagind = unique([flagA' flagC']);
flagind; %iterate through bins smh

east(flagind) = NaN;
plot(east(:,1))
