
% making new L0 data from the completed raw.mat

han = '';
inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat"];

releaseNum = input("Release Number: ");
version = input("Version: ");

Data = load(inputFiles(releaseNum));

if strcmp(version, 'raw')
    % do nothing
elseif strcmp(version, 'L0')
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
    han = 'L0';
else
    return
end


% Limit data to during dye release
TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)));

plotRelease_func
