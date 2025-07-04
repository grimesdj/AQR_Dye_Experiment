
% making new L0 data from the completed raw.mat
function Stats = aquastat(releaseNum, version)

inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
               "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat"];

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
    
elseif strcmp(version, 'M1')
    Data = fetch_M1(releaseNum);
    Data.east = Data.east(:,1:12); 
    Data.north = Data.north(:, 1:12);
else
    return
end


% Limit data to during dye release
TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.time <= datenum(TRange.EndTime_UTC_(releaseNum)));


statsPres = [mean(Data.pressure(dye, :),'all', 'omitnan'), std(Data.pressure(dye, :), 0, 'all', 'omitnan')]
statsEast = [mean(Data.east(dye, :),'all', 'omitnan'), std(Data.east(dye, :), 0, 'all', 'omitnan')]
statsNorth = [mean(Data.north(dye, :),'all', 'omitnan'), std(Data.north(dye, :), 0, 'all', 'omitnan')]
theta = atan2d(statsEast(1), statsNorth(1))
mag = (statsEast(1)^2 + statsNorth(1)^2)^0.5
figure, plot(Data.time(dye), Data.pressure(dye), '.')

%plotRelease_func

end










