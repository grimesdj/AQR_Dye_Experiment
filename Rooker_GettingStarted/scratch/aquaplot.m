
% making new L0 data from the completed raw.mat
function Stats = aquaplot(releaseNum, version, plots, binStart, binEnd)

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat"];

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
elseif strcmp(version, 'Vector')
    releaseNum = releaseNum+3;
    AQD = Data;
    Data = load(inputFiles(releaseNum));    
    rotation
    releaseNum = releaseNum-3;
else
    return
end



if nargin < 4
    binStart = 1;
    binEnd = size(Data.b1,2);
    
    if nargin < 3
        plots = 'on';
    
    
    end    
end


% Limit data to during dye release
TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.time <= datenum(TRange.EndTime_UTC_(releaseNum)));

% Collect Stats
Pres = [mean(Data.pressure(dye, :),'all', 'omitnan'), std(Data.pressure(dye, :), 0, 'all', 'omitnan')];
East = [mean(Data.east(dye, binStart:binEnd),'all', 'omitnan'), std(Data.east(dye, binStart:binEnd), 0, 'all', 'omitnan')];
North = [mean(Data.north(dye, binStart:binEnd),'all', 'omitnan'), std(Data.north(dye, binStart:binEnd), 0, 'all', 'omitnan')];
theta = atan2d(East(1), North(1));
Mag = [(East(1)^2 + North(1)^2)^0.5, (East(2)^2 + North(2)^2)^0.5]  ;

Stats = struct('Pres', Pres, 'East', East, 'North', North, 'theta', theta, 'Mag', Mag, 'version', version);

if strcmp(plots, 'on')
    plotRelease_func


end
end
