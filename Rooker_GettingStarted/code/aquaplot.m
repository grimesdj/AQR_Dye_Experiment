
% making new L0 data from the completed raw.mat
function Stats = aquaplot(releaseNum, version, plots, binStart, binEnd)

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\Release1\L0\KELP1_AquadoppHR_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Aquadopp_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Aquadopp_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat"];



if strcmp(version, 'raw')
    
    Data = load(inputFiles(releaseNum));
    Data.Time           = Data.time;
    Data.Velocity_X     = Data.v1;
    Data.Velocity_Y     = Data.v2;
    Data.Velocity_Z     = Data.v3;
    Data.Velocity_East  = Data.east;
    Data.Velocity_North = Data.north;
    Data.Velocity_Up    = Data.up;
    Data.Velocity_Beam1 = Data.b1;
    Data.Velocity_Beam2 = Data.b2;
    Data.Velocity_Beam3 = Data.b3;
    
    Data.Amplitude_Minimum = min(Data.a1, min(Data.a2, Data.a3));
    Data.Correlation_Minimum = min(Data.c1, min(Data.c2, Data.c3));
    
elseif strcmp(version, 'L0')
    Data = load(inputFiles(releaseNum+3));
elseif strcmp(version, 'M1')
    Data = fetch_M1(releaseNum);
elseif strcmp(version, 'Vector')
    releaseNum = releaseNum+6;
    AQD = Data;
    Data = load(inputFiles(releaseNum));    
    Data = Vector_rotation(Data, AQD);
    releaseNum = releaseNum-3;
    Data.Time           = Data.time;
    Data.Velocity_X     = Data.v1;
    Data.Velocity_Y     = Data.v2;
    Data.Velocity_Z     = Data.v3;
    Data.Velocity_East  = Data.east;
    Data.Velocity_North = Data.north;
    Data.Velocity_Up    = Data.up;
    Data.Velocity_Beam1 = Data.b1;
    Data.Velocity_Beam2 = Data.b2;
    Data.Velocity_Beam3 = Data.b3;
    
    Data.Amplitude_Minimum = min(Data.a1, min(Data.a2, Data.a3));
    Data.Correlation_Minimum = min(Data.c1, min(Data.c2, Data.c3));
else
    return
end



if nargin < 4
    binStart = 1;
    binEnd = size(Data.Velocity_Beam1,2);
    
    if nargin < 3
        plots = 'on';
    
    
    end    
end


% Limit data to during dye release
TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.Time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.Time <= datenum(TRange.EndTime_UTC_(releaseNum)));

% Collect Stats
Pres = [mean(Data.Pressure(dye, :),'all', 'omitnan'), std(Data.Pressure(dye, :), 0, 'all', 'omitnan')];
 East = [mean(Data.Velocity_East(dye, binStart:binEnd),'all', 'omitnan'), std(Data.Velocity_East(dye, binStart:binEnd), 0, 'all', 'omitnan')];
 North = [mean(Data.Velocity_North(dye, binStart:binEnd),'all', 'omitnan'), std(Data.Velocity_North(dye, binStart:binEnd), 0, 'all', 'omitnan')];
 theta = atan2d(East(1), North(1));
 Mag = [(East(1)^2 + North(1)^2)^0.5, (East(2)^2 + North(2)^2)^0.5]  ;
% 
 Stats = struct('Pres', Pres, 'East', East, 'North', North, 'theta', theta, 'Mag', Mag, 'version', version);


if strcmp(plots, 'on')
    plotRelease_Standardized


end
end
