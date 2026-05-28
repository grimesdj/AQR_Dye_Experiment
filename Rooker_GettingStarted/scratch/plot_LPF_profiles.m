%% Plot Profiles

% Load M1
M1.ADCP = load('../../../../Kelp_data/data/Release2/L0/ADCP/M1_ADCP.mat');
M1.LPF = load('../../../../Kelp_data/data/Release2/L1/ADCP/M1_LPF_600sec.mat');
M1.Moor = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L1/mooring_M1.mat');

releasenum = 2;
release = string(releasenum);


depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00') ;
            datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];

M1.ADCP.Time = M1.ADCP.Time(M1.ADCP.Time >= depTime(1, releasenum) & M1.ADCP.Time <= depTime(2, releasenum));
M1.Moor.Time = M1.Moor.Time(M1.Moor.Time >= depTime(1, releasenum) & M1.Moor.Time <= depTime(2, releasenum));
Time = interp1(M1.ADCP.Moor, M1.Moor.