
% Returns a Structure M1 with specified release

function Data = fetch_M1(releaseNum)


files = dir("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\2024_PROCESSED_DATA\M1\L0\ADCP\ADCP_M1_*.mat");

TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

for i = 1:length(files)
    
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        if any((Times.Time >= datenum(TRange.StartTime_UTC_(releaseNum))) & Times.Time<= datenum(TRange.EndTime_UTC_(releaseNum)))
            M1 =  load([files(i).folder, filesep, files(i).name]);
        end
    end
end
Data.time = M1.Time;
Data.east = M1.Velocity_East';
Data.north = M1.Velocity_North';
Data.up = M1.Velocity_Up';
Data.b1 = M1.Velocity_Beam(:,:,1)';
Data.b2 = M1.Velocity_Beam(:,:,2)';
Data.b3 = M1.Velocity_Beam(:,:,3)';
Data.a1 = M1.Amplitude_Beam(:,:,1)';
Data.a2 = M1.Amplitude_Beam(:,:,2)';
Data.a3 = M1.Amplitude_Beam(:,:,3)';
Data.c1 = M1.Correlation_Beam(:,:,1)';
Data.c2 = M1.Correlation_Beam(:,:,2)';
Data.c3 = M1.Correlation_Beam(:,:,3)';


end