
% Returns a Structure M2 with specified release

function Data = fetch_M2(releaseNum)


files = dir("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\2024_PROCESSED_DATA\M2\L0\ADCP\ADCP_M2_*.mat");

TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");
M2 = struct([]);
for i = 1:length(files)
    
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        if any((Times.Time >= datenum(TRange.StartTime_UTC_(releaseNum))) & Times.Time<= datenum(TRange.EndTime_UTC_(releaseNum)))
            M2 =  load([files(i).folder, filesep, files(i).name]);
        end
        if ~isempty(M2)
            break
        end
    end
end
Data.time = M2.Time;
Data.east = M2.Velocity_East';
Data.north = M2.Velocity_North';
Data.up = M2.Velocity_Up';
Data.b1 = M2.Velocity_Beam(:,:,1)';
Data.b2 = M2.Velocity_Beam(:,:,2)';
Data.b3 = M2.Velocity_Beam(:,:,3)';
Data.a1 = M2.Amplitude_Beam(:,:,1)';
Data.a2 = M2.Amplitude_Beam(:,:,2)';
Data.a3 = M2.Amplitude_Beam(:,:,3)';
Data.c1 = M2.Correlation_Beam(:,:,1)';
Data.c2 = M2.Correlation_Beam(:,:,2)';
Data.c3 = M2.Correlation_Beam(:,:,3)';
Data.pressure = M2.Pressure(2,:)';

end