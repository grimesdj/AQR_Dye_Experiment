
% Returns a Structure M1 with specified release

function Data = fetch_M1(releaseNum)


files = dir("../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_*.mat");

TRange = readtable("../../../../Kelp_data/info/dye_mixing_cals_and_releases/dye_release_times.csv");

for i = 1:length(files)
    
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        if any((Times.Time >= datenum(TRange.StartTime_UTC_(releaseNum))) & Times.Time<= datenum(TRange.EndTime_UTC_(releaseNum)))
            M1 =  load([files(i).folder, filesep, files(i).name]);
        end
    end
end
Data.Time = M1.Time;
Data.Velocity_East = M1.Velocity_East';
Data.Velocity_North = M1.Velocity_North';
Data.Velocity_Up = M1.Velocity_Up';
Data.Velocity_Beam1 = M1.Velocity_Beam(:,:,1)';
Data.Velocity_Beam2 = M1.Velocity_Beam(:,:,2)';
Data.Velocity_Beam3 = M1.Velocity_Beam(:,:,3)';

Data.Correlation_Minimum = min(M1.Correlation_Beam(:,:,1)', min(M1.Correlation_Beam(:,:,2)', M1.Correlation_Beam(:,:,3)'));
Data.Amplitude_Minimum = min(M1.Amplitude_Beam(:,:,1)', min(M1.Amplitude_Beam(:,:,2)', M1.Amplitude_Beam(:,:,3)'));
Data.Pressure = M1.Pressure(2,:)';
Data.Heading = M1.Heading(:, :)';
Data.Pitch = M1.Pitch(:,:)';
Data.Roll = M1.Roll(:,:)';
end