


function Data = fetch_sig1k(inputDir, startTime, endTime)
%%
% 
% USAGE: Opens signature1000 ADCP files based on deployment times
%
%   inputDir: Directory where the .mat's can be found
%   statTime: Start of deployment (I think it can be in any date format)
%   endTIme: End of deployment
%
%   % Returns: Structure 'Data' with sig1K data
%
%%

files = dir([inputDir, '*.mat']);

for i = 1:length(files)
    
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        if any((Times.Time >= datenum(startTime)) & Times.Time<= datenum(endTime))
            Sig =  load([files(i).folder, filesep, files(i).name]);
        end
    end
end
Data.Time = Sig.Time;
Data.Velocity_East = Sig.Velocity_East';
Data.Velocity_North = Sig.Velocity_North';
Data.Velocity_Up = Sig.Velocity_Up';
Data.Velocity_Beam1 = Sig.Velocity_Beam(:,:,1)';
Data.Velocity_Beam2 = Sig.Velocity_Beam(:,:,2)';
Data.Velocity_Beam3 = Sig.Velocity_Beam(:,:,3)';

Data.Correlation_Minimum = min(Sig.Correlation_Beam(:,:,1)', min(Sig.Correlation_Beam(:,:,2)', Sig.Correlation_Beam(:,:,3)'));
Data.Amplitude_Minimum = min(Sig.Amplitude_Beam(:,:,1)', min(Sig.Amplitude_Beam(:,:,2)', Sig.Amplitude_Beam(:,:,3)'));
Data.Pressure = Sig.Pressure(2,:)';
Data.Heading = Sig.Heading(:, :)';
Data.Pitch = Sig.Pitch(:,:)';
Data.Roll = Sig.Roll(:,:)';
end