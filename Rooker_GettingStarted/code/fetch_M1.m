
% Returns a Structure M1 with specified release

function Data = fetch_M1(releaseNum, mooringID)


files = dir("../../../../Kelp_data/Summer2025/Rooker/M1/L0/ADCP/ADCP_M1_*.mat");
files2 = dir("../../../../Kelp_data/Summer2025/Rooker/M2/L0/ADCP/ADCP_M2_*.mat");

if mooringID == 2
    files = files2;
end

TRange = readtable("../../../../Kelp_data/info/dye_mixing_cals_and_releases/dye_release_times.csv");
count = 0;
for i = 1:length(files)
    ind = sprintf('file%d', count+1);
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        if any((Times.Time >= datenum(TRange.StartTime_UTC_(releaseNum+3))) & Times.Time<= datenum(TRange.EndTime_UTC_(releaseNum+3)))
            Moor.(ind) =  load([files(i).folder, filesep, files(i).name]);
            count = count+1;
        end
    end
end


Data.Time                   = [];
Data.Velocity_East          = [];
Data.Velocity_North         = [];
Data.Velocity_Up            = [];
Data.Velocity_Beam1         = [];
Data.Velocity_Beam2         = [];
Data.Velocity_Beam3         = [];
Data.Correlation_Minimum    = [];
Data.Amplitude_Minimum      = [];
Data.Pressure               = [];
Data.Heading                = [];
Data.Pitch                  = [];
Data.Roll                   = [];

for i = 1:count
    ind = sprintf('file%d', i);

    Data.Time                   = [Data.Time;               Moor.(ind).Time'];
    Data.Velocity_East          = [Data.Velocity_East;      Moor.(ind).Velocity_East(:,:)'];
    Data.Velocity_North         = [Data.Velocity_North;     Moor.(ind).Velocity_North(:,:)'];
    Data.Velocity_Up            = [Data.Velocity_Up;        Moor.(ind).Velocity_Up(:,:)'];
    Data.Velocity_Beam1         = [Data.Velocity_Beam1;     Moor.(ind).Velocity_Beam(:,:,1)'];
    Data.Velocity_Beam2         = [Data.Velocity_Beam2;     Moor.(ind).Velocity_Beam(:,:,2)'];
    Data.Velocity_Beam3         = [Data.Velocity_Beam3;     Moor.(ind).Velocity_Beam(:,:,3)'];

    Data.Correlation_Minimum    = [Data.Correlation_Minimum;min(Moor.(ind).Correlation_Beam(:,:,1)', min(Moor.(ind).Correlation_Beam(:,:,2)', Moor.(ind).Correlation_Beam(:, :,3)'))];
    Data.Amplitude_Minimum      = [Data.Amplitude_Minimum;  min(Moor.(ind).Amplitude_Beam(:, :,1)', min(Moor.(ind).Amplitude_Beam(:, :,2)', Moor.(ind).Amplitude_Beam(:, :,3)'))];
    Data.Pressure               = [Data.Pressure;           Moor.(ind).Pressure(2,:)'];
    Data.Heading                = [Data.Heading;            Moor.(ind).Heading(:, :)'];
    Data.Pitch                  = [Data.Pitch;              Moor.(ind).Pitch(:, :)'];
    Data.Roll                   = [Data.Roll;               Moor.(ind).Roll(:, :)'];
   
end

dye = find(Data.Time >= datenum(TRange.StartTime_UTC_(releaseNum+3)) & Data.Time <= datenum(TRange.EndTime_UTC_(releaseNum+3)));

Data.Time                   = Data.Time(dye);
Data.Velocity_East          = Data.Velocity_East(dye,:);
Data.Velocity_North         = Data.Velocity_North(dye,:);
Data.Velocity_Up            = Data.Velocity_Up(dye,:);
Data.Velocity_Beam1         = Data.Velocity_Beam1(dye,:);
Data.Velocity_Beam2         = Data.Velocity_Beam2(dye,:);
Data.Velocity_Beam3         = Data.Velocity_Beam3(dye,:);
Data.Correlation_Minimum   = Data.Correlation_Minimum(dye,:);
Data.Amplitude_Minimum      = Data.Amplitude_Minimum(dye,:);
Data.Pressure               = Data.Pressure(dye,:);
Data.Heading                = Data.Heading(dye,:);
Data.Pitch                  = Data.Pitch(dye,:);
Data.Roll                   = Data.Roll(dye,:);
