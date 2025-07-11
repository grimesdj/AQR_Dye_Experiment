% converting Vector XYZ to ENU matching the orientation of Aquadopp

% get aquadopp orietation

% % % Lets turn this into a function % % %

% what does rotation return? -> Vector data with updated Parameters
function CorrectVec = Vector_rotation(Data, AQD)

Data.Heading = interp1(AQD.Time, AQD.Heading, Data.Time);
Data.Heading = Data.Heading(~isnan(Data.Heading));

Data.Pitch = interp1(AQD.Time, AQD.Pitch, Data.Time);
Data.Pitch = Data.Pitch(~isnan(Data.Pitch));

Data.Roll = interp1(AQD.Time, AQD.Roll, Data.Time);
Data.Roll = Data.Roll(~isnan(Data.Roll));

Data.Velocity_East = Data.Velocity_X;
Data.Velocity_North = Data.Velocity_Y;
Data.Velocity_Up = Data.Velocity_Z;
%overlap = find(Data.Time >= AQD.Time(1) & Data.Time <= AQD.Time(end));

for i = 1:length(Data.Heading)
    % if  ~overlap(i)
    % 
    %     Data.Velocity_East(i) = NaN;
    %     Data.Velocity_North(i) = NaN;        
    %     Data.Velocity_Up(i) = NaN;
    % 
    % else

        theta = Data.Heading(i) - 90;

        st = sind(theta);
        sp = sind(Data.Pitch(i));
        so = sind(Data.Roll(i));
        ct = cosd(theta);
        cp = cosd(Data.Pitch(i));
        co = cosd(Data.Roll(i));
        H =[ct st 0; -st ct 0; 0 0 1];
        P = [cp 0 -sp; 0 1 0; sp 0 cp];
        R = [1 0 0; 0 co -so; 0 so co];
        T = H*P*R;

        coords = [Data.Velocity_X(i); Data.Velocity_Y(i); Data.Velocity_Z(i)];

        ENU = T * coords;

        Data.Velocity_East(i) = ENU(1);
        Data.Velocity_North(i) = ENU(2);
        Data.Velocity_Up(i) = ENU(3);

       
    %end

end
CorrectVec = Data;
end
