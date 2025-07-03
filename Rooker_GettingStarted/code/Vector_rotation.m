% converting Vector XYZ to ENU matching the orientation of Aquadopp

% get aquadopp orietation

% % % Lets turn this into a function % % %

%what does rotation return? -> Vector data with updated Parameters
function CorrectVec = Vector_rotation(Data, AQD)

Data.heading = interp1(AQD.time, AQD.heading, Data.time);
Data.pitch = interp1(AQD.time, AQD.pitch, Data.time);
Data.roll = interp1(AQD.time, AQD.roll, Data.time);
overlap = find(Data.time >= AQD.time(1) & Data.time <= AQD.time(end));

for i = 1:length(Data.time)
    if  ~overlap(i)
        
        Data.east(i) = NaN;
        Data.north(i) = NaN;        
        Data.up(i) = NaN;
        
    else

        theta = Data.heading(i) - 90;

        st = sind(theta);
        sp = sind(Data.pitch(i));
        so = sind(Data.roll(i));
        ct = cosd(theta);
        cp = cosd(Data.pitch(i));
        co = cosd(Data.roll(i));
        H =[ct st 0; -st ct 0; 0 0 1];
        P = [cp 0 -sp; 0 1 0; sp 0 cp];
        R = [1 0 0; 0 co -so; 0 so co];
        T = H*P*R;

        coords = [Data.v1(i); Data.v2(i); Data.v3(i)];

        ENU = T * coords;

        Data.east(i) = ENU(1);
        Data.north(i) = ENU(2);
        Data.up(i) = ENU(3);

       CorrectVec = Data;
    end
end
end
