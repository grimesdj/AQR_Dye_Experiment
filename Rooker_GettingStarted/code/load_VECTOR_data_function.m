function load_VECTOR_data_function(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime)

%% load header data
%% data for each field starts at column 39 or 40
hdrFile = sprintf(['%s/%s.hdr'], inputDir, inputFile);
fid = fopen(hdrFile);% open file
while ~feof(fid);
    %grab line
    line = fgetl(fid);
    % 
    % some skipped lines are very short, or just have filename
    if length(line)<39
        continue
    elseif strcmp(line(1),'[')
        continue
    end
    %
    % each line has a field name and value
    string = deblank(line(1:36));
    value = line(39:end);
    %
    % select pertinent fields and allocate variables
    if     strcmp(string,'Number of measurements')
        nsamples = str2num(value);
    elseif strcmp(string,'Number of velocity checksum errors')
        nerrors = str2num(value);
% $$$     elseif strcmp(string,'Measurement/Burst interval')
% $$$         i = strfind(value,'sec');
% $$$         dt = str2num(value(1:i-1));
    elseif strcmp(string,'Number of cells')
        nbins = str2num(value);
    elseif strcmp(string,'Sampling rate')
        i = strfind(value,'Hz');
        sample_rate = str2num(value(1:i-1));
        dt          = 1./sample_rate;
    elseif strcmp(string,'Head frequency')
        i = strfind(value,'kHz');
        head_freq = str2num(value(1:i-1));
    elseif strcmp(string,'Sampling volume')
        i = strfind(value,'cm');
        cff=100;
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        sample_vol = str2num(value(1:i-1))/cff;% cm-->m
    elseif strcmp(string,'Transmit length')
        i = strfind(value,'cm');
        cff=100;
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        transmit_length = str2num(value(1:i-1))/cff;% cm-->m
    elseif strcmp(string,'Receive length')
        i = strfind(value,'m');
        cff=1;
        if isempty(i)
            i = strfind(value,'cm');
            cff=100;
        end
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        receive_length = str2num(value(1:i-1))/cff;% cm-->m
    elseif strcmp(string,'Distance between pings')
        i = strfind(value,'m');
        cff=1;
        if isempty(i)
            i = strfind(value,'cm');
            cff=100;
        end
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        pulse_distance = str2num(value(1:i-1))/cff;% cm-->m
    elseif strcmp(string,'Salinity')
        i = strfind(value,'ppt');
        salinity = str2num(value(1:i-1));% cm-->m
    elseif strcmp(string,'Cell size')
        i = strfind(value,'cm');
        cff=100;
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        binsize = str2num(value(1:i-1))/cff;% cm-->m
    elseif strcmp(string,'Velocity scaling')
        i = strfind(value,'cm');
        cff=100;
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        vel_scale = str2num(value(1:i-1))/cff;% cm/s,mm/s-->m/s
    elseif strcmp(string,'Blanking distance')
        i = strfind(value,'m');
        blank = str2num(value(1:i-1));
    elseif strcmp(string,'Coordinate system')
        coords = value;
    elseif strcmp(string,'Serial number') & ~exist('sn','var')
        sn = deblank(value);
    elseif strcmp(string,'Transformation matrix')
        temp   = textscan(value,'%f %f %f');
        T      = cell2mat(temp);
        line   = fgetl(fid);  value = line(39:end);
        temp   = textscan(value,'%f %f %f');
        T(2,:) = cell2mat(temp);        
        line   = fgetl(fid);  value = line(39:end);
        temp   = textscan(value,'%f %f %f');
        T(3,:) = cell2mat(temp);
        updwn  = input('was the instrument facing up (1) or down (0)?\n');
        if updwn
            T(2:3,:) = -T(2:3,:);
        end
    end
clear line string value i
end

meta_data = struct('SN',sn,'Nsamples',nsamples,'Nerrors',nerrors,'dt',dt, ...
                   'coords',coords,'sample_rate',sample_rate,'velocity_scale',vel_scale,...
                   'sample_volume',sample_vol,'transmit_length',transmit_length,...
                   'receive_length',receive_length,'pulse_distance',pulse_distance,...
                   'salinity',salinity,'head_frequency_kHz',head_freq,'transform_matrix',T);
%
fclose(fid);
%
%
% vhdFile = sprintf(['%s',filesep,'%s.hdr'], inputDir,inputFile);
senFile = sprintf(['%s/%s.sen'], inputDir,inputFile);
fid     = fopen(senFile,'r');
format  = '%f %f %f %f %f %f %s %s %f %f %f %f %f %f %d %d';
sen     = textscan(fid,format);
date    = datenum(sen{3},sen{1},sen{2},sen{4},sen{5},sen{6})+tos/24;
Ndt     = length(date);
% $$$ time    = [0:(meta_data.sample_rate-1)]'/86400+date';
% $$$ time    = reshape(time,meta_data.sample_rate*Ndt,1);
batt_volt= sen{9};
sspeed   = sen{10};
head     = sen{11};
pitch    = sen{12};
roll     = sen{13};
temperature= sen{14};
%
%
% get the pressure/velocity infor from .dat file
datFile = sprintf(['%s/%s.dat'], inputDir,inputFile);
fid     = fopen(datFile,'r');
format  = '%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d';
dat     = textscan(fid,format);
Seconds = (1:meta_data.Nsamples)'./meta_data.sample_rate;
Velocity_X    = dat{3};
Velocity_Y    = dat{4};
Velocity_Z    = dat{5};
Amplitude_Beam1    = dat{6};
Amplitude_Beam2    = dat{7};
Amplitude_Beam3    = dat{8};
SNR1  = dat{9};
SNR2  = dat{10};
SNR3  = dat{11};
Correlation_Beam1    = dat{12};
Correlation_Beam2    = dat{13};
Correlation_Beam3    = dat{14};
Pressure = dat{15};
%
SENseconds = 0:60:60*(Ndt-1);
headInterp = interp1(SENseconds,head,Seconds);
%
TStart = date(1);
%TEnd = date(end);
step = Seconds/86400;
A.Time_sensor = date;
A.Time = TStart+step;

l = find(A.Time>=atmTime(1) & A.Time<=atmTime(2));
A.pressureOffset = mean(Pressure(l(1):l(2)));
A.Pressure = Pressure - A.pressureOffset;

dep = find(A.Time>=depTime(1) & A.Time<=depTime(2));
dep_sensor = find(A.Time_sensor>=depTime(1) & A.Time_sensor<=depTime(2));

nsamples = length(dep);
A.Time = A.Time(dep);
A.Time_sensor = A.Time_sensor(dep_sensor);

%
switch coords
  case {'XYZ','ENU'}
    if strcmp(coords,'XYZ')
        shape = size(Velocity_X(dep));
        BEAM  = inv(T)*[Velocity_X(dep)'; Velocity_Y(dep)'; Velocity_Z(dep)'];
        Velocity_Beam1  = reshape(BEAM(1,:)',shape);
        Velocity_Beam2  = reshape(BEAM(2,:)',shape);
        Velocity_Beam3  = reshape(BEAM(3,:)',shape);
        
        Velocity_X = Velocity_X(dep);
        Velocity_Y = Velocity_Y(dep);
        Velocity_Z = Velocity_Z(dep);
    else        
    Velocity_East = Velocity_X(dep);
    Velocity_North= Velocity_Y(dep);
    Velocity_Up   = Velocity_Z(dep);
    end
  case {'BEA','BEAM'}
    Velocity_Beam1 = Velocity_X(dep);
    Velocity_Beam2 = Velocity_Y(dep);
    Velocity_Beam3 = Velocity_Z(dep);
    shape = size(Velocity_Beam1);
    XYZ= T*[Velocity_Beam1'; Velocity_Beam2'; Velocity_Beam3'];
    Velocity_X = reshape(XYZ(1,:)',shape);
    Velocity_Y = reshape(XYZ(2,:)',shape);
    Velocity_Z = reshape(XYZ(3,:)',shape);
end
%
if ~strcmp(coords,'ENU')
    fixedHead = input('Do you want to manually enter heading? (1=yes, 0=no) \n');
    if fixedHead
        if ~exist('headingOffset','var')
            headingOffset = input('Input a constant heading in degrees magnetic (nautical convention)\n');
        end
        hh = pi*headingOffset/180;
        H = [ cos(hh) sin(hh) 0*hh;...
             -sin(hh) cos(hh) 0*hh;...
              0*hh      0*hh  0*hh+1];
        ENU = H*[Velocity_X(dep)';Velocity_Y(dep)';Velocity_Z(dep)'];
        Velocity_East  = ENU(1,:)';
        Velocity_North = ENU(2,:)';
        Velocity_Up    = ENU(3,:)';    
        
    else
    % rotate to EW, pitch/roll matrices don't work w cabled head
    hh = reshape(pi*(headInterp(dep))/180,1,1,nsamples);
    % pp = reshape(pi*pitch/180,1,1,nsamples);
    % rr = reshape(pi*roll/180,1,1,nsamples);
    H = [ cos(hh) sin(hh) 0*hh;...
         -sin(hh) cos(hh) 0*hh;...
          0*hh      0*hh  0*hh+1];
% $$$     P = [cos(pp) -sin(pp).*sin(rr) -cos(rr).*sin(pp);...
% $$$           0*pp       cos(rr)         -sin(rr)   ;...
% $$$          sin(pp)  sin(rr).*cos(pp)  cos(pp).*cos(rr)];
    for j = 1:nsamples
     ENU = H(:,:,j)*[Velocity_X(j);Velocity_Y(j);Velocity_Z(j)];
     Velocity_East (j) = ENU(1);
     Velocity_North(j) = ENU(2);
     Velocity_Up   (j) = ENU(3);    
    end
    end
end
%
A.Config= meta_data;
sensor_data = struct('date',date,'battery_voltage',batt_volt,'sound_speed',sspeed,'heading',head,'pitch',pitch,'roll',roll,'temperature',temperature);
A.sensor    = sensor_data;
A.Sound_Speed    = sspeed(dep_sensor);
Range = (A.Sound_Speed.^2)/(8*6*1000*1.01);
A.VRange = interp1(A.Time_sensor, Range, A.Time);
A.Seconds   = Seconds(dep);
if fixedHead
    A.fixed_heading = headingOffset;
end
A.Velocity_X    = Velocity_X;
A.Velocity_Y    = Velocity_Y;
A.Velocity_Z    = Velocity_Z;
A.Velocity_Beam1    = Velocity_Beam1;
A.Velocity_Beam2    = Velocity_Beam2;
A.Velocity_Beam3    = Velocity_Beam3;
A.Velocity_East  = Velocity_East;
A.Velocity_North = Velocity_North;
A.Velocity_Up    = Velocity_Up;
A.Amplitude_Beam1    = Amplitude_Beam1(dep);
A.Amplitude_Beam2    = Amplitude_Beam2(dep);
A.Amplitude_Beam3    = Amplitude_Beam3(dep);
A.Correlation_Beam1    = Correlation_Beam1(dep);
A.Correlation_Beam2    = Correlation_Beam2(dep);
A.Correlation_Beam3    = Correlation_Beam3(dep);
A.SNR1  = SNR1(dep);
A.SNR2  = SNR2(dep);
A.SNR3  = SNR3(dep);
A.fname = fileName;


disp('Saving raw data')
save([outputFile,'.mat'],'-struct','A')
%

end