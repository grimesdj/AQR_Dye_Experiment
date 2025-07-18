%% load header data
%% data for each field starts at column 39 or 40
hdrFile = sprintf(['%s',filesep,'%s.hdr'], inputDir,inputFile);
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
    elseif strcmp(string,'Distance between pings')
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
    elseif strcmp(string,'Salinity') 
        salinity = deblank(value);
    elseif strcmp(string,'Head frequency')
        i = strfind(value,'kHz');        
        head_freq = str2num(value(1:i-1));
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
                   'salinity',salinity,'head_frequency_kHz',head_freq,'transformation_matrix',T);
%
fclose(fid);
%
%
% vhdFile = sprintf(['%s',filesep,'%s.hdr'], inputDir,inputFile);
senFile = sprintf(['%s',filesep,'%s.sen'], inputDir,inputFile);
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
datFile = sprintf(['%s',filesep,'%s.dat'], inputDir,inputFile);
fid     = fopen(datFile,'r');
format  = '%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d';
dat     = textscan(fid,format);
seconds = (1:meta_data.Nsamples)'./meta_data.sample_rate;
v1    = dat{3};
v2    = dat{4};
v3    = dat{5};
a1    = dat{6};
a2    = dat{7};
a3    = dat{8};
SNR1  = dat{9};
SNR2  = dat{10};
SNR3  = dat{11};
c1    = dat{12};
c2    = dat{13};
c3    = dat{14};
pressure = dat{15};
%
SENseconds = 0:60:60*(Ndt-1);
headInterp = interp1(SENseconds,head,seconds);
%
switch coords
  case {'XYZ','ENU'}
    if strcmp(coords,'XYZ')
        shape = size(v1);
        BEAM  = inv(T)*[v1(:)'; v2(:)'; v3(:)'];
        b1  = reshape(BEAM(1,:)',shape);
        b2  = reshape(BEAM(2,:)',shape);
        b3  = reshape(BEAM(3,:)',shape);
    end        
    east = v1;
    north= v2;
    up   = v3;
  case {'BEA','BEAM'}
    b1 = v1;
    b2 = v2;
    b3 = v3;
    shape = size(b1);
    XYZ= T*[b1(:)'; b2(:)'; b3(:)'];
    v1 = reshape(XYZ(1,:)',shape);
    v2 = reshape(XYZ(2,:)',shape);
    v3 = reshape(XYZ(3,:)',shape);
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
        ENU = H*[v1';v2';v3'];
        east  = ENU(1,:)';
        north = ENU(2,:)';
        up    = ENU(3,:)';    
        
    else
    % rotate to EW, pitch/roll matrices don't work w cabled head
    hh = reshape(pi*(headInterp)/180,1,1,nsamples);
    % pp = reshape(pi*pitch/180,1,1,nsamples);
    % rr = reshape(pi*roll/180,1,1,nsamples);
    H = [ cos(hh) sin(hh) 0*hh;...
         -sin(hh) cos(hh) 0*hh;...
          0*hh      0*hh  0*hh+1];
% $$$     P = [cos(pp) -sin(pp).*sin(rr) -cos(rr).*sin(pp);...
% $$$           0*pp       cos(rr)         -sin(rr)   ;...
% $$$          sin(pp)  sin(rr).*cos(pp)  cos(pp).*cos(rr)];
    for j = 1:nsamples
     ENU = H(:,:,jj)*[v1(j);v2(j);v3(j)];
     east (j) = ENU(1);
     north(j) = ENU(2);
     up   (j) = ENU(3);    
    end
    end
end
%
A.config= meta_data;
sensor_data = struct('date',date,'battery_voltage',batt_volt,'sound_speed',sspeed,'heading',head,'pitch',pitch,'roll',roll,'temperature',temperature);
A.sensor    = sensor_data;
A.sspeed    = sspeed;
A.seconds   = seconds;
A.pressure  = pressure;
if fixedHead
    A.fixed_heading = headingOffset;
end
A.v1    = v1;
A.v2    = v2;
A.v3    = v3;
A.b1    = b1;
A.b2    = b2;
A.b3    = b3;
A.east  = east;
A.north = north;
A.up    = up;
A.a1    = a1;
A.a2    = a2;
A.a3    = a3;
A.c1    = c1;
A.c2    = c2;
A.c3    = c3;
A.SNR1  = SNR1;
A.SNR2  = SNR2;
A.SNR3  = SNR3;
A.fname = fileName;
%

