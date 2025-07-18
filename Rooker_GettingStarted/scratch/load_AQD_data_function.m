


function A = load_AQD_data_function(inputDir, inputFile, fileName, tos)

%% load header data
%% data for each field starts at column 39 or 40
hdrFile = sprintf('%s/%s.hdr', inputDir, inputFile);
fid = fopen(hdrFile);% open file

lag1 = NaN;
lag2 = NaN;

while ~feof(fid);
    %grab line
    line = fgetl(fid);
    %some skipped lines are very short
    if length(line)<26
        continue
    end
    %
    % each line has a field name and value
    string = deblank(line(1:26));
    value = line(39:end);
    %
    % select pertinent fields and allocate variables
    if     strcmp(string,'Number of measurements')
        nsamples = str2num(value);
    elseif strcmp(string,'Number of checksum errors')
        nerrors = str2num(value);
    elseif strcmp(string,'Measurement/Burst interval')
        i = strfind(value,'sec');
        dt = str2num(value(1:i-1));
    elseif strcmp(string, 'Pulse distance (Lag1)')
        i = strfind(value, 'm');
        lag1 = str2num(value(1:i-1))
    elseif strcmp(string, 'Pulse distance (Lag2)')
        i = strfind(value, 'm');
        lag2 = str2num(value(1:i-1));
    elseif strcmp(string,'Number of cells')
        nbins = str2num(value);
    elseif strcmp(string,'Cell size')
        i = strfind(value,'cm');
        cff=100;
        if isempty(i)
            i = strfind(value,'mm');
            cff=1000;
        end
        binsize = str2num(value(1:i-1))/cff;% cm-->m
    elseif strcmp(string,'Blanking distance')
        i = strfind(value,'m');
        blank = str2num(value(1:i-1));
    elseif strcmp(string,'Coordinate system')
        coords = value(1:3);
    elseif strcmp(string,'Serial number') & ~exist('sn','var')
        sn = deblank(value);
    elseif strcmp(string,'Transformation matrix')
        T(1,1:3)  = str2num(value);
        line      = fgetl(fid);
        value     = line(39:end);
        T(2,1:3)  = str2num(value);
        line      = fgetl(fid);
        value     = line(39:end);
        T(3,1:3)  = str2num(value);
    end
clear line string value i
end


meta_data = struct('SN',sn,'Nsamples',nsamples,'Nerrors',nerrors,'dt',dt, ...
                   'Nbins',nbins,'binSize',binsize,'blank',blank,'coords',coords,'transform_matrix',T);
fclose(fid);
%
%% Extracts date, temp, pressure, Vr, heading, pitch, roll
senFile = sprintf(['%s/%s.sen'], inputDir,inputFile);
fid = fopen(senFile,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %*[^\n]');
fclose(fid)
A1 = [A{:,1}, A{:,2}, A{:,3}, A{:,4}, A{:,5}, A{:,6}, A{:,7}, A{:,8}, A{:,9}, A{:,10}, A{:,11}, A{:,12} A{:,13}, A{:,14}, A{:,15} A{:,16} A{:,17}];
time = datenum(A1(:,3),A1(:,1),A1(:,2),A1(:,4),A1(:,5),A1(:,6))+tos/24;
volt = A1(:,11);
sspeed = A1(:,12);
vrange = (sspeed.^2)/(8*1000^2*lag1);
heading = A1(:,13);pitch = A1(:,14);roll = A1(:,15);
pressure = A1(:,16);
temperature = A1(:,17);
clear A A1
%
A.config= meta_data;
A.date  = datestr(time(1));
A.time  = time; A.volt = volt;
A.seconds= (time-time(1))*86400;
A.sspeed= sspeed;A.vrange = vrange;
A.heading  = heading;A.pitch = pitch;A.roll = roll;
A.pressure  = pressure;A.temperature = temperature;
A.fname = fileName;

%
%% load beam amplitudes and velocities (may be in beam coords or ENU, see A.config)
a1 = load(strcat(fileName,'.a1'));
[Na,Ma] = size(a1);
if Ma>nbins% sometimes there are extra columns (beam#, ensemble#)
    bin1 = Ma-nbins+1;
else
    bin1 = 1;
end
a2 = load(strcat(fileName,'.a2'));
a3 = load(strcat(fileName,'.a3'));
A.a1 = a1(:,bin1:end);
A.a2 = a2(:,bin1:end);
A.a3 = a3(:,bin1:end);
%
v1 = load(strcat(fileName,'.v1'));
v2 = load(strcat(fileName,'.v2'));
v3 = load(strcat(fileName,'.v3'));
switch coords
  case {'XYZ','ENU'}
    A.v1 = v1(:,bin1:end);
    A.v2 = v2(:,bin1:end);
    A.v3 = v3(:,bin1:end);
    if strcmp(coords,'XYZ')
        shape = size(A.v1);
        BEAM  = inv(T)*[A.v1(:)'; A.v2(:)'; A.v3(:)'];
        A.b1  = reshape(BEAM(1,:)',shape);
        A.b2  = reshape(BEAM(2,:)',shape);
        A.b3  = reshape(BEAM(3,:)',shape);
    end        
    A.east = A.v1;
    A.north= A.v2;
    A.up   = A.v3;
  case {'BEA'}
    b1 = v1(:,bin1:end);
    b2 = v2(:,bin1:end);
    b3 = v3(:,bin1:end);
    shape = size(b1);
    XYZ= T*[b1(:)'; b2(:)'; b3(:)'];
    A.v1 = reshape(XYZ(1,:)',shape);
    A.v2 = reshape(XYZ(2,:)',shape);
    A.v3 = reshape(XYZ(3,:)',shape);
    A.b1 = b1;
    A.b2 = b2;
    A.b3 = b3;
end
%
if ~strcmp(coords,'ENU')
    % rotate to EW, need to work out the pitch/roll matrices
    hh = reshape(pi*(heading-90)/180,1,1,nsamples);
    pp = reshape(pi*pitch/180,1,1,nsamples);
    rr = reshape(pi*roll/180,1,1,nsamples);
    H = [ cos(hh), sin(hh), 0*hh;...
         -sin(hh), cos(hh), 0*hh;...
          0*hh,      0*hh,  0*hh+1];
    P = [cos(pp), 0*pp, -sin(pp);...
          0*pp,  0*pp+1, 0*pp   ;...
         sin(pp),  0*pp,  cos(pp)];
     
    O = [1+0*rr 0*rr 0*rr;...
        0*rr cos(rr) -sin(rr);...
        0*rr sin(rr) cos(rr)];
     
     
    shape = size(A.v1);
    for j = 1:nsamples
     R   = H(:,:,j)*P(:,:,j)*O(:,:,j);
     ENU = R*[A.v1(j,:);A.v2(j,:);A.v3(j,:)];
     A.east (j,:) = ENU(1,:);
     A.north(j,:) = ENU(2,:);
     A.up   (j,:) = ENU(3,:);    
    end
end
%
c1 = load(strcat(fileName,'.c1'));
c2 = load(strcat(fileName,'.c2'));
c3 = load(strcat(fileName,'.c3'));
A.c1 = c1(:,bin1:end);
A.c2 = c2(:,bin1:end);
A.c3 = c3(:,bin1:end);
%
A.dbins = blank + binsize*[1:nbins];
% 

end
