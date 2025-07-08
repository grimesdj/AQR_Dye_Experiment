

function A = loadAQD(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime)


% pull data

HRflag = 0;

%% load header data
%% data for each field starts at column 39 or 40
hdrFile = sprintf('%s/%s.hdr', inputDir, inputFile);
fid = fopen(hdrFile);% open file

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
        lag1 = str2num(value(1:i-1));
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
    elseif strcmp(string, 'Extended velocity range')
        disp('Instrument being processed is HR')
        HRflag = 1;
        
    end
clear line string value i;
end

meta_data = struct('SN',sn,'Nsamples',nsamples,'Nerrors',nerrors,'dt',dt, ...
                   'Nbins',nbins,'binSize',binsize,'blank',blank,'coords',coords,'transform_matrix',T);
fclose(fid);
%
%% Extracts date, temp, pressure, Vr, heading, pitch, roll
senFile = sprintf(['%s/%s.sen'], inputDir,inputFile);
fid = fopen(senFile,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %*[^\n]');
fclose(fid);
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

% Here is where to set pressure offset and trim time


if ~exist('atmTime','var')
    disp('pick 2 points bounding when out of water for ATM pressure offset')
    plot(pressure)
    l = ginput(2);
    l = round(l(:,1));
    atmTime = [time(l(1)), time(l(2))];
    fprintf('atmTime = \n')
    fprintf('%s --- %s', datestr(atmTime(1)), datestr(atmTime(2)));
else
    l = find(time>=atmTime(1) & time<=atmTime(2));
end
A.pressureOffset = mean(pressure(l(1):l(2)));
A.Pressure = pressure - A.pressureOffset;
%
if ~exist('depTime','var')
    % now trim the data to when it was in the water
    disp('pick start/end points of deployment')
    l = ginput(2);
    l = round(l(:,1));
    depTime = [time(l(1)), time(l(2))];
    fprintf('depTime = \n')
    fprintf('%s --- %s', datestr(depTime(1)), datestr(depTime(2)));
else
    dep = find(time>=depTime(1) & time<=depTime(2));
end
nsamples = length(dep);
meta_data.Nsamples = nsamples;


A.Config= meta_data;
A.Date  = datestr(time(1));
A.Time  = time(dep); A.Volt = volt;
A.Seconds= (A.Time-A.Time(1))*86400;
A.Sound_Speed= sspeed(dep);A.VRange = vrange(dep);
A.Heading  = heading(dep);A.Pitch = pitch(dep);A.Roll = roll(dep);
A.Pressure  = pressure(dep);A.Temperature = temperature(dep);
A.fileName = fileName;
%%%%%%%

% Heres where we get amp and cor data

%% load beam amplitudes and velocities (may be in beam coords or ENU, see A.config)

disp('Trimming data to deployment time while loading')

a1 = load(strcat(fileName,'.a1'));
[Na,Ma] = size(a1);
if Ma>nbins% sometimes there are extra columns (beam#, ensemble#)
    bin1 = Ma-nbins+1;
else
    bin1 = 1;
end
a2 = load(strcat(fileName,'.a2'));
a3 = load(strcat(fileName,'.a3'));
A.Amplitude_Beam1 = a1(dep,bin1:end);
A.Amplitude_Beam2 = a2(dep,bin1:end);
A.Amplitude_Beam3 = a3(dep,bin1:end);
%
v1 = load(strcat(fileName,'.v1'));
v2 = load(strcat(fileName,'.v2'));
v3 = load(strcat(fileName,'.v3'));

A.Velocity_X = v1(dep,bin1:end);
A.Velocity_Y = v2(dep,bin1:end);
A.Velocity_Z = v3(dep,bin1:end);

switch coords
  case {'XYZ','ENU'}
        if strcmp(coords,'XYZ')
        shape = size(A.Velocity_X);
        BEAM  = inv(T)*[A.Velocity_X(:)'; A.Velocity_Y(:)'; A.Velocity_Z(:)'];
        A.Velocity_Beam1  = reshape(BEAM(1,:)',shape);
        A.Velocity_Beam2  = reshape(BEAM(2,:)',shape);
        A.Velocity_Beam3  = reshape(BEAM(3,:)',shape);
        end        
    A.Velocity_East  = A.Velocity_X;
    A.Velocity_North = A.Velocity_Y;
    A.Velocity_Up    = A.Velocity_Z;
  case {'BEA'}
    b1 = v1(dep,bin1:end);
    b2 = v2(dep,bin1:end);
    b3 = v3(dep,bin1:end);
    shape = size(b1);
    XYZ= T*[b1(:)'; b2(:)'; b3(:)'];
    A.Velocity_X = reshape(XYZ(1,:)',shape);
    A.Velocity_Y = reshape(XYZ(2,:)',shape);
    A.Velocity_Z = reshape(XYZ(3,:)',shape);
    A.Velocity_Beam1 = b1;
    A.Velocity_Beam2 = b2;
    A.Velocity_Beam3 = b3;
end
%
if ~strcmp(coords,'ENU')
    % rotate to EW, need to work out the pitch/roll matrices
    hh = reshape(pi*(A.Heading-90)/180,1,1,nsamples);
    pp = reshape(pi*A.Pitch/180,1,1,nsamples);
    rr = reshape(pi*A.Roll/180,1,1,nsamples);
    H = [ cos(hh), sin(hh), 0*hh;...
         -sin(hh), cos(hh), 0*hh;...
          0*hh,      0*hh,  0*hh+1];
    P = [cos(pp), 0*pp, -sin(pp);...
          0*pp,  0*pp+1, 0*pp   ;...
         sin(pp),  0*pp,  cos(pp)];
     
    O = [1+0*rr 0*rr 0*rr;...
        0*rr cos(rr) -sin(rr);...
        0*rr sin(rr) cos(rr)];
     
     
    shape = size(A.Velocity_X);
    for j = 1:nsamples
     R   = H(:,:,j)*P(:,:,j)*O(:,:,j);
     ENU = R*[A.Velocity_X(j,:);A.Velocity_Y(j,:);A.Velocity_Z(j,:)];
     A.Velocity_East  (j,:) = ENU(1,:);
     A.Velocity_North (j,:) = ENU(2,:);
     A.Velocity_Up    (j,:) = ENU(3,:);    
    end
end

% Find correlations
if HRflag
    %
    c1 = load(strcat(fileName,'.c1'));
    c2 = load(strcat(fileName,'.c2'));
    c3 = load(strcat(fileName,'.c3'));
    A.Correlation_Beam1 = c1(dep,bin1:end);
    A.Correlation_Beam2 = c2(dep,bin1:end);
    A.Correlation_Beam3 = c3(dep,bin1:end);

else
    Correlation_Beam1 = NaN;
    Correlation_Beam2 = NaN;
    Correlation_Beam3 = NaN;
end



 % vars  = {'Time','Volt','Seconds','Sound_Speed','Heading','Pitch','Roll','Pressure','Temperature','Amplitude_Beam1','Amplitude_Beam2','Amplitude_Beam3','Velocity_X','Velocity_Y','Velocity_Z','Correlation_Beam1','Correlation_Beam2','Correlation_Beam3','Velocity_Beam1','Velocity_Beam2','Velocity_Beam3','Velocity_East','Velocity_North','Velocity_Up'};
 % for jj = 1:length(vars)
 %     eval(['A.',vars{jj},' = A.',vars{jj},'(:,:);'])
 % end

disp('Saving raw data')
save([outputFile,'.mat'],'-struct','A')



end

