% mainVEC.m
% 
%   USAGE: loads Vector data from textfiles into .mat files
% 
%   (script) -> requires user to enter:
%       inputDir = Directory where textfiles are
%       inputFile = file root for AQD files
%       outputDir = Directory to save raw .mat
%       L0Dir = Directory to save L0 .mat
%       atmTime = two datetimes when instrument was in the air
%       depTime = start and end times of deployment
%       

addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\code'


clear all
close all

% Release number
releasenum = 2; 
release = string(releasenum);

% Does the signal need to be unwrapped?
unwrap = 0; % 1 for unwrap, 0 for not unwrap




% Enter input /directory/ and fileName root without file extension
inputDir  = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', "Release" + release, 'raw');
inputFile = "KELP" + release + "_Vector";
headingOffset = 331.4;% based on heading from AquodoppHR_KELP1
fileName  = fullfile(inputDir,inputFile);
% Enter raw output /directory/ and fileName without .mat
outputDir = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', "Release" + release, 'raw');
outputName= inputFile + "_raw";
outputFile= fullfile(outputDir, outputName);
% Enter processed output fileName without .mat
L0Dir   = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', "Release" + release, 'L0', 'ADCP');
L0Name  = 'VEC_ADCP';
% Enter path to save figures
figDir = fullfile(inputDir, '..', 'figures');
if ~exist(figDir,'dir'), mkdir(figDir), end
%
% Enter time-period for estimating the atmospheric pressure offset and deployment
if releasenum == 1
    atmTime = [datenum('03-Jul-2024 17:30:00'), datenum('03-Jul-2024 18:10:00')];
    depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00')];
elseif releasenum == 2
    atmTime = [datenum('08-Jul-2024 16:00:00'), datenum('08-Jul-2024 16:30:00')];
    depTime  = [datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];
end

%
% time offset if necessary
tos = 0;
%
% returns structure A with all vector data
%if ~exist([outputDir,'/',outputName,'.mat'],'file')
disp('Loading raw data...')
loadVEC(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime);
    
disp('Generating L0 data...')
L0_Vector(outputFile, L0Dir, L0Name, unwrap, releasenum);

%% Functions

function loadVEC(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime)
% 
%   USAGE: loadVEC(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime)
%       inputDir   = Directory where textfiles are
%       inputFile  = file root for AQD files
%       fileName   = Directory and file root
%       outputFile = File to save raw .mat
%       atmTime    = [Datenum Datenum] in air
%       depTime    = [Start_time , End_time] in water
%       tos        = time offset for time zone correction

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
    elseif strcmp(string, 'Nominal velocity range')
        i = strfind(value,'m/s');
        VRange = str2num(value(1:i-1));
        % disp('Nominal Velocity Range Seems to be in dm/s rather than m/s? -> dividing by 10')
        % VRange = VRange / 10;
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
        ENU = H*[Velocity_X';Velocity_Y';Velocity_Z'];
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
Range = ((A.Sound_Speed).^2)./(8*A.Config.head_frequency_kHz*1000*0.1);
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
save([outputFile + ".mat"],'-struct','A')
%


% Summary figure
figure
ax1 = subplot(3, 1, 1);
plot(ax1, A.Time, A.Velocity_Beam1, '.')
ylabel({'Beam 1', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
hold on 
vr = plot(ax1, A.Time, A.VRange, 'r', A.Time, -1*A.VRange, 'r');
grid minor
legend(vr, 'Maximum Velocity Range', 'Location','northeast')


ax2 = subplot(3, 1, 2);
plot(ax2, A.Time, A.Velocity_Beam2, '.')
ylabel({'Beam 2', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor
hold on 
plot(ax2, A.Time, A.VRange, 'r', A.Time, -1*A.VRange, 'r');

ax3 = subplot(3, 1, 3);
plot(ax3, A.Time, A.Velocity_Beam3, '.')
ylabel({'Beam 3', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
grid minor
hold on 
plot(ax3, A.Time, A.VRange, 'r', A.Time, -1*A.VRange, 'r');

sgtitle('VEC Raw Beam Velocities', 'Fontsize', 25)

end

% EOF

function L0_Vector(outputFile, L0Dir, L0Name, unwrap, releasenum)
% 
%   USAGE: L0_Vector(outputFile, L0Dir, L0Name)
%       outputFile = folder path and filename (no extension) to raw data
%       L0Dir = directory to save finished L0 data
%       L0Name = L0 filename (without extension) to be saved
% 
%       takes raw Vector data and performs L0 QA/QC
% 

A = load([outputFile + ".mat"]);
% temporary addpath for testing :(
addpath '/Users/jasonrooker/Library/CloudStorage/OneDrive-UNC-Wilmington/Kelp_repo/AQR_Dye_Experiment/Rooker_GettingStarted/code'
%addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\Rooker_GettingStarted\code'

%fprintf('\n============================\n\nDo you want to unwrap beam Velocities?')
%unwrap = input('\n(1 = yes; 0 = no)');
if unwrap == 1
    fprintf('unwrapping!')
    for beam = 1:3
        Vel = (sprintf('Velocity_Beam%d',beam));
        vwrap = reshape(A.(Vel), [], 24);
        Vr = reshape(A.VRange, [], 24);
        Vr = Vr.*0.01;
        [unwrap, A.(sprintf('Suspect_Beam%d', beam))] = unwrap_VEC(vwrap, Vr);
        A.(Vel) = unwrap(:);
    end

    A.Correlation_Beam1(find(A.Suspect_Beam1)) = 999;
    A.Correlation_Beam2(find(A.Suspect_Beam2)) = 999;
    A.Correlation_Beam3(find(A.Suspect_Beam3)) = 999;

% rotate unwrapped data!
    % XYZ
    disp('rotating unwrapped data to XYZ')
    shape = size(A.Velocity_Beam1);
    XYZ= A.Config.transform_matrix*[A.Velocity_Beam1(:)'; A.Velocity_Beam2(:)'; A.Velocity_Beam3(:)'];
    A.Velocity_X = reshape(XYZ(1,:)',shape);
    A.Velocity_Y = reshape(XYZ(2,:)',shape);
    A.Velocity_Z = reshape(XYZ(3,:)',shape);

    % ENU
    disp('not rotating unwrapped data to ENU')
   
    
    
    
    % hh = reshape(pi*(A.Heading-90)/180,1,1,A.Config.Nsamples);
    % pp = reshape(pi*A.Pitch/180,1,1,A.Config.Nsamples);
    % rr = reshape(pi*A.Roll/180,1,1,A.Config.Nsamples);
    % H = [ cos(hh), sin(hh), 0*hh;...
    %      -sin(hh), cos(hh), 0*hh;...
    %       0*hh,      0*hh,  0*hh+1];
    % P = [cos(pp), 0*pp, -sin(pp);...
    %       0*pp,  0*pp+1, 0*pp   ;...
    %      sin(pp),  0*pp,  cos(pp)];
    % 
    % O = [1+0*rr 0*rr 0*rr;...
    %     0*rr cos(rr) -sin(rr);...
    %     0*rr sin(rr) cos(rr)];
    % 
    % 
    % 
    % for j = 1:A.Config.Nsamples
    %     R   = H(:,:,j)*P(:,:,j)*O(:,:,j);
    %     ENU = R*[A.Velocity_X(j,:);A.Velocity_Y(j,:);A.Velocity_Z(j,:)];
    %     A.Velocity_East  (j,:) = ENU(1,:)';
    %     A.Velocity_North (j,:) = ENU(2,:)';
    %     A.Velocity_Up    (j,:) = ENU(3,:)';
    % end

end

fprintf('\n============================\n\nDo you want to use external heading?')
headcorrect = input('\n(1 = yes; 0 = no)');
if headcorrect == 1
    %head = load(input('\nEnter path for correct heading:'));
    disp('Using AquaDopp HPR for now')
    if releasenum == 1
        HPRDir = fullfile(L0Dir, '..', '..', 'raw', '*AquadoppHR_raw.mat');
    elseif releasenum == 2
        HPRDir = fullfile(L0Dir, '..', '..', 'raw', '*Aquadopp_raw.mat');
    end
    HPRfiles = dir(HPRDir);
    HPR = load(fullfile(HPRfiles(1).folder, HPRfiles(1).name));
    A = Vector_rotation(A, HPR);
end
% %
% %
% % plot some stuff
% if ~exist('atmTime','var')
%     disp('pick 2 points bounding when out of water for ATM pressure offset')
%     plot(A.pressure)
%     l = ginput(2);
%     l = round(l(:,1));
%     atmTime = [A.time(l(1)), A.time(l(2))];
%     fprintf('atmTime = \n')
%     fprintf('%s --- %s', datestr(atmTime(1)), datestr(atmTime(2)));
% else
%     l = find(A.time>=atmTime(1) & A.time<=atmTime(2));
% end
% A.pressureOffset = mean(A.pressure(l(1):l(2)));
% %
% if ~exist('depTime','var')
%     % now trim the data to when it was in the water
%     disp('pick start/end points of deployment')
%     l = ginput(2);
%     l = round(l(:,1));
%     depTime = [A.time(l(1)), A.time(l(2))];
%     fprintf('depTime = \n')
%     fprintf('%s --- %s', datestr(depTime(1)), datestr(depTime(2)));
% else
%     valid = find(A.time>=depTime(1) & A.time<=depTime(2));
% end
% %
% vars  = {'time','volt','seconds','sspeed','heading','pitch','roll','pressure','temperature','a1','a2','a3','v1','v2','v3','c1','c2','c3','b1','b2','b3','east','north','up'};
% for jj = 1:length(vars)
%     eval(['A.',vars{jj},' = A.',vars{jj},'(valid,:);'])
% end
nsamples = length(A.Time);
%A.dbins = (A.Config.blank + A.Config.binSize*A.Config.Nbins);
%
%A.maxRange = (A.Pressure-A.pressureOffset);
%A.ylims      = [0 min(max(A.maxRange))];
%dum1       = A.maxRange.*ones(1,A.Config.Nbins);
%dum2       = ones(nsamples,1)*A.dbins;
%qcFlag0    =  (dum2<=dum1);
A.Amplitude_Minimum   = min(A.Amplitude_Beam1, min(A.Amplitude_Beam2, A.Amplitude_Beam3));
A.Correlation_Minimum = min(A.Correlation_Beam1, min(A.Correlation_Beam2, A.Correlation_Beam3, 'omitnan'), 'omitnan');   
A.qcFlag              =  double( A.Amplitude_Minimum > 20 & A.Correlation_Minimum > 40 );
%
Time = datetime(A.Time,'convertFrom','datenum');
%
%
disp('Applying qcFlag to NaN data A < 20 and C < 40')
%

A.Velocity_East = A.Velocity_East.*A.qcFlag;
A.Velocity_North = A.Velocity_North.*A.qcFlag;
A.Velocity_Up = A.Velocity_Up.*A.qcFlag;

% A.Velocity_East(~A.qcFlag')=nan;
% A.Velocity_North(~A.qcFlag')=nan;
% A.Velocity_Up(~A.qcFlag')=nan;
%
%
%
A.Velocity_Beam1 = A.Velocity_Beam1.*A.qcFlag;
A.Velocity_Beam2 = A.Velocity_Beam2.*A.qcFlag;
A.Velocity_Beam3 = A.Velocity_Beam3.*A.qcFlag;

% A.Velocity_Beam1(~A.qcFlag')=nan;
% A.Velocity_Beam2(~A.qcFlag')=nan;
% A.Velocity_Beam3(~A.qcFlag')=nan;



%
A.Velocity_X = A.Velocity_X.*A.qcFlag;
A.Velocity_Y = A.Velocity_Y.*A.qcFlag;
A.Velocity_Z = A.Velocity_Z.*A.qcFlag;
%
%

%
% add the config info to the structure A to quick save as netcdf4
fieldNames = fields(A.Config);
originalFields = fields(A);

for j = 1:length(fieldNames)
A.(fieldNames{j}) = A.Config.(fieldNames{j});
end
A = orderfields(A,cat(1,fieldNames,originalFields));
ncfile = fullfile(L0Dir,L0Name + ".nc");
if exist(ncfile,'file')
    delete(ncfile)
end

struct2nc(A,ncfile,'NETCDF4');
%
%

% Make the names match convention
% A.Time = A.time;
% A.Config = A.config;
% A.Pressure = A.pressure;
%A.Bins = A.dbins;
% fieldsToKeep = {'Time', 'Velocity_East', 'Velocity_North', 'Velocity_Up', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Velocity_Beam1', 'Velocity_Beam2', 'Velocity_Beam3', 'Amplitude_Minimum', 'Correlation_Minimum', 'Config', 'Pressure','Heading', 'Pitch', 'Roll', 'Bins', 'Temperature', 'u1', 'u2', 'u3'};
% L0 = rmfield(A, setdiff(fieldnames(A), fieldsToKeep));

disp('Saving L0 data')
out = fullfile(L0Dir, L0Name + ".mat");
save(out,'-struct','A')

% Summary figure
figure
ax1 = subplot(3, 1, 1);
plot(ax1, A.Time, A.Velocity_Beam1, '.')
ylabel({'Beam 1', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor

ax2 = subplot(3, 1, 2);
plot(ax2, A.Time, A.Velocity_Beam2, '.')
ylabel({'Beam 2', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor

ax3 = subplot(3, 1, 3);
plot(ax3, A.Time, A.Velocity_Beam3, '.')
ylabel({'Beam 3', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
grid minor

sgtitle('VEC L0 Beam Velocities', 'Fontsize', 25)


function CorrectVec = Vector_rotation(Data, AQD)
% converting Vector XYZ to ENU matching the orientation of Aquadopp

% get aquadopp orietation

% % % Lets turn this into a function % % %

% what does rotation return? -> Vector data with updated Parameters

Data.Heading = interp1(AQD.Time, AQD.heading, Data.Time);
Data.Heading = Data.Heading(~isnan(Data.Heading));

Data.Pitch = interp1(AQD.Time, AQD.pitch, Data.Time);
Data.Pitch = Data.Pitch(~isnan(Data.Pitch));

Data.Roll = interp1(AQD.Time, AQD.roll, Data.Time);
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


function [v_unwrap, suspect_pts] = unwrap_VEC(v_wrapped, Vr)

disp('using Shcherbina et al 2018 unwrapper')

[nbins, nt] = size(v_wrapped);
v_unwrap = v_wrapped;
filt_diff = (v_wrapped - medfilt1(v_wrapped,150))./Vr;
filt_diff1 = filt_diff;
suspect_pts = abs(filt_diff)>100; % seems that you have to input the threshold based on deployment

figure, %plot(v_wrapped(:, 1), '.')
hold on, plot(find(suspect_pts(:,1)),v_wrapped(suspect_pts(:,1)), 'r.', 'MarkerSize', 10)
hold on, plot(find(~suspect_pts(:,1)),v_wrapped(~suspect_pts(:, 1)), 'g.', 'MarkerSize', 10)

ylabel('Velocity, [m/s]', 'FontSize', 16)
xlabel('Point Index #', 'FontSize', 16)
lgd = legend('Velocity Wrapped Points', 'Non-Wrapped Points', 'fontsize', 16);
disp('Unwrapping...')
for ncol = 1:nt
prog = ncol/nt * 100;
fprintf('%.2f%% Complete\r', prog)
%ncol = 30;

si = find(suspect_pts(:, ncol));
%create difference operator
E = eye(nbins);
D = diff(E);
E = E(:,si);
v_prime = D*v_wrapped(:,ncol);

%solve least squares problem, and correct profile
r = round( (2*Vr(ncol)*D*E)\v_prime );
v_unwrap(si, ncol) = v_wrapped(si, ncol) - r*2*Vr(ncol);

end

fprintf('\n --- Unwrapped! --- \n')
%plots;
figure
ax1 = subplot(3,1,1);
pcolor((v_wrapped(:,1:ncol)./mean(Vr))');
shading flat
title('Original')
colorbar
caxis([-1,1])

ax2 =subplot(3,1,2);
pcolor((v_wrapped(:,1:ncol)./mean(Vr) - filt_diff1(:,1:ncol))');
shading flat
title('Filtered')
colorbar
caxis([-1,1])

ax3 = subplot(3,1,3);
pcolor((v_unwrap(:,1:ncol)./mean(Vr))');
shading flat
title('Unwrap')
colorbar
caxis([-1,1])
linkaxes([ax1 ax2 ax3],'x')


end


end


% EOF





