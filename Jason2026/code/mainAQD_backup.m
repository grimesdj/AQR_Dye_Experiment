% mainAQD.m
% 
%   USAGE: loads AQD data from textfiles into .mat files
% 
%   (script) -> requires user to enter:
%       inputDir = Directory where textfiles are
%       inputFile = file root for AQD files
%       outputDir = Directory to save raw .mat
%       L0Dir = Directory to save L0 .mat
%       atmTime = two datetimes when instrument was in the air
%       depTime = start and end times of deployment
%       


clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '../../../../Kelp_data/data/Release1/raw';
inputFile = 'KELP1_AquadoppHR';
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '../../../../Kelp_data/Summer2025/Rooker/Release1/raw';
outputName= [inputFile,'_raw'];
outputFile= [outputDir, '/', outputName];
% Enter processed output fileName without .mat
L0Dir   = '../../../../Kelp_data/Summer2025/Rooker/Release1/L0';
L0Name  = [inputFile,'_L0'];
% Enter time when instrument was in air for pressure offset
atmTime = [datenum('03-Jul-2024 14:00:00'), datenum('03-Jul-2024 18:00:00')];
depTime = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00')];
% Enter path to save figures
figDir = [outputDir,'/../figures/'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-offset (UTC->EDT) tos = -4 hours
tos = 0;
%
% Generate and save raw data
disp('Generating raw data')
loadAQD(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime);

% Generate and save L0
disp('Generating L0 data')
L0_AQD(outputFile, L0Dir, L0Name);

% % Create Plots
% disp('Generating figures')
% L0_plots(L0Dir, L0Name, figDir, inputFile)


fprintf('\n====================\n\nDone!\n\n====================\n')

%% Functions

function loadAQD(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime)
% 
%   USAGE: loadAQD(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime)
%       inputDir  = Directory where textfiles are
%       inputFile = file root for AQD files
%       fileName  = Directory and file root
%       outputDir = Directory to save raw .mat
%       atmTime   = [Datenum Datenum] in air
%       depTime   = [Start_time , End_time] in water
%       tos       = time offset for time zone correction


% pull data

HRflag = 0;

% load header data
% data for each field starts at column 39 or 40
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
% Extracts date, temp, pressure, Vr, heading, pitch, roll
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

% load beam amplitudes and velocities (may be in beam coords or ENU, see A.config)

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
    A.Correlation_Beam1 = NaN;
    A.Correlation_Beam2 = NaN;
    A.Correlation_Beam3 = NaN;
end

%
disp('Saving raw data')
save([outputFile,'.mat'],'-struct','A')


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
ylim([-1 1])

ax2 = subplot(3, 1, 2);
plot(ax2, A.Time, A.Velocity_Beam2, '.')
ylabel({'Beam 2', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
hold on 
plot(ax2, A.Time, A.VRange, 'r', A.Time, -1*A.VRange, 'r');
grid minor
ylim([-1 1])

ax3 = subplot(3, 1, 3);
plot(ax3, A.Time, A.Velocity_Beam3, '.')
ylabel({'Beam 3', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
hold on 
plot(ax3, A.Time, A.VRange, 'r', A.Time, -1*A.VRange, 'r');
grid minor
ylim([-1 1])

sgtitle('AQD Raw Beam Velocities', 'Fontsize', 25)

end

% EOF
function L0_AQD(outputFile, L0Dir, L0Name)
% 
%   USAGE: L0_AQD(outputFile, L0Dir, L0Name)
%       outputFile = folder path and filename (no extension) to raw data
%       L0Dir = directory to save finished L0 data
%       L0Name = L0 filename (without extension) to be saved
% 
%       takes raw AQD data and performs L0 QA/QC
% 
A = load([outputFile, '.mat']);
fprintf('\n============================\nDo you want to unwrap beam Velocities?')
unwrap = input('(1 = yes; 0 = no)');
if unwrap == 1
    for beam = 1:3 
        [A.(sprintf('Velocity_Beam%d',beam)), A.(sprintf('Suspect_Beam%d', beam))] = unwrap_AQD(A.(sprintf('Velocity_Beam%d',beam)), A.VRange);
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
    disp('rotating unwrapped data to ENU')
    hh = reshape(pi*(A.Heading-90)/180,1,1,A.Config.Nsamples);
    pp = reshape(pi*A.Pitch/180,1,1,A.Config.Nsamples);
    rr = reshape(pi*A.Roll/180,1,1,A.Config.Nsamples);
    H = [ cos(hh), sin(hh), 0*hh;...
         -sin(hh), cos(hh), 0*hh;...
          0*hh,      0*hh,  0*hh+1];
    P = [cos(pp), 0*pp, -sin(pp);...
          0*pp,  0*pp+1, 0*pp   ;...
         sin(pp),  0*pp,  cos(pp)];
     
    O = [1+0*rr 0*rr 0*rr;...
        0*rr cos(rr) -sin(rr);...
        0*rr sin(rr) cos(rr)];
     
     
    
    for j = 1:A.Config.Nsamples
        R   = H(:,:,j)*P(:,:,j)*O(:,:,j);
        ENU = R*[A.Velocity_X(j,:);A.Velocity_Y(j,:);A.Velocity_Z(j,:)];
        A.Velocity_East  (j,:) = ENU(1,:);
        A.Velocity_North (j,:) = ENU(2,:);
        A.Velocity_Up    (j,:) = ENU(3,:);
    end

end

%

A.dbins = (A.Config.blank + A.Config.binSize*A.Config.Nbins);
A.maxRange = (A.Pressure-A.pressureOffset).*cosd(20)-1*A.Config.binSize;
A.ylims      = [0 min(max(A.maxRange),max(A.dbins))];
dum1       = A.maxRange.*ones(1,A.Config.Nbins);
dum2       = ones(A.Config.Nsamples,1)*A.dbins;
qcFlag0    =  (dum2<=dum1);
A.Amplitude_Minimum   = min(A.Amplitude_Beam1, min(A.Amplitude_Beam2, A.Amplitude_Beam3));
A.Correlation_Minimum = min(A.Correlation_Beam1, min(A.Correlation_Beam2, A.Correlation_Beam3, 'omitnan'), 'omitnan');   
A.qcFlag              =  double( qcFlag0 & A.Amplitude_Minimum > 20 & A.Correlation_Minimum > 40 );

% QCFlag
disp('Applying qcFlag to NaN data A < 20 and C < 40')

%

A.Velocity_East = A.Velocity_East.*A.qcFlag;
A.Velocity_North = A.Velocity_North.*A.qcFlag;
A.Velocity_Up = A.Velocity_Up.*A.qcFlag;

A.Velocity_East(~A.qcFlag)=nan;
A.Velocity_North(~A.qcFlag)=nan;
A.Velocity_Up(~A.qcFlag)=nan;

A.Velocity_Beam1 = A.Velocity_Beam1.*A.qcFlag;
A.Velocity_Beam2 = A.Velocity_Beam2.*A.qcFlag;
A.Velocity_Beam3 = A.Velocity_Beam3.*A.qcFlag;

A.Velocity_Beam1(~A.qcFlag)=nan;
A.Velocity_Beam2(~A.qcFlag)=nan;
A.Velocity_Beam3(~A.qcFlag)=nan;

A.Velocity_X = A.Velocity_X.*A.qcFlag;
A.Velocity_Y = A.Velocity_Y.*A.qcFlag;
A.Velocity_Z = A.Velocity_Z.*A.qcFlag;

A.Velocity_X(~A.qcFlag)=nan;
A.Velocity_Y(~A.qcFlag)=nan;
A.Velocity_Z(~A.qcFlag)=nan;

%

disp('skipping nc file for now')

%

% Summary figure
figure
ax1 = subplot(3, 1, 1);
plot(ax1, A.Time, A.Velocity_Beam1, '.')
ylabel({'Beam 1', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor
ylim([-1 1])

ax2 = subplot(3, 1, 2);
plot(ax2, A.Time, A.Velocity_Beam2, '.')
ylabel({'Beam 2', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor
ylim([-1 1])

ax3 = subplot(3, 1, 3);
plot(ax3, A.Time, A.Velocity_Beam3, '.')
ylabel({'Beam 3', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
ylim([-1 1])

sgtitle('AQD L0 Beam Velocities', 'Fontsize', 25)

disp('Saving L0 data')
save([L0Dir,'/',L0Name,'.mat'],'-struct','A')

end

% EOF

function [v_unwrap, suspect_pts] = unwrap_AQD(v_wrapped, Vr)

disp('using Shcherbina et al 2018 unwrapper')

[nbins, nt] = size(v_wrapped);
v_unwrap = v_wrapped;
filt_diff = (v_wrapped - medfilt1(v_wrapped,150))./Vr;
filt_diff1 = filt_diff;
suspect_pts = abs(filt_diff)>1; % seems that you have to input the threshold based on deployment

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

% EOF