clear all
close all

%% User Input Data
addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\code'
addpath /Users/jasonrooker/Library/CloudStorage/OneDrive-UNC-Wilmington/Kelp_repo/AQR_Dye_Experiment/code


releasenum = 1; % Enter Release number here
release = string(releasenum);
unwrap = 1; % 1 for unwrap, 0 for not unwrap

atmTimes = [datenum('03-Jul-2024 14:00:00'), datenum('03-Jul-2024 18:00:00'); ...
            datenum('08-Jul-2024 16:00:00'), datenum('08-Jul-2024 16:30:00')];

depTimes = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00'); ...
            datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];


% Enter input /directory/ and fileName root without file extension
inputDir  = '../../../../Kelp_data/data/Release' + release + '/raw/';
inputFile = 'KELP1_AquadoppHR';
fileName  = [inputDir + inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '../../../../Kelp_data/Summer2025/Rooker/Release' + release + '/raw/';
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = '../../../../Kelp_data/Summer2025/Rooker/Release' + release + '/L0/';
L0Name  = [inputFile,'_L0'];
% Enter time when instrument was in air for pressure offset
atmTime = atmTimes(releasenum, :);
depTime = depTimes(releasenum, :);
% Enter path to save figures
figDir = [inputDir + '/../figures/'];
if ~exist(figDir,'dir'), mkdir(figDir), end
%
% Enter time-offset (UTC->EDT) tos = -4 hours
tos = 0;
%
% returns structure A with all aquadopp data
fprintf('loading AQD release %s...\n', release)
A = loadAQD(inputDir, inputFile, tos, fileName);
%disp('saving disabled')
save([outputDir,'/',outputName,'.mat'],'-struct','A')
%
%
% plot some stuff
fprintf('Barometric Compensation...')
if ~exist('atmTime','var')
    disp('pick 2 points bounding when out of water for ATM pressure offset')
    plot(A.Pressure)
    l = ginput(2);
    l = round(l(:,1));
    atmTime = [A.Time(l(1)), A.Time(l(2))];
    fprintf('atmTime = \n')
    fprintf('%s --- %s', datestr(atmTime(1)), datestr(atmTime(2)));
else
    l = find(A.Time>=atmTime(1) & A.Time<=atmTime(2));
end
A.PressureOffset = mean(A.Pressure(l(1):l(2)));
%
if ~exist('depTime','var')
    % now trim the data to when it was in the water
    disp('pick start/end points of deployment')
    l = ginput(2);
    l = round(l(:,1));
    depTime = [A.Time(l(1)), A.Time(l(2))];
    fprintf('depTime = \n')
    fprintf('%s --- %s', datestr(depTime(1)), datestr(depTime(2)));
else
    valid = find(A.Time>=depTime(1) & A.Time<=depTime(2));
end
%
vars  = {'Time','volt','seconds','sspeed','heading','pitch','roll','Pressure','Temperature','Amplitude_Beam1','Amplitude_Beam2','Amplitude_Beam3','Velocity_X','Velocity_Y','Velocity_Z','Correlation_Beam1','Correlation_Beam2','Correlation_Beam3','Velocity_Beam1','Velocity_Beam2','Velocity_Beam3','Velocity_East','Velocity_North','Velocity_Up', 'VRange'};
for jj = 1:length(vars)
    var = vars{jj};
    fprintf('storing valid vars: %s\r', var)
    dum = A.(var);
    A.(var) = dum(valid,:);
    %eval(['A.',vars{jj},' = A.',vars{jj},'(valid,:);'])
end
nsamples = length(valid);
%
A.maxRange = (A.Pressure-A.PressureOffset).*cosd(20)-1*A.Config.binSize;
ylims      = [0 min(max(A.maxRange),max(A.dbins))];
dum1       = A.maxRange.*ones(1,A.Config.Nbins);
dum2       = ones(nsamples,1)*A.dbins;
qcFlag0    =  (dum2<=dum1);
A.qcFlag   =  double( (dum2<=dum1) & min(A.Amplitude_Beam1,min(A.Amplitude_Beam2,A.Amplitude_Beam3))>75 & min(A.Correlation_Beam1,min(A.Correlation_Beam2,A.Correlation_Beam3))>30 );
%
time = datetime(A.Time,'convertFrom','datenum');

%% Fix Velocity Wrapping

if ~exist('unwrap')
    fprintf('\n============================\nDo you want to unwrap beam Velocities?')
    unwrap = input('(1 = yes; 0 = no)');
end
if unwrap == 1
    for beam = 1:3 
        fprintf('unwrapping beam %d\n', beam)
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

%% QA/QC
%
% use acceleration and jolt to filter bad data
u1   = A.Velocity_Beam1;
d1   = gradient(u1)/A.Config.dt;
dd1  = gradient(d1)/A.Config.dt;
u2   = A.Velocity_Beam2;
d2   = gradient(u2)/A.Config.dt;
dd2  = gradient(d2)/A.Config.dt;
u3   = A.Velocity_Beam3;
d3   = gradient(u3)/A.Config.dt;
dd3  = gradient(d3)/A.Config.dt;
%
r01  =  nanstd(u1(:));
r02  =  nanstd(u2(:));
r03  =  nanstd(u3(:));
R0   = (u1./r01).^2 + (u2./r02).^2 + (u3./r03).^2;
%
r11  = nanstd(d1(:));
r12  = nanstd(d2(:));
r13  = nanstd(d3(:));
R1   = (d1./r11).^2 + (d2./r12).^2 + (d3./r13).^2;
%
r21  = nanstd(dd1(:));
r22  = nanstd(dd2(:));
r23  = nanstd(dd3(:));
R2   = (dd1./r21).^2 + (dd2./r22).^2 + (dd3./r23).^2;
%
valid = (R0<0.5);% & (R1<5) & (R2<10);
A.qcFlag = A.qcFlag & valid;
%
% make a few quick plots
fig0 = figure;
ax1 = subplot(2,1,1);
plot(time,A.Temperature)
ylabel(ax1,'$T$ [$^\circ$]','interpreter','latex')
set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
ax2 = subplot(2,1,2);
plot(time,A.Pressure)
ylabel(ax2,'$P$ [m]','interpreter','latex')
xlabel(ax2,'time [s]','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out')
figName = [figDir + filesep + inputFile + '_temperature_pressure.png'];
exportgraphics(fig0,figName)
%
%
% quick convolution running mean filter
np1 = round(0.3/A.Config.binSize);% 10 cm vertical 
np2 = 31;% 5min for 1Hz data
f1 = hamming(np1);f1 = f1./sum(f1);
f2 = hamming(np2);f2 = f2./sum(f2);
% do a nan-mean filter, keep track of normalization
on = conv2(f1,f2,A.qcFlag','same');
%
A1 = conv2(f1,f2,(A.Amplitude_Beam1.*A.qcFlag)','same')./on;
A2 = conv2(f1,f2,(A.Amplitude_Beam2.*A.qcFlag)','same')./on;
A3 = conv2(f1,f2,(A.Amplitude_Beam3.*A.qcFlag)','same')./on;
%
fig1 = figure;
ax1 = subplot(3,1,1);
imagesc(time,A.dbins',A1.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'X')
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
title(ax1,'Amplitude')
%
ax2 = subplot(3,1,2);
imagesc(time,A.dbins',A2.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Y')
ylabel('mab','interpreter','latex')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
%
ax3 = subplot(3,1,3);
imagesc(time,A.dbins',A3.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir + filesep + inputFile + '_amplitude.png'];
exportgraphics(fig1,figName)
%
%
C1 = conv2(f1,f2,(A.Correlation_Beam1.*A.qcFlag)','same')./on;
C2 = conv2(f1,f2,(A.Correlation_Beam2.*A.qcFlag)','same')./on;
C3 = conv2(f1,f2,(A.Correlation_Beam3.*A.qcFlag)','same')./on;
%
fig1 = figure;
ax1 = subplot(3,1,1);
imagesc(time,A.dbins',C1.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'X')
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
title(ax1,'Correlation')
%
ax2 = subplot(3,1,2);
imagesc(time,A.dbins',C2.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Y')
ylabel('mab','interpreter','latex')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
%
ax3 = subplot(3,1,3);
imagesc(time,A.dbins',C3.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir + filesep + inputFile + '_correlation.png'];
exportgraphics(fig1,figName)
%
%
% Now plot currents
V1 = conv2(f1,f2,(A.Velocity_X.*A.qcFlag)','same')./on;
V2 = conv2(f1,f2,(A.Velocity_Y.*A.qcFlag)','same')./on;
V3 = conv2(f1,f2,(A.Velocity_Z.*A.qcFlag)','same')./on;

%
fig2 = figure;
ax1 = subplot(3,1,1);
imagesc(time,A.dbins',V1.*qcFlag0'),caxis([-0.25 0.25]),colormap(cmocean('balance')),colorbar
text(time(1),ylims(2),'X')
%
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
ax2 = subplot(3,1,2);
imagesc(time,A.dbins',V2.*qcFlag0'),caxis([-0.25 0.25]),colormap(cmocean('balance')),colorbar
text(time(1),ylims(2),'Y')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
ylabel('mab','interpreter','latex')
%
ax3 = subplot(3,1,3);
imagesc(time,A.dbins',V3.*qcFlag0'),caxis([-0.125 0.125]),colormap(cmocean('balance')),colorbar
text(time(1),ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir + filesep + inputFile + '_velocity.png'];
exportgraphics(fig2,figName)
%
% get the time-averaged current.
% Note, this is not in depth normalized (sigma) coordinates.
% first, nan any values that don't pass QC.
V1(~A.qcFlag')=nan;
V2(~A.qcFlag')=nan;
V3(~A.qcFlag')=nan;
flag = sum(A.qcFlag,1) > 0.50*nsamples;
% Uz = time averegd; U = depth & time averaged
Uz = nanmean(V1,2); U = nanmean(Uz); Uz(~flag)=nan;
Vz = nanmean(V2,2); V = nanmean(Vz); Vz(~flag)=nan;
Wz = nanmean(V3,2); W = nanmean(Wz); Wz(~flag)=nan;
%
fig3 = figure;
plot(Uz,A.dbins,'.-r',Vz,A.dbins,'.-b',Wz,A.dbins,'.-k')
hh = legend('X','Y','Z');
set(hh,'fontsize',9)
xlabel('m/s','interpreter','latex')
ylabel('mab','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
figName = [figDir + filesep + inputFile + '_mean_velocity.png'];
exportgraphics(fig3,figName)
%
A.VelX = V1;
A.VelY= V2;
A.VelZ   = V3;
%
A.VelXAvg = Uz;
A.VelYAvg= Vz;
A.VelZAvg   = Wz;

%
%% QC velocities
A.Velocity_X(~A.qcFlag') = nan;
A.Velocity_Y(~A.qcFlag') = nan;
A.Velocity_Z(~A.qcFlag') = nan;

A.Velocity_Beam1(~A.qcFlag') = nan;
A.Velocity_Beam2(~A.qcFlag') = nan;
A.Velocity_Beam3(~A.qcFlag') = nan;

A.Velocity_East(~A.qcFlag') = nan;
A.Velocity_North(~A.qcFlag') = nan;
A.Velocity_Up(~A.qcFlag') = nan;


save(fullfile(L0Dir, L0Name, '.mat'),'-struct','A')
%
% add the config info to the structure A to quick save as netcdf4
fieldNames = fields(A.Config);
originalFields = fields(A);
%
for j = 1:length(fieldNames)
 A.(fieldNames{j}) = A.Config.(fieldNames{j});
end
A = orderfields(A,cat(1,fieldNames,originalFields));
ncfile = [L0Dir + L0Name + '.nc'];
if exist(ncfile,'file')
    delete(ncfile)
end
struct2nc(A,ncfile,'NETCDF4');
%
%
v1bar = nansum( A.Velocity_X.*A.qcFlag ,2)./sum(A.qcFlag,2);
v2bar = nansum( A.Velocity_Y.*A.qcFlag ,2)./sum(A.qcFlag,2);
v3bar = nansum( A.Velocity_Z.*A.qcFlag ,2)./sum(A.qcFlag,2);
% rotate 60 degrees:
vrot = [cosd(60) sind(60);-sind(60) cosd(60)]*[v1bar';v2bar'];
v1rot = vrot(1,:)';
v2rot = vrot(2,:)';
%
%
dv    = 0.01;
vbins = [-1:dv:1];
v1H   = hist(v1rot,vbins);
v2H   = hist(v2rot,vbins);
v3H   = hist(v3bar,vbins);
%
figure, plot(v1rot,v2rot,'.')
hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
xlabel('$\bar{v}_X$','interpreter','latex')
ylabel('$\bar{v}_Y$','interpreter','latex')
figName = [figDir + filesep + inputFile + '_depth_avg_XYvel_checker_pattern.png'];
exportgraphics(gcf,figName)
%
%
figure,
subplot(3,1,1)
bar(vbins, v1H./length(v1bar))
ylabel('$p(\bar{v}_{X''})$','interpreter','latex')
subplot(3,1,2)
bar(vbins, v2H./length(v1bar))
ylabel('$p(\bar{v}''_{Y''})$','interpreter','latex')
subplot(3,1,3)
bar(vbins, v3H./length(v1bar))
ylabel('$p(\bar{v}_Z)$','interpreter','latex')
xlabel('velocity')
figName = [figDir + filesep + inputFile + '_depth_avg_vel_XYrot60deg_histograms.png'];
exportgraphics(gcf,figName)
%
%
%
%
b1bar = nansum( A.Velocity_Beam1.*A.qcFlag ,2)./sum(A.qcFlag,2);
b2bar = nansum( A.Velocity_Beam2.*A.qcFlag ,2)./sum(A.qcFlag,2);
b3bar = nansum( A.Velocity_Beam3.*A.qcFlag ,2)./sum(A.qcFlag,2);
figure, plot(b1bar,b2bar,'.')
hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
xlabel('$\bar{v}_X$','interpreter','latex')
ylabel('$\bar{v}_Y$','interpreter','latex')
figName = [figDir + filesep + inputFile + '_depth_avg_Beam1Beam2_checker_pattern.png'];
exportgraphics(gcf,figName)

%% Functions

function A = loadAQD(inputDir, inputFile, tos, fileName)

%% load header data
%% data for each field starts at column 39 or 40
hdrFile = sprintf(['%s',filesep,'%s.hdr'], inputDir,inputFile);
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
%% Extracts date, temp, pressure, heading, pitch, roll
senFile = sprintf(['%s',filesep,'%s.sen'], inputDir,inputFile);
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
A.Config= meta_data;
A.date  = datestr(time(1));
A.Time  = time; A.volt = volt;
A.seconds= (time-time(1))*86400;
A.sspeed= sspeed;
A.VRange = vrange;
A.heading  = heading;A.pitch = pitch;A.roll = roll;
A.Pressure  = pressure;A.Temperature = temperature;
A.fname = fileName;
%
%% load beam amplitudes and velocities (may be in beam coords or ENU, see A.Config)
a1 = load(strcat(fileName,'.a1'));
[Na,Ma] = size(a1);
if Ma>nbins% sometimes there are extra columns (beam#, ensemble#)
    bin1 = Ma-nbins+1;
else
    bin1 = 1;
end
a2 = load(strcat(fileName,'.a2'));
a3 = load(strcat(fileName,'.a3'));
A.Amplitude_Beam1 = a1(:,bin1:end);
A.Amplitude_Beam2 = a2(:,bin1:end);
A.Amplitude_Beam3 = a3(:,bin1:end);
%
v1 = load(strcat(fileName,'.v1'));
v2 = load(strcat(fileName,'.v2'));
v3 = load(strcat(fileName,'.v3'));
switch coords
  case {'XYZ','ENU'}
    A.Velocity_X = v1(:,bin1:end);
    A.Velocity_Y = v2(:,bin1:end);
    A.Velocity_Z = v3(:,bin1:end);
    if strcmp(coords,'XYZ')
        shape = size(A.Velocity_X);
        BEAM  = inv(T)*[A.Velocity_X(:)'; A.Velocity_Y(:)'; A.Velocity_Z(:)'];
        A.Velocity_Beam1  = reshape(BEAM(1,:)',shape);
        A.Velocity_Beam2  = reshape(BEAM(2,:)',shape);
        A.Velocity_Beam3  = reshape(BEAM(3,:)',shape);
    end        
    A.Velocity_East = A.Velocity_X;
    A.Velocity_North= A.Velocity_Y;
    A.Velocity_Up   = A.Velocity_Z;
  case {'BEA'}
    b1 = v1(:,bin1:end);
    b2 = v2(:,bin1:end);
    b3 = v3(:,bin1:end);
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
    hh = reshape(pi*(heading-90)/180,1,1,nsamples);
    pp = reshape(pi*pitch/180,1,1,nsamples);
    rr = reshape(pi*roll/180,1,1,nsamples);
    H = [ cos(hh) sin(hh) 0*hh;...
         -sin(hh) cos(hh) 0*hh;...
          0*hh      0*hh  0*hh+1];
    P = [cos(pp) -sin(pp).*sin(rr) -cos(rr).*sin(pp);...
          0*pp       cos(rr)         -sin(rr)   ;...
         sin(pp)  sin(rr).*cos(pp)  cos(pp).*cos(rr)];
    shape = size(A.Velocity_X);
    for j = 1:nsamples
     R   = H(:,:,j)*P(:,:,j);
     ENU = R*[A.Velocity_X(j,:);A.Velocity_Y(j,:);A.Velocity_Z(j,:)];
     A.Velocity_East (j,:) = ENU(1,:);
     A.Velocity_North(j,:) = ENU(2,:);
     A.Velocity_Up   (j,:) = ENU(3,:);    
    end
end
%
c1 = load(strcat(fileName,'.c1'));
c2 = load(strcat(fileName,'.c2'));
c3 = load(strcat(fileName,'.c3'));
A.Correlation_Beam1 = c1(:,bin1:end);
A.Correlation_Beam2 = c2(:,bin1:end);
A.Correlation_Beam3 = c3(:,bin1:end);
%
A.dbins = blank + binsize*[1:nbins];
% 
end

%EOF 


function [v_unwrap, suspect_pts] = unwrap_AQD(v_wrapped, Vr)
%   unwrap_AQD
%       Corrects ADCP velocity-wrapped data by idententifying suprious
%       points using a median filter, estimating the error using a least-sqraures-regression, 
%       and correcting by shifting the
%       appropriate integer multiples of the maximum velocity range
%
%        -- Originally intended for AquaDopp but may be suitable for other
%        pulse-coherent ADCPs -- 
%
%   Written by:
%       Jason Rooker, Summer 2025
%       Derek Grimes
%       Based on the method originally devoloped by
%                   Andrey Shcherbina, Eric D'Asaro, and Sven Nylund (2018)
%       
%
%       USAGE: [v_unwrap, suspect_pts] = unwrap_AQD(v_wrapped, Vr) 
%
%       Inputs:
%               v_wrapped:  vector of wrapped velocity measurements
%               Vr:         Float, Maximum Velocity Range
%
%       Outputs:
%               v_unwrap:     vector of unwrapped velocities
%               suspect_pts:  indicies of points that were determined to be
%                             wrapped

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
