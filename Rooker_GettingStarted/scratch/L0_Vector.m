
function A = L0_Vector(A, atmosphTime, deployTime, inputFile, figDir);

%
date = datetime(A.sensor.date,'convertFrom','datenum');
time = datetime(A.sensor.date(1)+A.seconds/86400,'convertFrom','datenum');
%
if ~exist('atmosphTime','var')
% plot some stuff
disp('pick 2 points bounding when out of water for ATM pressure offset')
plot(pressure)
l = ginput(2);
atmosphTime = [time(l(1,1)) time(l(2,1))]
end
if ~exist('deployTime','var')
% now trim the data to when it was in the water
disp('pick start/end points of deployment')
plot(pressure)
l = ginput(2);
deployTime = [time(l(1,1)) time(l(2,1))]
end
iATM = find(time>=atmosphTime(1) & time<=atmosphTime(2));
iDEP = find(time>=deployTime(1) & time<=deployTime(2));
%    
A.pressureOffset = mean(A.pressure(iATM));
%
vel_valid = iDEP;
vars  = {'seconds','pressure','a1','a2','a3','v1','v2','v3','c1','c2','c3','SNR1','SNR2','SNR3','b1','b2','b3','east','north','up'};
for jj = 1:length(vars)
    eval(['A.',vars{jj},' = A.',vars{jj},'(vel_valid);'])
end
nsamples = length(vel_valid);
time     = time(vel_valid);
A.time   = datenum(time);
%
sen_valid = find(date>=deployTime(1) &...
                  date<=deployTime(2));
vars  = {'date','battery_voltage','sound_speed','heading','pitch','roll','temperature'};
for jj = 1:length(vars)
    eval(['A.',vars{jj},' = A.sensor.',vars{jj},'(sen_valid);'])
end
A.seconds = A.seconds-A.seconds(1);
date = date(sen_valid);
%
A.qcFlag   =  min(A.a1,min(A.a2,A.a3))>60 & min(A.c1,min(A.c2,A.c3))>50 & min(A.SNR1,min(A.SNR2,A.SNR3))>15;
%
% make a few quick plots
fig0 = figure;
ax1 = subplot(2,1,1);
plot(date,A.temperature)
ylabel(ax1,'$T$ [$^\circ$]','interpreter','latex')
set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
ax2 = subplot(2,1,2);
plot(time,A.pressure)
ylabel(ax2,'$P$ [m]','interpreter','latex')
xlabel(ax2,'time','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out','xlim',get(ax1,'xlim')) 
figName = [figDir,inputFile,'_temp_pres.png'];
exportgraphics(fig0,figName)
%
%
% use acceleration and jolt to filter bad data
u1   = A.b1;
d1   = gradient(u1)/A.config.dt;
dd1  = gradient(d1)/A.config.dt;
u2   = A.b2;
d2   = gradient(u2)/A.config.dt;
dd2  = gradient(d2)/A.config.dt;
u3   = A.b3;
d3   = gradient(u3)/A.config.dt;
dd3  = gradient(d3)/A.config.dt;
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
A.A1 = (A.a1.*A.qcFlag)';
A.A2 = (A.a2.*A.qcFlag)';
A.A3 = (A.a3.*A.qcFlag)';
%
%
A.C1 = (A.c1.*A.qcFlag)';
A.C2 = (A.c2.*A.qcFlag)';
A.C3 = (A.c3.*A.qcFlag)';
%
for i = size(A.A1)
    Amin(i) = min(A.A1(i), min(A.A2(i), A.A3(i)));
    Cmin(i) = min(A.C1(i), min(A.C2(i), A.C3(i)));
end

A.Amplitude_Minimum = Amin';
A.Correlation_Minimum = Cmin';

%
% Now plot currents
V1 = (A.v1.*A.qcFlag)';
V2 = (A.v2.*A.qcFlag)';
V3 = (A.v3.*A.qcFlag)';

A.Velocity_East = (A.east.*A.qcFlag);
A.Velocity_North = (A.north.*A.qcFlag);
A.Velocity_Up = (A.up.*A.qcFlag);

A.Velocity_East(~A.qcFlag')=nan;
A.Velocity_North(~A.qcFlag')=nan;
A.Velocity_Up(~A.qcFlag')=nan;
%
%
%
A.Velocity_Beam1 = (A.b1.*A.qcFlag);
A.Velocity_Beam2 = (A.b2.*A.qcFlag);
A.Velocity_Beam3 = (A.b3.*A.qcFlag);

A.Velocity_Beam1(~A.qcFlag')=nan;
A.Velocity_Beam2(~A.qcFlag')=nan;
A.Velocity_Beam3(~A.qcFlag')=nan;

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
A.Velocity_X = V1';
A.Velocity_Y = V2';
A.Velocity_Z = V3';
%
A.VelXAvg = Uz;
A.VelYAvg = Vz;
A.VelZAvg = Wz;
%

%
% add the config info to the structure A to quick save as netcdf4
fieldNames = fields(A.config);
originalFields = fields(A);
%
for j = 1:length(fieldNames)
 A.(fieldNames{j}) = A.config.(fieldNames{j});
end
A = orderfields(A,cat(1,fieldNames,originalFields));
%ncfile = [L0Dir,'/',L0Name,'.nc'];
%if exist(ncfile,'file')
    %eval(['!rm ',ncfile])
%end
disp('skipping nc file')
%struct2nc(A,ncfile,'NETCDF4');
%
%

% Make the names match convention
A.Time = A.time;
A.Config = A.config;
fieldsToKeep = {'Time', 'Velocity_East', 'Velocity_North', 'Velocity_Up', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Velocity_Beam1', 'Velocity_Beam2', 'Velocity_Beam3', 'Amplitude_Minimum', 'Correlation_Minimum', 'Config'};
A.L0 = rmfield(A, setdiff(fieldnames(A), fieldsToKeep));
end