
function A = L0_AQD(A, atmTime, depTime)


%
%
% plot some stuff
if ~exist('atmTime','var')
    disp('pick 2 points bounding when out of water for ATM pressure offset')
    plot(A.pressure)
    l = ginput(2);
    l = round(l(:,1));
    atmTime = [A.time(l(1)), A.time(l(2))];
    fprintf('atmTime = \n')
    fprintf('%s --- %s', datestr(atmTime(1)), datestr(atmTime(2)));
else
    l = find(A.time>=atmTime(1) & A.time<=atmTime(2));
end
A.pressureOffset = mean(A.pressure(l(1):l(2)));
%
if ~exist('depTime','var')
    % now trim the data to when it was in the water
    disp('pick start/end points of deployment')
    l = ginput(2);
    l = round(l(:,1));
    depTime = [A.time(l(1)), A.time(l(2))];
    fprintf('depTime = \n')
    fprintf('%s --- %s', datestr(depTime(1)), datestr(depTime(2)));
else
    valid = find(A.time>=depTime(1) & A.time<=depTime(2));
end
%
vars  = {'time','volt','seconds','sspeed','heading','pitch','roll','pressure','temperature','a1','a2','a3','v1','v2','v3','c1','c2','c3','b1','b2','b3','east','north','up'};
for jj = 1:length(vars)
    eval(['A.',vars{jj},' = A.',vars{jj},'(valid,:);'])
end
nsamples = length(valid);
%
A.maxRange = (A.pressure-A.pressureOffset).*cosd(20)-1*A.config.binSize;
A.ylims      = [0 min(max(A.maxRange),max(A.dbins))];
dum1       = A.maxRange.*ones(1,A.config.Nbins);
dum2       = ones(nsamples,1)*A.dbins;
qcFlag0    =  (dum2<=dum1);
A.qcFlag   =  double( (dum2<=dum1) & min(A.a1,min(A.a2,A.a3))>30 & min(A.c1,min(A.c2,A.c3))>40 );
%
time = datetime(A.time,'convertFrom','datenum');
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

% quick convolution running mean filter
np1 = round(0.3/A.config.binSize);% 10 cm vertical 
np2 = 31;% 5min for 1Hz data
f1 = hamming(np1);f1 = f1./sum(f1);
f2 = hamming(np2);f2 = f2./sum(f2);
% do a nan-mean filter, keep track of normalization
on = conv2(f1,f2,A.qcFlag','same');
%
A.A1 = conv2(f1,f2,(A.a1.*A.qcFlag)','same')./on;
A.A2 = conv2(f1,f2,(A.a2.*A.qcFlag)','same')./on;
A.A3 = conv2(f1,f2,(A.a3.*A.qcFlag)','same')./on;
%
%
A.C1 = conv2(f1,f2,(A.c1.*A.qcFlag)','same')./on;
A.C2 = conv2(f1,f2,(A.c2.*A.qcFlag)','same')./on;
A.C3 = conv2(f1,f2,(A.c3.*A.qcFlag)','same')./on;
%
%for i = size(A.A1,1)
%   for j = size(A.A1, 2)
    Amin = min(A.A1, min(A.A2, A.A3));
    Cmin = min(A.C1, min(A.C2, A.C3));
%    end
%end

A.Amplitude_Minimum = Amin';
A.Correlation_Minimum = Cmin';

%
% Now plot currents
V1 = conv2(f1,f2,(A.v1.*A.qcFlag)','same')./on;
V2 = conv2(f1,f2,(A.v2.*A.qcFlag)','same')./on;
V3 = conv2(f1,f2,(A.v3.*A.qcFlag)','same')./on;

A.Velocity_East = (conv2(f1,f2,(A.east.*A.qcFlag)','same')./on)';
A.Velocity_North = (conv2(f1,f2,(A.north.*A.qcFlag)','same')./on)';
A.Velocity_Up = (conv2(f1,f2,(A.up.*A.qcFlag)','same')./on)';

A.Velocity_East(~A.qcFlag')=nan;
A.Velocity_North(~A.qcFlag')=nan;
A.Velocity_Up(~A.qcFlag')=nan;
%
%
%
A.Velocity_Beam1 = (conv2(f1,f2,(A.b1.*A.qcFlag)','same')./on)';
A.Velocity_Beam2 = (conv2(f1,f2,(A.b2.*A.qcFlag)','same')./on)';
A.Velocity_Beam3 = (conv2(f1,f2,(A.b3.*A.qcFlag)','same')./on)';

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
A.Pressure = A.pressure;
fieldsToKeep = {'Time', 'Velocity_East', 'Velocity_North', 'Velocity_Up', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Velocity_Beam1', 'Velocity_Beam2', 'Velocity_Beam3', 'Amplitude_Minimum', 'Correlation_Minimum', 'Config', 'Pressure'};
A.L0 = rmfield(A, setdiff(fieldnames(A), fieldsToKeep));
end
