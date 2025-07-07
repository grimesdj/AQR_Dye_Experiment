
function A = L0_AQD(A)


% temporary addpath for testing :(
%addpath '/Users/jasonrooker/Library/CloudStorage/OneDrive-UNC-Wilmington/Kelp_repo/AQR_Dye_Experiment/Rooker_GettingStarted/code'
addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\Rooker_GettingStarted\code'

fprintf('\n============================\nDo you want to unwrap beam Velocities?')
unwrap = input('(1 = yes; 0 = no)');
if unwrap == 1
    for beam = 1:3
        A.L0.(sprintf('u%d',beam)) = aquawrap(A.(sprintf('Velocity_Beam%d',beam)), A.VRange);
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
A.dbins = (A.Config.blank + A.Config.binSize*A.Config.Nbins);
%
A.maxRange = (A.Pressure-A.pressureOffset).*cosd(20)-1*A.Config.binSize;
A.ylims      = [0 min(max(A.maxRange),max(A.dbins))];
dum1       = A.maxRange.*ones(1,A.Config.Nbins);
dum2       = ones(nsamples,1)*A.dbins;
qcFlag0    =  (dum2<=dum1);
A.qcFlag   =  double( (dum2<=dum1) & min(A.Amplitude_Beam1,min(A.Amplitude_Beam2,A.Amplitude_Beam3))>20 & min(A.Correlation_Beam1,min(A.Correlation_Beam2,A.Correlation_Beam3))>40 );
%
Time = datetime(A.Time,'convertFrom','datenum');
%
%
%for i = size(A.A1,1)
%   for j = size(A.A1, 2)
    %Amin = min(A.a1, min(A.a2, A.a3));
    %Cmin = min(A.c1, min(A.c2, A.c3));
%    end
%end

%A.Amplitude_Minimum = Amin;
%A.Correlation_Minimum = Cmin;

disp('Applying qcFlag to NaN data A < 20 and C < 40')
%
% Now plot currents

A.Velocity_East = A.Velocity_East.*A.qcFlag;
A.Velocity_North = A.Velocity_North.*A.qcFlag;
A.Velocity_Up = A.Velocity_Up.*A.qcFlag;

A.Velocity_East(~A.qcFlag')=nan;
A.Velocity_North(~A.qcFlag')=nan;
A.Velocity_Up(~A.qcFlag')=nan;
%
%
%
A.Velocity_Beam1 = A.Velocity_Beam1.*A.qcFlag;
A.Velocity_Beam2 = A.Velocity_Beam2.*A.qcFlag;
A.Velocity_Beam3 = A.Velocity_Beam3.*A.qcFlag;

A.Velocity_Beam1(~A.qcFlag')=nan;
A.Velocity_Beam2(~A.qcFlag')=nan;
A.Velocity_Beam3(~A.qcFlag')=nan;



%
A.Velocity_X = A.Velocity_X.*A.qcFlag;
A.Velocity_Y = A.Velocity_Y.*A.qcFlag;
A.Velocity_Z = A.Velocity_Z.*A.qcFlag;
%
%

%
% add the config info to the structure A to quick save as netcdf4
%fieldNames = fields(A.Config);
%originalFields = fields(A);
%
%for j = 1:length(fieldNames)
% A.(fieldNames{j}) = A.Config.(fieldNames{j});
%end
%A = orderfields(A,cat(1,fieldNames,originalFields));
%ncfile = [L0Dir,'/',L0Name,'.nc'];
%if exist(ncfile,'file')
    %eval(['!rm ',ncfile])
%end
disp('skipping nc file for now')
%struct2nc(A,ncfile,'NETCDF4');
%
%

% Make the names match convention
%A.Time = A.time;
%A.Config = A.config;
%A.Pressure = A.pressure;
%A.Bins = A.dbins;
fieldsToKeep = {'Time', 'Velocity_East', 'Velocity_North', 'Velocity_Up', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Velocity_Beam1', 'Velocity_Beam2', 'Velocity_Beam3', 'Amplitude_Minimum', 'Correlation_Minimum', 'Config', 'Pressure', 'Bins'};
A.L0 = rmfield(A, setdiff(fieldnames(A), fieldsToKeep));


end
