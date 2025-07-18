
function L0_Vector(outputFile, L0Dir, L0Name)

A = load([outputFile, '.mat']);
% temporary addpath for testing :(
addpath '/Users/jasonrooker/Library/CloudStorage/OneDrive-UNC-Wilmington/Kelp_repo/AQR_Dye_Experiment/Rooker_GettingStarted/code'
%addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\Rooker_GettingStarted\code'

fprintf('\n============================\n\nDo you want to unwrap beam Velocities?')
unwrap = input('\n(1 = yes; 0 = no)');
if unwrap == 1
    for beam = 1:3
        Vel = (sprintf('Velocity_Beam%d',beam));
        vwrap = reshape(A.(Vel), [], 24);
        Vr = reshape(A.VRange, [], 24);
        Vr = Vr.*0.01;
        [unwrap, A.(sprintf('Suspect_Beam%d', beam))] = aquawrap(vwrap, Vr);
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
    HPRfiles = dir([L0Dir, '/../raw/*AquadoppHR_raw.mat']);
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
% fieldsToKeep = {'Time', 'Velocity_East', 'Velocity_North', 'Velocity_Up', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Velocity_Beam1', 'Velocity_Beam2', 'Velocity_Beam3', 'Amplitude_Minimum', 'Correlation_Minimum', 'Config', 'Pressure','Heading', 'Pitch', 'Roll', 'Bins', 'Temperature', 'u1', 'u2', 'u3'};
% L0 = rmfield(A, setdiff(fieldnames(A), fieldsToKeep));

disp('Saving L0 data')
save([L0Dir,'/',L0Name,'.mat'],'-struct','A')

end
