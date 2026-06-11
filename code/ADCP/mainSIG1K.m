clear all
close all

%% User Input Data

% addpath for struct2nc because it fails as an internal function?
addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\code'


releasenum = 1; % Enter Release number here
releasenum = string(releasenum);

% stages of processing
% 1) define deployment number:
adcp_ID  = 1;
adcp_file_roots = {'S103071A014_KELP1','S104339A001_KELP1'};
adcp_mooring_ID = {'M2', 'M1'};
echo_mode = 0;
% used to shift timezone,e.g. to convert EST to UTC
time_shift = 0/24;
% 2) raw data input directory & filename convention:
rootDIR = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', 'FullExperiment', 'raw', sprintf('%s',adcp_file_roots{adcp_ID}));
fRoot   = [adcp_file_roots{adcp_ID},'_'];
% 3) output directory:
outRoot = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA');
% 4) output data file prefix:
filePrefix= sprintf('ADCP_%s_',adcp_mooring_ID{adcp_ID});
%
% 4a) current/echo average interval (seconds)
dtAvg     = 600;
L0dir     = [outRoot, filesep, adcp_mooring_ID{adcp_ID},filesep,'L0',filesep,'ADCP',filesep];
L0FRoot   = sprintf('%sL0_%dmin',filePrefix,dtAvg/60);
L1dir     = [outRoot, filesep, adcp_mooring_ID{adcp_ID},filesep,'L1',filesep,'ADCP',filesep];
L1FRoot   = sprintf('%sL1',filePrefix);
%
% $$$ atmosphTime = [datenum('13-May-2024 11:30:54'), datenum('13-May-2024 15:46:17');...
% $$$                datenum('03-Jun-2024 16:39:11'), datenum('04-Jun-2024 02:41:20')];%[datenum('15-Feb-2024 13:12:05'), datenum('15-Feb-2024 14:08:32')];
% $$$ 
% $$$ deployTime  = [datenum('13-May-2024 16:00:00')];%[datenum('15-Feb-2024 15:00:00')]; %datenum('09-Oct-2023 16:00:00');
% $$$ recoverTime = [datenum('03-Jun-2024 16:00:00')];%[datenum('17-Mar-2024 15:00:00')]; %datenum('30-Oct-2023 14:00:00');
%
files = dir([rootDIR,filesep,fRoot,'*.mat']);
Nf    = length(files);
%
% height of transducer off of the bottom;
hab = 0.25;% need to measure this!
HeadingOffset = 0;
%
% 5) time-periods when instrument was air (leave times empty to manually reselect them)
% 6) deploy/recovery times
switch adcp_ID
  case 1
    offset_file = fullfile(rootDIR, sprintf('%s%d.mat', fRoot, 2));
    load(offset_file,'Config','Data','Descriptions');
    atmosphTime = [datenum('02-Jul-2024 18:45:00'), datenum('02-Jul-2024 19:15:00')];
    it          = find(Data.Burst_Time>=atmosphTime(1) & Data.Burst_Time<=atmosphTime(2));
    ATM_Time    = nanmean(atmosphTime);
    ATM_Pressure = nanmean(Data.Burst_Pressure(it));
  case 2
    offset_file = fullfile(rootDIR, sprintf('%s%d.mat', fRoot, 1));
    load(offset_file,'Config','Data','Descriptions');
    atmosphTime = [datenum('02-Jul-2024 19:15:00'), datenum('02-Jul-2024 19:45:00')];
    it          = find(Data.Burst_Time>=atmosphTime(1) & Data.Burst_Time<=atmosphTime(2));
    ATM_Time    = nanmean(atmosphTime);
    ATM_Pressure = nanmean(Data.Burst_Pressure(it));
end
%
deployTime  =  datenum('03-Jul-2024 00:00:00');
recoverTime  = datenum('25-Jul-2024 00:00:00');
%
% shift time limits
atmosphTime = atmosphTime + time_shift;
deployTime  = deployTime  + time_shift;
recoverTime = recoverTime + time_shift;
%
% load and pre-process data.
fprintf('loading %s data\n', adcp_mooring_ID{adcp_ID})
loadSIG1K(rootDIR, fRoot, ATM_Time, ATM_Pressure, Config, hab, L0dir, filePrefix, deployTime, recoverTime, echo_mode, Descriptions, HeadingOffset)
%
% make time-averages
fprintf('time averaging %s\n', adcp_mooring_ID{adcp_ID})
time_average_and_rotate_sig1000_RDI_matrix_format_function(Config, L0dir, filePrefix, dtAvg, echo_mode, L0FRoot)
%
% estimate hourly wave stats
fprintf('estimating wave stats for %s\n', adcp_mooring_ID{adcp_ID})
waves = estimate_wave_bulk_stats_SIG1000_RDI_matrix_format_function(Config, L0dir, filePrefix, L0FRoot, hab, L1dir, L1FRoot);

%
%
fig = figure;
imagesc(datetime(waves.Time','convertFrom','datenum'),waves.frequency,log10(waves.Spp)),colormap(cmocean('thermal')),caxis([-2 1])
ylabel('$f$ [Hz]','interpreter','latex')
cb = colorbar;
ylabel(cb,'[m$^2$/Hz]','interpreter','latex')
ax = gca;
set(ax,'ydir','normal','ticklabelinterpreter','latex','tickdir','out','plotboxaspectratio',[1.5 1 1])
figname = sprintf('%s/figures/%s_spectra.pdf',L1dir,L1FRoot);
exportgraphics(fig,figname)
% 





% %% Release specific times


fprintf('compiling release times\n')
depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00') ;
            datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];

for r = 1:size(depTime, 1)
    fprintf('fetching release %d\n', r)
    R = fetch_sig1k(L0dir, depTime(r, 1), depTime(r, 2));

    cfgFile = dir(fullfile(L0dir, '*config.mat'));

    if isempty(cfgFile)
        error('No config file found in %s', L0dir);
    end
    
    cFig = load(fullfile(cfgFile(1).folder, cfgFile(1).name));
    R.Config = cFig;
    fprintf('saving %s...\n',  "/Release" + r + '/L0/ADCP/', adcp_mooring_ID(adcp_ID) + "_ADCP.mat")
    save(fullfile(outRoot, '..', "Release" + r, 'L0', 'ADCP', adcp_mooring_ID(adcp_ID) + "_ADCP.mat"), '-struct', 'R')
end



% L1 DATA IS SMOOTHED 
% 
% 
% Moor = load(fullfile(L1dir, L1FRoot));
% data = Moor.currents;
% fields = fieldnames(data);
% for r = 1:2
%     dep(r, :) = find(data.Time >= depTime(r,1) & data.Time <= depTime(r,2));
%     for i = 1:length(fields)
%         field = fields{i};
%         len = size(data.(field), 2);
%         if len == length(data.Time)
%             dum = data.(field);
%             dum = dum(:, dep(r,:));
%             trim.(field) = dum;
%         else
%             trim.(field) = data.(field);
%         end
%     end
%     trimDir = fullfile(outRoot, '..', sprintf('Release%d', r), 'L0', 'ADCP', [adcp_mooring_ID{adcp_ID} '_trimmed.mat']);
%     fprintf('saving %s...\n', trimDir)
%     save(trimDir, '-struct', "trim")
% end

% % Fetch Release times
% disp('fetching data to fit release 1 times')
% R1 = fetch_sig1k(L0dir, R1depTime(1), R1depTime(2));
% 
% disp('Saving Release 1')
% save([release_outRoot,'/Release1/', adcp_mooring_ID{adcp_ID}, '.mat'],'-struct','R1')
% 
% disp('fetching data to fit release 2 times')
% R23 = fetch_sig1k(L0dir, R23depTime(1), R23depTime(2));
% 
% disp('Saving Release 2')
% save([release_outRoot,'/Release2/', adcp_mooring_ID{adcp_ID}, '.mat'],'-struct','R23')
% 
% 









%% Functions

function loadSIG1K(rootDIR, fRoot, ATM_Time, ATM_Pressure, Config, hab, L0dir, filePrefix, deployTime, recoverTime, echo_mode, Descriptions, HeadingOffset)

files = dir([rootDIR,filesep,fRoot,'*.mat']);
Nf    = length(files);
%
Config.ATM_Time = ATM_Time;
Config.ATM_Pressure=ATM_Pressure;
%
% average the start/end pressures
ATM_Pressure = mean(ATM_Pressure);
%
% get bin sizing etc for currents
blnk = Config.Burst_BlankingDistance;
binw = Config.Burst_CellSize;
Nc   = double(Config.Burst_NCells);
mab  = hab+blnk+binw.*[1:Nc]';
%
%  get bin sizing etc for echo
if echo_mode
blnkE= Config.EchoSounder_BlankingDistance;
binwE= Config.EchoSounder_CellSize;
NcE  = double(Config.EchoSounder_NCells);
mabE = hab+blnkE+binwE.*[1:NcE]';
end
%
% save the config info
outFile = [L0dir, filePrefix, 'config.mat'];
if ~exist(L0dir,'dir')
    mkdir(L0dir)
end
save(outFile,'Config','Descriptions')
outFile = [L0dir, filePrefix, 'config.nc'];
if exist(outFile,'file')
    delete(outFile)
end
struct2nc(Config,outFile,'NETCDF4')
%
% limit archive file to 24hr length
maxDuration = 24*3600;
%
% get sample rate
fs = double(Config.Burst_SamplingRate);
%
% maximum number of samples per file
Nmax = fs*maxDuration;
% number of measurements in current file
N    = 0;
%
% flag to load next file
loadFlag = 1;
ii = 0;
outNf = 0;
fprintf('\n \n')
while ii<=Nf
    %
    if loadFlag
        ii  = ii+1;
        loadFlag=0;
        fileName = [fRoot,num2str(ii),'.mat'];
        fin = [files(ii).folder,'/',fileName];
        fprintf(['pre-processing:       %s \n'], fileName)
        load(fin)
        % use 5th beam time; the other beams are offset by dt = 1/(2*fs)
        t    = Data.IBurst_Time;
        nt   = length(t);
        % check if instrument is deployed
        is   = find(t>=deployTime,1,'first');
        if isempty(is)
            disp(['--> instrument is not yet deployed, skipping this file.'])
            loadFlag=1;
            continue
        end
        % check if instrument was pulled
        iStop= find(t>=recoverTime,1,'first');
        if ~isempty(iStop)
            disp(['trimming data at recovery time: nt = ',num2str(nt),'--> ',num2str(iStop-1)]);
            nt = iStop-1;
        end
    end
    %
    if N==0
        % initialize output structure
        out  = struct('Time',[],'Heading',[],'Pitch',[],'Roll',[],'Temperature',[],'Pressure',[],'Battery',[],...
                      'Velocity_Beam',[],...,'VelBeam2',[],'VelBeam3',[],'VelBeam4',[],'VelBeam5',[],...
                      'Amplitude_Beam',[],...'AmpBeam2',[],'AmpBeam3',[],'AmpBeam4',[],'AmpBeam5',[],...
                      'Correlation_Beam',[],...'CorBeam2',[],'CorBeam3',[],'CorBeam4',[],'CorBeam5',[],...,
                      'Velocity_East',[],'Velocity_North',[],'Velocity_Up',[],'Velocity_Error',[],...
                      'bin_mab',mab);
        if echo_mode
            out.bin_mab_Echo = mabE;
            out.Echo1=[];
            out.Echo2=[];
        end
    end
    %
    if N+nt-(is-1)<=Nmax
        ie   = nt;
        loadFlag = 1;
    else
        ie   = Nmax-N-(is-1);
        loadFlag = 0;
    end
    %
    out.Time       (1,N+1:N+(ie-(is-1)))     = Data.IBurst_Time(is:ie,1);
    out.Heading    (1,N+1:N+(ie-(is-1)))     = Data.Burst_Heading(is:ie,1);
    out.Pitch      (1,N+1:N+(ie-(is-1)))     = Data.Burst_Pitch(is:ie,1);
    out.Roll       (1,N+1:N+(ie-(is-1)))     = Data.Burst_Roll(is:ie,1);
    out.Temperature(1:2,N+1:N+(ie-(is-1)))   = [Data.IBurst_Temperature(is:ie,1)'; Data.Burst_Temperature(is:ie,1)'];
    out.Pressure   (1:2,N+1:N+(ie-(is-1)))   = [Data.IBurst_Pressure(is:ie,1)';    Data.Burst_Pressure(is:ie,1)'] - ATM_Pressure;
    out.Battery    (1,N+1:N+(ie-(is-1)))     = Data.IBurst_Battery(is:ie,1);
    out.Velocity_Beam    (1:Nc,N+1:N+(ie-(is-1)),1:5) = permute(cat(2,Data.Burst_Velocity_Beam(is:ie,:,1:Nc) ,Data.IBurst_Velocity_Beam(is:ie,1,1:Nc)), [3 1 2]) ;
    out.Amplitude_Beam   (1:Nc,N+1:N+(ie-(is-1)),1:5) = permute(cat(2,Data.Burst_Amplitude_Beam(is:ie,:,1:Nc),Data.IBurst_Amplitude_Beam(is:ie,:,1:Nc)), [3 1 2]);
    out.Correlation_Beam (1:Nc,N+1:N+(ie-(is-1)),1:5) = permute(cat(2,Data.Burst_Correlation_Beam(is:ie,:,1:Nc),Data.IBurst_Correlation_Beam(is:ie,:,1:Nc)), [3 1 2]);
% $$$     out.VelBeam1   (N+1:N+(ie-(is-1)),1:4,1:Nc)  = Data.Burst_Velocity_Beam(is:ie,1,1:Nc);
% $$$     out.VelBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_Beam(is:ie,2,1:Nc);
% $$$     out.VelBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_Beam(is:ie,3,1:Nc);
% $$$     out.VelBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_Beam(is:ie,4,1:Nc);    
% $$$     out.VelBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_Velocity_Beam(is:ie,1,1:Nc);
% $$$     out.AmpBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,1,1:Nc);
% $$$     out.AmpBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,2,1:Nc);
% $$$     out.AmpBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,3,1:Nc);
% $$$     out.AmpBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,4,1:Nc);    
% $$$     out.AmpBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_Amplitude_Beam(is:ie,1,1:Nc);
% $$$     out.CorBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,1,1:Nc);
% $$$     out.CorBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,2,1:Nc);
% $$$     out.CorBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,3,1:Nc);
% $$$     out.CorBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,4,1:Nc);    
% $$$     out.CorBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_Correlation_Beam(is:ie,1,1:Nc);
    out.Velocity_East    (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,1,1:Nc),[3,1,2]);
    out.Velocity_North   (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,2,1:Nc),[3,1,2]);
    out.Velocity_Up      (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,3,1:Nc),[3,1,2]);
    out.Velocity_Error   (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,4,1:Nc),[3,1,2]);
    if echo_mode
    out.Echo1            (1:NcE,N+1:N+(ie-(is-1))) = Data.Echo1Bin1_1000kHz_Echo(is:ie,1:NcE)';
    out.Echo2            (1:NcE,N+1:N+(ie-(is-1))) = Data.Echo2Bin1_1000kHz_Echo(is:ie,1:NcE)';
    end
    %
    N = N+ie-(is-1);
    is= ie+1;
    %
    if N==Nmax || (ii==Nf & loadFlag)
        % process data
        minC  = min( out.Correlation_Beam, [], 3);%min( cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        minA  = min( out.Amplitude_Beam  , [], 3);  %cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        maxRNG= out.Pressure(2,:)*cosd(25)-3*binw;
        qcFlag= (minC>50 & minA>30 & maxRNG>mab);
        %
        out.Correlation_Minimum = minC;
        out.Amplitude_Minimum = minA;
        out.maxRNG = maxRNG;
        out.qcFlag = qcFlag;
        %
        disp('not applying heading offset correction, velocities are relative to magnetic north/south')
% $$$         % now rotate to xyz and ENU
% $$$         [out,Config2,beam2xyz] = beam2earth_sig1000_DG_FRFcoords(out,Config,'');        
        out.HeadingOffset = HeadingOffset;
% $$$         out.beam2xyz = beam2xyz;
        %
        % save output
        outNf = outNf+1;
        outFileName = sprintf([filePrefix,'%03d.mat'],outNf);
        fout  = [L0dir,outFileName];
        fprintf(['saving output file:   %s \n'],outFileName)
        save(fout,'-struct','out')
        outFileName = sprintf([filePrefix,'%03d.nc'],outNf);
        fout = [L0dir,outFileName];
        if exist(fout,'file')
            delete(fout)
        end
        struct2nc(out,fout,'NETCDF4')
        %
        N = 0;
        %
        if (ii==Nf & loadFlag )
        fprintf('done! \n')
        break
        end
    end
    %
    %
end
end

% EOF


function time_average_and_rotate_sig1000_RDI_matrix_format_function(Config, L0dir, filePrefix, dtAvg, echo_mode, L0FRoot)

%
load([L0dir,filesep,filePrefix,'config.mat'])
fs = double(Config.Burst_SamplingRate);
Nc = Config.Burst_NCells;
%
files    = dir([L0dir,filesep,filePrefix,'*.mat']);
fNameCell=extractfield(files,'name');
files    = files(~contains(fNameCell,'config') & ~contains(fNameCell,'min.mat'));
Nf       = length(files);
%
fprintf('\n \n')
%
% define filter parameters
fc = 1/dtAvg;% frequency cutoff
Ns = fs/fc;% filter half-width
if mod(Ns,2), Ns=Ns+1; end
Fw = 2*Ns+1;% filter width
F  = hanning(Ns+1);% Hann-function filter
F  = F./sum(F);% normalize
%
for ii= 1:Nf
    %
    % load the raw data
    inFileName = sprintf([filePrefix,'%03d.mat'],ii);
    fin  = [L0dir,inFileName];
    fprintf(['loading file:   %s \n'],inFileName)
    in = load(fin,'Velocity_East','Velocity_North','Velocity_Up','Velocity_Error','Correlation_Minimum','Amplitude_Minimum','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','Temperature','bin_mab');
    %
    if echo_mode
        echo = load(fin,'bin_mab_Echo','Echo1','Echo2');
        fieldNames = fields(echo);
        for jj = 1:length(fieldNames)
            in.(fieldNames{jj}) = echo.(fieldNames{jj});
        end
    end

    % correct any time gaps
    struct_fieldnames = fields(in);
  
    for fieldnum = 1:length(struct_fieldnames)
        sff = struct_fieldnames{fieldnum};
        if strcmp(sff, 'Time')
            continue
        else
            dum = in.(sff);
            if max(size(dum)) == length(in.Time)
                [t_fix, x_fix] = correct_burst(in.Time', dum');
                in.(sff) = x_fix';
            end
        end
    end
        
    in.Time = t_fix';
    %
    if ii==1
        % no padding
        tP   = [];
        qcP  = [];
        vb1p = [];
        vb2p = [];
        vb3p = [];
        vb4p = [];
        aP   = [];
        cP   = [];
        e1p  = [];
        e2p  = [];
        hP   = [];
        pP   = [];
        rP   = [];
        Pp   = [];
        Tp   = [];
    end
    %
    np = length(tP);
    nt = length(in.Time);
    % fractional number of full filter segements
    ns = (nt+np)/Ns;
    % number of whole segments
    N  = floor( ns );
    % remainder
    ns = rem(ns,1);
    %
    % average Iburst and burst Temp/Pres
    in.Pressure    = nanmean(in.Pressure,1);
    in.Temperature = nanmean(in.Temperature,1);
    %
    %
    % perform time-convolution (filling nan's in the mask), conv2 x 3 (nans=0, qcFlag, ones)
    on  = [qcP, in.qcFlag];
    o   = ones(1,size(on,2));
    %
    % convolve NAN matrix for normalization
    norm0   = conv2(o,F','same');
    norm    = conv2(on,F','same');
    % throw out regions with <%25 coverage
    avgFlag = norm./norm0;
    norm(norm<=0.25)=inf;
    %
    % convolve signals
    t1    = [tP, in.Time];
    vb1   = conv2(1, F, [vb1p, in.Velocity_East ].*on,'same')./norm;
    vb2   = conv2(1, F, [vb2p, in.Velocity_North].*on,'same')./norm;
    vb3   = conv2(1, F, [vb3p, in.Velocity_Up   ].*on,'same')./norm;
    vb4   = conv2(1, F, [vb4p, in.Velocity_Error].*on,'same')./norm;
    a     = conv2(1, F, [aP  , in.Amplitude_Minimum].*on,'same')./norm;
    c     = conv2(1, F, [cP  , in.Correlation_Minimum].*on,'same')./norm;        
    head  = conv2(1, F, [hP  , in.Heading ]    ,'same')./norm0;
    pitch = conv2(1, F, [pP  , in.Pitch   ]    ,'same')./norm0;
    roll  = conv2(1, F, [rP  , in.Roll    ]    ,'same')./norm0;
    P     = conv2(1, F, [Pp  , in.Pressure]    ,'same')./norm0;
    T     = conv2(1, F, [Tp  , in.Temperature] ,'same')./norm0;
    if echo_mode
        echo1 = conv2(1, F, [e1p , in.Echo1]       ,'same')./norm0;
        echo2 = conv2(1, F, [e2p , in.Echo2]       ,'same')./norm0;
    end
    %
    %
    avg = struct('Time',t1(Ns:Ns:Ns*(N-1)),'Velocity_East',vb1(:,Ns:Ns:Ns*(N-1)),'Velocity_North',vb2(:,Ns:Ns:Ns*(N-1)),'Velocity_Up',vb3(:,Ns:Ns:Ns*(N-1)),'Velocity_Error',vb4(:,Ns:Ns:Ns*(N-1)),'Amplitude_Minimum',a(:,Ns:Ns:Ns*(N-1)),'Correlation_Minimum',c(:,Ns:Ns:Ns*(N-1)),'Heading',head(:,Ns:Ns:Ns*(N-1)),'Pitch',pitch(:,Ns:Ns:Ns*(N-1)),'Roll',roll(:,Ns:Ns:Ns*(N-1)),'Pressure',P(Ns:Ns:Ns*(N-1)),'Temperature',T(Ns:Ns:Ns*(N-1)),'qcFlag',avgFlag(:,Ns:Ns:Ns*(N-1)),'HeadingOffset',in.HeadingOffset,'bin_mab',in.bin_mab);    
    %
    if echo_mode
        avg(1).bin_mab_Echo=in.bin_mab_Echo;
        avg(1).Echo1       =echo1(:,Ns:Ns:Ns*(N-1));
        avg(1).Echo2       =echo2(:,Ns:Ns:Ns*(N-1));
    end
        
    % 
    % need to rotate from magnetic to true?
% $$$     [avg,Config2,beam2xyz] = beam2earth_sig1000_DG(avg,Config,'',in.HeadingOffset);
% $$$     avg.beam2xyz=beam2xyz;
    %
    if ii==1
        out = avg;
    else
        vars = fields(out);
        for jj = 1:length(vars);
            var = vars{jj};
            if ismember(var,{'HeadingOffset','bin_mab','bin_mab_Echo','beam2xyz'})
                continue
            end
            eval(['out.',var,'= cat(2,out.',var,',avg.',var,');'])
        end
    end
  
    %
    % get pad data for next file
    tP   = in.Time       (:,nt-Ns*(1+ns)+1:nt);
    qcP  = in.qcFlag     (:,nt-Ns*(1+ns)+1:nt);
    vb1p = in.Velocity_East    (:,nt-Ns*(1+ns)+1:nt);
    vb2p = in.Velocity_North   (:,nt-Ns*(1+ns)+1:nt);
    vb3p = in.Velocity_Up      (:,nt-Ns*(1+ns)+1:nt);
    vb4p = in.Velocity_Error   (:,nt-Ns*(1+ns)+1:nt);
    aP   = in.Amplitude_Minimum  (:,nt-Ns*(1+ns)+1:nt);
    cP   = in.Correlation_Minimum(:,nt-Ns*(1+ns)+1:nt);    
    hP   = in.Heading    (nt-Ns*(1+ns)+1:nt);
    pP   = in.Pitch      (nt-Ns*(1+ns)+1:nt);
    rP   = in.Roll       (nt-Ns*(1+ns)+1:nt);
    Tp   = in.Temperature(nt-Ns*(1+ns)+1:nt);
    Pp   = in.Pressure   (nt-Ns*(1+ns)+1:nt);
    if echo_mode
        e1p  = in.Echo1      (:,nt-Ns*(1+ns)+1:nt);
        e2p  = in.Echo2      (:,nt-Ns*(1+ns)+1:nt);    
    end
    %
end
fprintf(['saving: %s \n'],[L0dir,L0FRoot])
if ~exist(L0dir,'dir')
    mkdir(L0dir)
end
save([L0dir,L0FRoot],'-struct','out')
fprintf('done! \n')
outFileName = [L0dir,L0FRoot,'.nc'];
if exist(outFileName,'file')
    delete(outFileName)
end
struct2nc(out,outFileName,'NETCDF4')
end

% EOF


function waves = estimate_wave_bulk_stats_SIG1000_RDI_matrix_format_function(Config, L0dir, filePrefix, L0FRoot, hab, L1dir, L1FRoot)
%
dtBurst = 1800;% seconds
dtEns   = 512 ;% seconds
rho0    = 1027.5;% kg/m^3
g       = 9.81;% m/s^2
%
%
% get sampling information from config file
load([L0dir,filePrefix,'config.mat'])
fs = double(Config.Burst_SamplingRate);
Nc = Config.Burst_NCells;
%
% load L0-file with time-averages
currents = load([L0dir,filesep,L0FRoot,'.mat']);
%
% sample time step
dt = 1/fs;
%
% number of samples to average
Na = dtBurst*fs;
% number of samples in each ensemble
Ne    = dtEns*fs;
olap  = 2/3;
chnks = (Na-Ne*olap-1)/(Ne*(1-olap));
% $$$ Nfreq =floor(Na/(chnks-(chnks-1)*olap))/2;
%
% get structure with all files in archive
files = dir([L0dir,filePrefix,'*.mat']);
fNameCell=extractfield(files,'name');
files = files(~contains(fNameCell,'config') & ~contains(fNameCell,'min.mat'));
Nf    = length(files);
%
% initialize counter...
fprintf('\n \n')
ensNum = 1;
waves  = struct();
for ii= 1:Nf
    %
    % load the raw data
    inFileName = sprintf([filePrefix,'%03d.mat'],ii);
    fin  = [L0dir,inFileName];
    fprintf(['loading file:   %s \n'],inFileName)
% $$$     in = load(fin,'VelEast','VelNorth','VelUp1','VelUp2','VelBeam5','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','mab');
    in = load(fin,'Velocity_East','Velocity_North','Velocity_Up','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','bin_mab');
    %
    if ii==1
        mab = in.bin_mab;
    end
    %
    % get number of ensembles in current file:
    nt = length(in.Time);
    N = floor(nt/Na);
    % 
    % loop over ensembles
    for jj = 1:N
        %
        % check to see which bins are in the water
        QC   = in.qcFlag(:,(jj-1)*Na+1:Na*jj);
        fGood= sum(QC,2)/Na;
        % only use bins in water >75 % of the time
        bins = find(fGood>0.75);
        Nb   = length(bins);
        %
        % mean time to nearest minute
        t = in.Time((jj-1)*Na+1:Na*jj);
        tavg = round(mean(t)*1440)/1440;
        t = (t-t(1))*86400;
        %
        waves.Time(1,ensNum) = tavg;
        if Nb<2
            disp('not enough good data')
            waves.Hs(1,ensNum) = [];
            waves.Tm(1,ensNum) = [];
            waves.Suu(:,ensNum)= [];
            waves.Svv(:,ensNum)= [];
            waves.Spp(:,ensNum)= [];
            waves.Spu(:,ensNum)= [];
            waves.Spp(:,ensNum)= [];
            waves.Spv(:,ensNum)= [];
% $$$             waves.mean_dir(1:2,ensNum)   =nan;
% $$$             waves.mean_spread(1:2,ensNum)=nan;
            waves.mSxx(1,ensNum) = [];
            waves.mSxy(1,ensNum) = [];
            waves.mSyy(1,ensNum) = [];
            waves.Z2(:,ensNum)   = [];
            ensNum = ensNum+1;
            continue
        end
        %
        % allocate current variables: (note, bad QC are zeros not nans, and not interpolated)
        U = in.Velocity_East (bins,(jj-1)*Na+1:Na*jj);
        V = in.Velocity_North(bins,(jj-1)*Na+1:Na*jj);
        W = in.Velocity_Up   (bins,(jj-1)*Na+1:Na*jj);
        P = in.Pressure(2,(jj-1)*Na+1:Na*jj)*1e4/(rho0*g);
        %
        % depth of pressure sensor
        dpthP= mean(P,2);% not accounting for sensor height above bed
        % assuming sensor is at "hab" above bed level        
        dpth = dpthP + hab;
        % depth of bins
        dpthU= dpthP-mab(bins);
        %
% $$$         % Use WAFOS package:
% $$$         % first pass, just use (U,V,W,P):
% $$$         xn  = [t', (P'-dpthP)*1e4, U', V', W'];
% $$$         pos0= [0, 0,-dpthP, 9, 1];
% $$$         posU= [zeros(Nb,1),zeros(Nb,1),-dpthU,ones(Nb,1)*10, ones(Nb,1)];
% $$$         posV= [zeros(Nb,1),zeros(Nb,1),-dpthU,ones(Nb,1)*11, ones(Nb,1)];
% $$$         posW= [zeros(Nb,1),zeros(Nb,1),-dpthU,ones(Nb,1)*12, zeros(Nb,1)];        
% $$$         pos = [pos0;posU;posV;posW];
% $$$         [Sd,D,Sw,Fcof,Gwt,Sxy,Sxy1] = dat2dspec_DG(xn,pos,dpth,Ne/4,90,'MLM');% 256 frequency bins & 90 directional bins
% $$$         %
% $$$         % quick plot in polar coords
% $$$         x = Sd.f'.*cos(Sd.theta);
% $$$         y = Sd.f'.*sin(Sd.theta);
% $$$         figure, contourf(x,y,log10(Sd.S),[-4:0.05:3],'edgecolor','none'),colormap(cmocean('thermal'))
% $$$         %
% $$$         %
% $$$         % estimate mean direction
% $$$         df = Sd.f(2)-Sd.f(1);
% $$$         dT = D.theta(2)-D.theta(1);
% $$$         I  = find(Sd.f>=1/20 & Sd.f<=1/4);
% $$$         A1 = sum( D.S.*cos(D.theta)*dT, 1);
% $$$         B1 = sum( D.S.*sin(D.theta)*dT, 1);
% $$$         %
        %
        % estimate bulk wave statistics on ~7 min chuncks, 2/3 overlap
        [Suu,fq]=welch_method(U',dt,chnks,olap);
        [Svv,fq]=welch_method(V',dt,chnks,olap);
        [Spp,fq]=welch_method(P',dt,chnks,olap);
        [Suv,fq]=welch_cospec(U',V',dt,chnks,olap);        
        [Spu,fq]=welch_cospec(repmat(P',1,Nb),U',dt,chnks,olap);
        [Spv,fq]=welch_cospec(repmat(P',1,Nb),V',dt,chnks,olap);
        %
        % note: still has zero frequency!
        fq = fq(2:end);
        Suu= Suu(2:end,:);
        Svv= Svv(2:end,:);
        Spp= Spp(2:end,:);
        Suv= Suv(2:end,:);
        Spu= Spu(2:end,:);
        Spv= Spv(2:end,:);
        %
        % depth correction
        om = 2*pi*fq;
        lom=length(om);
        %       
        L = disper(om.', dpth);
        k = 2.*pi ./ L;% don't have the wavenumber() funciton but i included my disper in the end
        k = k';

        %k=wavenumber(om.',dpth); % these are the radian wave numbers (rad/m)
        %
        % convert pressure to surface elevation
        cP     =  cosh(k.*dpth)./ cosh(k.*(dpth-dpthP)); % SePP=Spp.*cP.^2
        % convert velocity to surface elevation
        cU2eta = (om./(g*k)).*cosh(k.*dpth)./cosh(k.*(dpth-dpthU'));% SeUU = SeUU.*cU.^2        
        % convert velocity to surface velocity
        cU2srf =              cosh(k.*dpth)./cosh(k.*(dpth-dpthU'));% SeUU = SeUU.*cU.^2        
        %
        % estimate max resolvable frequency where depth bin is below wave-length/(pi)
        fmax_vs_bin = fq(Ne/2-sum(dpthU'>=1./(k)))';
        % nan surface transformed fields above fmax
        %        Spp(fq>fmax(1))=nan;
        for kk=1:Nb
        Suu(fq>fmax_vs_bin(kk),kk)=nan;
        Svv(fq>fmax_vs_bin(kk),kk)=nan;
        Suv(fq>fmax_vs_bin(kk),kk)=nan;
        Spu(fq>fmax_vs_bin(kk),kk)=nan;
        Spv(fq>fmax_vs_bin(kk),kk)=nan;
        end
        % surface velocity spectra
        SUU = Suu.*cU2srf.^2;
        SVV = Svv.*cU2srf.^2;
        SUV = Suv.*cU2srf.^2;
        SPU = repmat(cP,1,Nb).*Spu.*cU2srf;
        SPV = repmat(cP,1,Nb).*Spv.*cU2srf;        
        %
        % map to surface elevation spectra
        SeUU = Suu.*cU2eta.^2;
        SeVV = Svv.*cU2eta.^2;
        SeUV = Suv.*cU2eta.^2;
        SePP = Spp.*cP.^2;
        SePU = repmat(cP,1,Nb).*Spu.*cU2eta;
        SePV = repmat(cP,1,Nb).*Spv.*cU2eta;
        %
        % estimate best guess for surface elevation spectrum
% $$$         % first, average the velocity spectra
% $$$         SeAVG = nanmean(SeUU+SeVV,2);
% $$$         % use z-score as a weighting function
% $$$         norm   = [tanh(pi*SePP./SeAVG), ones(Ne/2,1)];
% $$$         norm(2:end-1,:) = (norm(1:end-2,:) + norm(2:end-1,:) + norm(3:end,:))/3;
% $$$         SeBest = nansum([SeAVG SePP].*norm,2)./sum(norm,2);
        %
        if hour(tavg)==12
            fig1 = figure;
            plt = semilogy(fq,nanmean(SeUU+SeVV,2),'-r',fq,SePP,'--k','linewidth',3);
            set(gca,'xlim',[1/(30*60) 1/3])
            grid on
            xlabel('cycles per second','interpreter','latex')
            ylabel('m$^2$/Hz','interpreter','latex')
            title(datestr(tavg),'interpreter','latex')
            h = legend(plt,'$S_{uu}+S_{vv}$','$S_{\eta\eta}$');
            set(h,'location','northeast','interpreter','latex')
            figname = sprintf('%s/figures/%s_spectra_%s.pdf',L1dir,L1FRoot,datestr(tavg,'mmddHHMM'));
            if ~exist([L1dir,filesep,'figures'])
                eval(['!mkdir -p ',[L1dir,filesep,'figures/']])
            end
            exportgraphics(fig1,figname)
            close(fig1)
        end
        %
        % surface corrected total velocity spectra
        SUtot = SeUU+SeVV;
        Z2 = nanmean(SePP./SUtot,2);
        %
        %
% $$$         % create fake (P,U,V), rotated by -70deg=160 nautical, or +70deg=20nautical
% $$$         Nb = 1;
% $$$         a  = 1;% 1 meter
% $$$         ang= 70*pi/180+pi;% direction from--> direction to
% $$$         P  = a*cos(2*pi/10.*(t-t(1))*86400);
% $$$         U  = P.*cos(ang) ;
% $$$         V  = P.*sin(ang);
% $$$         %
% $$$         % estimate bulk wave statistics on 10 min chuncks, no overlap
% $$$         [Suu,fq]=welch_method(U,dt,chnks,olap);
% $$$         [Svv,fq]=welch_method(V,dt,chnks,olap);
% $$$         [Spp,fq]=welch_method(P,dt,chnks,olap);
% $$$         [Suv,fq]=welch_cospec(U,V,dt,chnks,olap);        
% $$$         [Spu,fq]=welch_cospec(repmat(P,1,Nb),U,dt,chnks,olap);
% $$$         [Spv,fq]=welch_cospec(repmat(P,1,Nb),V,dt,chnks,olap);
% $$$         %
% $$$         % note: still has zero frequency!
% $$$         fq  = fq(2:end);
% $$$         SUU = Suu(2:end,:);
% $$$         SVV = Svv(2:end,:);
% $$$         SePP= Spp(2:end,:);
% $$$         SUV = Suv(2:end,:);
% $$$         SPU = Spu(2:end,:);
% $$$         SPV = Spv(2:end,:);
% $$$         %
        % estimate bulk stats
        % I think there is a sign issue here
        coPU =  SPU;
        coPV =  SPV;
        coUV =  SUV;
        r2d = 180/pi;
        %
        % get a1
        a1      = coPU ./ sqrt( SePP .* ( SUU + SVV ) );
        b1      = coPV ./ sqrt( SePP .* ( SUU + SVV ) );
        dir1    = r2d* ( atan2(b1,a1) );          
        spread1 = r2d* ( sqrt( 2 .* ( 1-sqrt(a1.^2 + b1.^2) ) ) );
        %
        % average over wind-wave band
        df= fq(2)-fq(1);
        I = find(fq>=1/20 & fq<=1/4);
        m0 = nansum(SePP(I)*df);
        m1 = nansum(fq(I).*SePP(I)*df);        
        ma1= nansum(a1(I,:).*SePP(I)*df,1)/m0;
        mb1= nansum(b1(I,:).*SePP(I)*df,1)/m0;
        mdir1=rem(90+180-r2d*atan2(mb1,ma1),360);% nautical 
        mspread1 = r2d*sqrt(abs(2*(1-(ma1.*cos(mdir1/r2d) + mb1.*sin(mdir1/r2d)))));
        %
        % get a2 b2
        a2 = (SUU - SVV) ./ (SUU + SVV);
        b2 = 2 .* coUV   ./ (SUU + SVV);
        spread2 = r2d*sqrt(abs(0.5-0.5*(a2.*cos(2.*dir1/r2d)+b2.*sin(2.*dir1/r2d))));
        %
        ma2      = nansum(a2(I,:).*SePP(I)*df,1)/m0;
        mb2      = nansum(b2(I,:).*SePP(I)*df,1)/m0;
        dir2     = r2d/2*atan2(b2,a2);
        mdir2    = 90-(r2d/2*atan2(mb2,ma2));% nautical
        mspread2 = r2d*sqrt(abs(0.5-0.5*(ma2.*cos(2.*mdir1/r2d)+mb2.*sin(2.*mdir1/r2d))));
        %
        % radiation stresses
        % Radiation Stress Estimates
        C = om./k;
        Cg = 0.5*(g*tanh(k.*dpth)+g*k.*dpth.*(sech(k.*dpth).^2))./sqrt(g*k.*tanh(k.*dpth));
        %
        % cartesian
        Sxx = ( (1.5 + 0.5*a2) .* (Cg./C) - 0.5 ) .* SePP;
        Syy = ( (1.5 - 0.5*a2) .* (Cg./C) - 0.5 ) .* SePP;
        Sxy = 0.5*b2 .* (Cg./C) .* SePP;
        %
        mSxx = sum( SePP(I).*Sxx(I,:).*df )./m0;
        mSyy = sum( SePP(I).*Syy(I,:).*df )./m0;
        mSxy = sum( SePP(I).*Sxy(I,:).*df )./m0;                
        %
        %
        Hs = 4*sqrt(nansum(df*SePP(I)));
        Tm = m0/m1;
        %
        % log the hourly averaged
        if ii==1 & jj==1
            waves.frequency  = fq(1:I(end))';
            waves.wavenumber = k(1:I(end))';            
        end
        waves.Hs     (1,ensNum)   = Hs;
        waves.Tm     (1,ensNum)   = Tm;
        waves.depth  (1,ensNum)   = dpth;
        waves.Spp    (:,ensNum)   = SePP(1:I(end));
        waves.Suu    (:,ensNum)   = nanmean(SUU(1:I(end),:),2);
        waves.Svv    (:,ensNum)   = nanmean(SVV(1:I(end),:),2);
        waves.Spu    (:,ensNum)   = nanmean(SPU(1:I(end),:),2);
        waves.Spv    (:,ensNum)   = nanmean(SPV(1:I(end),:),2);
        waves.mdir   (1:2,ensNum) = [nanmean(mdir1)   ,nanmean(mdir2)];
        waves.mspread(1:2,ensNum) = [nanmean(mspread1),nanmean(mspread2)];
        waves.a1     (:,ensNum)   = nanmean(a1,2);
        waves.b1     (:,ensNum)   = nanmean(b1,2);        
        waves.a2     (:,ensNum)   = nanmean(a2,2);
        waves.b2     (:,ensNum)   = nanmean(b2,2);        
        waves.mSxx   (1,ensNum)     = nanmean(mSxx);
        waves.mSxy   (1,ensNum)     = nanmean(mSxy);
        waves.mSyy   (1,ensNum)     = nanmean(mSyy);
        waves.Z2     (:,ensNum)    = Z2;
        ensNum = ensNum+1;
    end
    %
end
fprintf(['saving: %s \n'],[L1dir,L1FRoot])
if ~exist(L1dir,'dir')
    mkdir(L1dir)
end
save([L1dir,L1FRoot,'_currents.mat'],'currents')
save([L1dir,L1FRoot,'_waves.mat'],'waves')
fprintf('done! \n')

currentsFileName = [L1dir,L1FRoot,'_currents.nc'];
if exist(currentsFileName,'file')
    delete(currentsFileName)
end
struct2nc(currents,currentsFileName,'NETCDF4')
wavesFileName = [L1dir,L1FRoot,'_waves.nc'];
if exist(wavesFileName,'file')
    delete(wavesFileName)
end
struct2nc(waves,wavesFileName,'NETCDF4')
end

% EOF

function L = disper(T, h)
%
%   USAGE: L = disper(T, h)
%   
%   Calculates Wavelength for intermdiate depth surface gravity waves
%   
%   Inputs:
%            Wave period, T [s] (float or vector)
%            Water depth, h [m] (float or vector)
%               
%   Outputs:
%           Wavelength, L [m] (scalar)
%               IF T and h are different sizes ( [N x 1] and [M x 1] )
%               THEN L:
%                           Rows --> T
%                           Cols --> h

disp('Estimating L using linear dispersion relation')


if ~(isvector(T) && isvector(h))
    error('T and h must be vectors or scalars')
end


if isscalar(T) || isscalar(h)
    % scalar (boring)

elseif isequal(size(T), size(h))
    % dont need to do anything

else
    T = T(:);
    h = h(:)';
    disp('performing implicit expansion')
end



w2 = (2*pi./T).^2; % linear dispersion relation
g = 9.81;

k = w2 ./ g; % deep water initial guess

for i = 1:20
    f = g.*k.*tanh(k.*h) - w2;
    dfdk = g.*tanh(k.*h) + g.*k.*h.*sech(k.*h).^2;
    k = k - f./dfdk;
end

L = 2*pi./k;
end
% EOF

function Data = fetch_sig1k(inputDir, startTime, endTime)
%%
% 
% USAGE: Opens signature1000 ADCP files based on deployment times
%
%   inputDir: Directory where the .mat's can be found
%   statTime: Start of deployment (I think it can be in any date format)
%   endTIme: End of deployment
%
%   % Returns: Structure 'Data' with sig1K data
%
%%

files = dir([inputDir, '*.mat']);

% fetch files
for i = 1:length(files)
    % skip config and 5min avg
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        % find files with appropriate times
        if any((Times.Time >= datenum(startTime)) & Times.Time<= datenum(endTime))
            fprintf('loading file %s...\n', files(i).name)
            Sig =  load([files(i).folder, filesep, files(i).name]);
            % look for QC flag
            fields = fieldnames(Sig);
            for j = 1:length(fields)
                    field = fields{j};
             
                if any(strcmp(fields, 'qcFlag'))
                    % apply QC flag
                    if contains(field, 'Velocity') && isequal(size(Sig.(field)), size(Sig.qcFlag))
                        dum = Sig.(field);
                        dum(~Sig.qcFlag) = NaN;
                        Sig.(field) = dum;
                    end
                end
                if ~exist('Data','var')
                Data = struct();
                end
                
                if ~isfield(Data, field)
                    Data.(field) = Sig.(field);
                elseif size(Sig.(field), 2) == size(Sig.Time, 2)
                    Data.(field) = [Data.(field) Sig.(field)];
                end
            end
        end
    end
end

% Data.Time = Sig.Time;
% Data.Velocity_East = Sig.Velocity_East';
% Data.Velocity_North = Sig.Velocity_North';
% Data.Velocity_Up = Sig.Velocity_Up';
% Data.Velocity_Beam1 = Sig.Velocity_Beam(:,:,1)';
% Data.Velocity_Beam2 = Sig.Velocity_Beam(:,:,2)';
% Data.Velocity_Beam3 = Sig.Velocity_Beam(:,:,3)';
% 
% Data.Correlation_Minimum = min(Sig.Correlation_Beam(:,:,1)', min(Sig.Correlation_Beam(:,:,2)', Sig.Correlation_Beam(:,:,3)'));
% Data.Amplitude_Minimum = min(Sig.Amplitude_Beam(:,:,1)', min(Sig.Amplitude_Beam(:,:,2)', Sig.Amplitude_Beam(:,:,3)'));
% Data.Pressure = Sig.Pressure(2,:)';
% Data.Heading = Sig.Heading(:, :)';
% Data.Pitch = Sig.Pitch(:,:)';
% Data.Roll = Sig.Roll(:,:)';
% Data.bin_mab = Sig.Config;



% Summary figure
figure
img = imagesc(datetime(Data.Time, 'convertfrom', 'datenum'), Data.bin_mab, Data.Velocity_East);
mask = isnan(Data.Velocity_East);
set(img, 'Alphadata', ~mask)
set(gca, 'YDir', 'normal')
cb = colorbar;
ylabel('Meters above bottom')
ylabel(cb, 'East Velocity, [m/s]', 'FontSize', 18)
set(gca, 'fontsize', 18)
clim([-0.1 0.1])
colormap(cmocean('balance'))
grid minor


figure
img = imagesc(datetime(Data.Time, 'convertfrom', 'datenum'), Data.bin_mab, Data.Velocity_North);
mask = isnan(Data.Velocity_North);
set(img, 'Alphadata', ~mask)
set(gca, 'YDir', 'normal')
cb = colorbar;
ylabel('Meters above bottom')
ylabel(cb, 'North Velocity, [m/s]', 'FontSize', 18)
set(gca, 'fontsize', 18)
clim([-0.1 0.1])
colormap(cmocean('balance'))
grid minor


figure
img = imagesc(datetime(Data.Time, 'convertfrom', 'datenum'), Data.bin_mab, Data.Velocity_Up);
mask = isnan(Data.Velocity_Up);
set(img, 'Alphadata', ~mask)
set(gca, 'YDir', 'normal')
cb = colorbar;
ylabel('Meters above bottom')
ylabel(cb, 'Up Velocity, [m/s]', 'FontSize', 18)
set(gca, 'fontsize', 18)
clim([-0.1 0.1])
colormap(cmocean('balance'))
grid minor

Data.Velocity_East = Data.Velocity_East';
Data.Velocity_North = Data.Velocity_North';
Data.Time = Data.Time';

end



% function struct2nc(x,ncfile,ncfiletype,deflate_lev)
% % STRUCT2NC writes all float,double and character vars to netcdf
% % Usage: struct2nc(x,ncfile,[ncfiletype],[deflate_lev])
% % x = structure
% % ncfile = name of netcdf output file (e.g. 'test.nc')
% % ncfiletype = netcdf file type (e.g. 'classic','netcdf4_classic')
% % deflate_lev = deflate level (0-9, 0 is none)
% %
% % This function writes all 'double','single' and 'char' variables
% % to NetCDF using the native Matlab NetCDF interface.  It skips all
% % other classes in the struct (e.g. structs, cell arrays, etc).  It
% % also only handles scalar, 1D, 2D, and 3D arrays currently, although
% % this could easily be extended.
% 
% if nargin==2
%     ncfiletype='classic';
%     deflate_lev=0;
% elseif nargin==3
%     switch ncfiletype
%         case {'netcdf4','netcdf4_classic'}
%             deflate_lev=6;
%         otherwise
%             deflate_lev=0;
%     end
% end
% 
% % Remove old output file
% if exist(ncfile,'file')
%     delete(ncfile)
% end
% 
% s      = fieldnames(x);
% layers = length(x);
% k=0;
% % create variables first, but don't write data
% for l=1:layers
%     for i=1:length(s)
%         vname=char(s(i));
%         if l>1
%            uname=[vname,num2str(l-1)];
%         else
%            uname=vname;
%         end
%             var   =x(l).(vname);
%             vtype = class(var);
%             vshape= size(var);
%             ndims = length(vshape);
%             vlen  = length(var(:));
%             if sum(vshape)==0
%                 continue
%             end
%     switch vtype;
%       case {'double','single','int32','uint8','uint64','logical'},
%         if strcmp(vtype,'logical')
%             x(l).(vname) = double(x(l).(vname));
%             vtype = 'double';
%         end
%             if vlen==1,
%                 nccreate(ncfile,uname,...
%                     'Datatype',vtype,'format',ncfiletype);
%                 k=k+1;
%                 vnames{k}=vname;
%                 unames{k}=uname;
%             else
%                 if min(vshape)==1,
%                     nccreate(ncfile,uname,...
%                         'Datatype',vtype,...
%                         'DeflateLevel',deflate_lev,...
%                         'Dimensions',{[uname '1'] vlen},...
%                         'format',ncfiletype);
%                     k=k+1;
%                     vnames{k}=vname;
%                     unames{k}=uname;                    
%                 elseif ndims==2,
%                     nccreate(ncfile,uname,...
%                         'Datatype',vtype,...
%                         'DeflateLevel',deflate_lev,...
%                         'Dimensions',{[uname '1'] vshape(1) [uname '2'] vshape(2)},...
%                         'format',ncfiletype);
%                     k=k+1;
%                     vnames{k}=vname;
%                     unames{k}=uname;                    
%                 elseif ndims==3,
%                     nccreate(ncfile,uname,...
%                         'Datatype',vtype,...
%                         'DeflateLevel',deflate_lev,...
%                         'Dimensions',...
%                         {[uname '1'] vshape(1) [uname '2'] vshape(2) [uname '3'] vshape(3)},...
%                         'format',ncfiletype);
%                     k=k+1;
%                     vnames{k}=vname;
%                     unames{k}=uname;                    
%                 else,
%                     disp('Skipping variable with more than 3 dimensions');
%                 end
%             end
%        case {'char'},
%          if min(vshape)==1,
%             nccreate(ncfile,uname,...
%                 'Datatype',vtype,...
%                 'Dimensions',{[uname '1'] vlen},.....
%                 'format',ncfiletype);
%             k=k+1;
%             vnames{k}=vname;
%             unames{k}=uname;
%          elseif ndims==2,
%              nccreate(ncfile,uname,...
%                       'Datatype',vtype,...
%                       'Dimensions',{[uname '1'] vshape(1) [uname '2'] vshape(2)},...
%                       'format',ncfiletype);
%              k=k+1;
%              vnames{k}=vname;
%              unames{k}=uname;                    
%          end
%         otherwise,
%             disp(['skipping ' vname])
%     end
% end
% end
% %write all the data at the end
% for l=1:layers
%     for i=1:length(vnames)
%         ncwrite(ncfile,unames{i},x(l).(vnames{i}));
%     end
% end
% end