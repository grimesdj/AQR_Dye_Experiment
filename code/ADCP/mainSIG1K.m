
% mainSIG1K.m
% 
%   USAGE: (script) loads and processes raw sig1K data into L0
%   '.mat' files
% 
%    REQUIRES USER TO ENTER:
%       adcp_ID         = index of ADCP to be used
%       adcp_file_roots = file roots for target ADCPs
%       adcp_mooring_ID = mooring ID for selected ADCPs
%       echo_mode       = (logical) Echo Mode status
%       time_shift      = time zone compensation
%       rootDIR         = Raw data directory
%       fRoot           = Raw data file root
%       outRoot         = parent directory of L0 and L1 dircetories
%       dtAvg           = Average sample interval (seconds)
%       atmosphTime     = two times when the instrument was in the air
%       deployTime      = time the instrument was deployed
%       recoverTime     = time the instrument was recovered
% 
%
%
%



clear all
close all
% stages of processing
% 1) define deployment number:
adcp_ID  = 1;
adcp_file_roots = {'S103071A014_KELP1','S104339A001_KELP1'};
adcp_mooring_ID = {'M2', 'M1'};
echo_mode = 0;
% used to shift timezone,e.g. to convert EST to UTC
time_shift = 0/24;
% 2) raw data input directory & filename convention:
rootDIR = sprintf('../../../../Kelp_data/data/FullExperiment/raw/%s',adcp_file_roots{adcp_ID});
fRoot   = [adcp_file_roots{adcp_ID},'_'];
% 3) output directory:
outRoot = '../../../../Kelp_data/Summer2025/Rooker/';
% 4) output data file prefix:
filePrefix= sprintf('ADCP_%s_',adcp_mooring_ID{adcp_ID});
%
% 4a) current/echo average interval (seconds)
dtAvg     = 300;
L0dir     = [outRoot, adcp_mooring_ID{adcp_ID}, '/','L0', '/','ADCP','/'];
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
    offset_file = sprintf([rootDIR,'/',fRoot,'%d.mat'],2);
    load(offset_file,'Config','Data','Descriptions');
    atmosphTime = [datenum('02-Jul-2024 18:45:00'), datenum('02-Jul-2024 19:15:00')];
    it          = find(Data.Burst_Time>=atmosphTime(1) & Data.Burst_Time<=atmosphTime(2));
    ATM_Time    = mean(atmosphTime, 'omitnan');
    ATM_Pressure = mean(Data.Burst_Pressure(it), 'omitnan');
  case 2
    offset_file = sprintf([rootDIR,'/',fRoot,'%d.mat'],1);
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


Config.ATM_Time = ATM_Time;
Config.ATM_Pressure=ATM_Pressure;
Config.Descriptions = Descriptions;
%
% load and pre-process data.
disp('load and pre-process data')
loadSIG1K(Config, rootDIR, fRoot, L0dir, filePrefix, hab, echo_mode, deployTime, recoverTime, HeadingOffset);
%

R1atmTime = [datenum('03-Jul-2024 17:30:00'), datenum('03-Jul-2024 18:10:00')];
R1depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00')];
R23atmTime = [datenum('08-Jul-2024 16:00:00'), datenum('08-Jul-2024 16:30:00')];
R23depTime  = [datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];

% Fetch Release times
disp('fetching data to fit release 1 times')
R1 = fetch_sig1k(L0dir, R1depTime(1), R1depTime(2));

disp('Saving Release 1')
save([outRoot,'/Release1/', adcp_mooring_ID{adcp_ID}, '.mat'],'-struct','R1')

disp('fetching data to fit release 2 times')
R23 = fetch_sig1k(L0dir, R23depTime(1), R23depTime(2));

disp('Saving Release 2')
save([outRoot,'/Release2/', adcp_mooring_ID{adcp_ID}, '.mat'],'-struct','R23')





% make time-averages
%disp('make time-averages')
%time_average_and_rotate_sig1000_RDI_matrix_format_function(L0dir, L0FRoot, filePrefix, dtAvg, echo_mode)
%

% estimate hourly wave stats
%disp('estimate hourly wave stats')
%waves = estimate_wave_bulk_stats_SIG1000_RDI_matrix_format_function(L0dir, filePrefix, L0FRoot, hab);
%
%
% disp('make some pictures')
% fig = figure;
% imagesc(datetime(waves.Time','convertFrom','datenum'),waves.frequency,log10(waves.Spp)),colormap(cmocean('thermal')),caxis([-2 1])
% ylabel('$f$ [Hz]','interpreter','latex')
% cb = colorbar;
% ylabel(cb,'[m$^2$/Hz]','interpreter','latex')
% ax = gca;
% set(ax,'ydir','normal','ticklabelinterpreter','latex','tickdir','out','plotboxaspectratio',[1.5 1 1])
% figname = sprintf('%s/figures/%s_spectra.pdf',L1dir,L1FRoot);
% exportgraphics(fig,figname)






%% Functions

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

for i = 1:length(files)
    
    if (~contains(files(i).name, 'config') && ~contains(files(i).name, 'min'))
        Times = load([files(i).folder, filesep, files(i).name], 'Time');
        if any((Times.Time >= datenum(startTime)) & Times.Time<= datenum(endTime))
            Sig =  load([files(i).folder, filesep, files(i).name]);
        end
    end
end
Data.Time = Sig.Time;
Data.Velocity_East = Sig.Velocity_East';
Data.Velocity_North = Sig.Velocity_North';
Data.Velocity_Up = Sig.Velocity_Up';
Data.Velocity_Beam1 = Sig.Velocity_Beam(:,:,1)';
Data.Velocity_Beam2 = Sig.Velocity_Beam(:,:,2)';
Data.Velocity_Beam3 = Sig.Velocity_Beam(:,:,3)';

Data.Correlation_Minimum = min(Sig.Correlation_Beam(:,:,1)', min(Sig.Correlation_Beam(:,:,2)', Sig.Correlation_Beam(:,:,3)'));
Data.Amplitude_Minimum = min(Sig.Amplitude_Beam(:,:,1)', min(Sig.Amplitude_Beam(:,:,2)', Sig.Amplitude_Beam(:,:,3)'));
Data.Pressure = Sig.Pressure(2,:)';
Data.Heading = Sig.Heading(:, :)';
Data.Pitch = Sig.Pitch(:,:)';
Data.Roll = Sig.Roll(:,:)';

% Summary figure
figure
ax1 = subplot(3, 1, 1);
plot(ax1, Data.Time, Data.Velocity_East, '.')
ylabel({'East', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor

ax2 = subplot(3, 1, 2);
plot(ax2, Data.Time, Data.Velocity_North, '.')
ylabel({'North', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor

ax3 = subplot(3, 1, 3);
plot(ax3, Data.Time, Data.Velocity_Up, '.')
ylabel({'Up', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
grid minor

sgtitle('Sig1K L0 Velocities', 'Fontsize', 25)
end

% EOF


function A = loadSIG1K(Config, rootDIR, fRoot, L0dir, filePrefix, hab, echo_mode, deployTime, recoverTime, HeadingOffset)

% 
% 
%   USAGE: A = load_and_process_sig1000_to_RDI_matrix_format_function(Config, rootDIR, fRoot, L0dir, filePrefix, hab, echo_mode, deployTime, recoverTime, HeadingOffset)
%       Config          = Structure with instrument configuration data
%       rootDIR         = directory for raw sig1K .mat file
%       fRoot           = root of file name for raw sig1K .mat file
%       L0dir           = Directory to save L0 data
%       filePrefix      = L0 file root
%       hab             = Instrument Height Above Bottom (m)
%       echo_mode       = (logical) Echo Mode status (0 == off)
%       deployTime      = time the instrument was deployed
%       recoverTime     = time the instrument was recovered
%       HeadingOffset   = Heading declination in nautical degrees
%       
% 

files = dir([rootDIR,filesep,fRoot,'*.mat']);
Nf    = length(files);
%
%
% average the start/end pressures
ATM_Pressure = mean(Config.ATM_Pressure);
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
%save the config info
outFile = [L0dir, filePrefix, 'config.mat'];
if ~exist(L0dir,'dir')
    eval(['!mkdir ',L0dir])
end
save(outFile,'Config')
% outFile = [L0dir, filePrefix, 'config.nc'];
% if exist(outFile,'file')
%     eval(['!rm ', outFile])
% end
%struct2nc(Config,outFile,'NETCDF4')


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
    out.AST              (1,N+1:N+(ie-(is-1)))     = Data.Burst_AltimeterDistanceAST(is:ie,1);
    out.AST_Offset       (1,N+1:N+(ie-(is-1)))     = Data.Burst_AltimeterTimeOffsetAST(is:ie,1);

    if echo_mode
    out.Echo1            (1:NcE,N+1:N+(ie-(is-1))) = Data.Echo1Bin1_1000kHz_Echo(is:ie,1:NcE)';
    out.Echo2            (1:NcE,N+1:N+(ie-(is-1))) = Data.Echo2Bin1_1000kHz_Echo(is:ie,1:NcE)';
    end
    %
    N = N+ie-(is-1);
    is = ie+1;
    %
    if N==Nmax || (ii==Nf & loadFlag)
        % process data
        minC  = min( out.Correlation_Beam, [], 3);%min( cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        minA  = min( out.Amplitude_Beam  , [], 3);  %cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        maxRNG= out.Pressure(2,:)*cosd(25)-binw;
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
        A = out;
        % save output
         outNf = outNf+1;
         outFileName = sprintf([filePrefix,'%03d.mat'],outNf);
         fout  = [L0dir,outFileName];
         fprintf(['saving output file:   %s \n'],outFileName)
         save(fout,'-struct','out')
%         outFileName = sprintf([filePrefix,'%03d.nc'],outNf);
 %        fout = [L0dir,outFileName];
  %       if exist(fout,'file')
   %         eval(['!rm ', outFileName])
    %     end
%         struct2nc(out,outFileName,'NETCDF4')
%         %
         N = 0;
%         %
         if (ii==Nf & loadFlag )
         fprintf('done! \n')
         break
         end
    end
    %
    %
end


% Summary figure
figure
ax1 = subplot(3, 1, 1);
plot(ax1, A.Time, A.Velocity_East, '.')
ylabel({'East', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor

ax2 = subplot(3, 1, 2);
plot(ax2, A.Time, A.Velocity_North, '.')
ylabel({'North', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor

ax3 = subplot(3, 1, 3);
plot(ax3, A.Time, A.Velocity_Up, '.')
ylabel({'Up', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
grid minor

sgtitle('Sig1K L0 Velocities', 'Fontsize', 25)

end
