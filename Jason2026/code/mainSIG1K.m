
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
adcp_ID  = 2;
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
save([outRoot,'/Release1/', adcp_mooring_ID, '.mat'],'-struct','R1')

disp('fetching data to fit release 2 times')
R23 = fetch_sig1k(L0dir, R23depTime(1), R23depTime(2));

disp('Saving Release 2')
save([outRoot,'/Release2/', adcp_mooring_ID, '.mat'],'-struct','R23')





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
end

% EOF



