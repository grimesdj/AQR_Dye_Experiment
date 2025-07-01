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
rootDIR = sprintf('/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/data/FullExperiment/raw/%s',adcp_file_roots{adcp_ID});
fRoot   = [adcp_file_roots{adcp_ID},'_'];
% 3) output directory:
outRoot = '/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/';
% 4) output data file prefix:
filePrefix= sprintf('ADCP_%s_',adcp_mooring_ID{adcp_ID});
%
% 4a) current/echo average interval (seconds)
dtAvg     = 300;
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
    offset_file = sprintf([rootDIR,'/',fRoot,'%d.mat'],2);
    load(offset_file,'Config','Data','Descriptions');
    atmosphTime = [datenum('02-Jul-2024 18:45:00'), datenum('02-Jul-2024 19:15:00')];
    it          = find(Data.Burst_Time>=atmosphTime(1) & Data.Burst_Time<=atmosphTime(2));
    ATM_Time    = nanmean(atmosphTime);
    ATM_Pressure = nanmean(Data.Burst_Pressure(it));
  case 2
    offset_file = sprintf([rootDIR,filesep,fRoot,'%d.mat'],1);
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
%
% load and pre-process data.
load_and_process_sig1000_to_RDI_matrix_format_function(Config, rootDIR, fRoot, L0dir, filePrefix, ATM_Time, ATM_Pressure, hab, echo_mode, 'Descriptions')
%
return
% make time-averages
time_average_and_rotate_sig1000_RDI_matrix_format
%
% estimate hourly wave stats
estimate_wave_bulk_stats_SIG1000_RDI_matrix_format
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



