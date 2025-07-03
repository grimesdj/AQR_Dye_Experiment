clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/data/Release2/raw';
inputFile = 'KELP2_Aquadopp';
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release2/raw';
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = '/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release2/L0';
L0Name  = [inputFile,'_L0'];
% Enter time when instrument was in air for pressure offset
atmTime = [datenum('08-Jul-2024 16:00:00'), datenum('08-Jul-2024 16:30:00')];
depTime = [datenum('08-Jul-2024 17:30:00'), datenum('11-Jul-2024 19:30:00')];
% Enter path to save figures
figDir = [inputDir,'/../figures/'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-offset (UTC->EDT) tos = -4 hours
tos = 0;
%
% returns structure A with all aquadopp data
A = load_AQD_data_function(inputDir, inputFile, fileName, tos);
% Save raw data
save([outputDir,'/',outputName,'.mat'],'-struct','A')

% Generate L0
A = L0_AQD(A, atmTime, depTime);
% Save L0 data
L0 = A.L0;

save([L0Dir,'/',L0Name,'.mat'],'-struct','L0')  

% Create Plots
%L0_plots(A, inputFile, figDir)


