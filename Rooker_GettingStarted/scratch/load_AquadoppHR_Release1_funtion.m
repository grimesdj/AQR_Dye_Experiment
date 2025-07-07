
clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '../../../../Kelp_data/data/Release1/raw';
inputFile = 'KELP1_AquadoppHR';
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '../../../../Kelp_data/Summer2025/Rooker/Release1/raw';
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = '../../../../Kelp_data/Summer2025/Rooker/Release1/L0';
L0Name  = [inputFile,'_L0'];
% Enter time when instrument was in air for pressure offset
atmTime = [datenum('03-Jul-2024 14:00:00'), datenum('03-Jul-2024 18:00:00')];
depTime = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00')];
% Enter path to save figures
figDir = [inputDir,'/../figures/'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-offset (UTC->EDT) tos = -4 hours
tos = 0;
%
% returns structure A with all aquadopp data
disp('Generating raw data')
A = loadAQD(inputDir, inputFile, fileName, tos, depTime, atmTime);
% Save raw data
disp('Saving raw data')
save([outputDir,'/',outputName,'.mat'],'-struct','A')

%aquawrap(A.Velocity_Beam1, A.VRange)

% Generate L0
disp('Generating L0 data')
A = L0_AQD(A, atmTime, depTime);
% Save L0 data
L0 = A.L0;

save([L0Dir,'/',L0Name,'.mat'],'-struct','L0')

% Create Plots
%L0_plots(A, inputFile, figDir)


