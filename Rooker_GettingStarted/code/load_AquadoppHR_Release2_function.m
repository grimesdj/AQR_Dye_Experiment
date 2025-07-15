
clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '../../../../Kelp_data/data/Release2/raw';
inputFile = 'KELP2_Aquadopp';
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '../../../../Kelp_data/Summer2025/Rooker/Release2/raw';
outputName= [inputFile,'HR_raw'];
outputFile= [outputDir, '/', outputName];
% Enter processed output fileName without .mat
L0Dir   = '../../../../Kelp_data/Summer2025/Rooker/Release2/L0';
L0Name  = [inputFile,'HR_L0'];
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
% Generate and save raw data
disp('Generating raw data')
loadAQD(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime);

% Generate and save L0
disp('Generating L0 data')
L0_AQD(outputFile, L0Dir, L0Name);


% Create Plots
%L0_plots(A, inputFile, figDir)


