
clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '../../../../Kelp_data/data/Release1/raw';
inputFile = 'KELP1_AquadoppHR';
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '../../../../Kelp_data/Summer2025/Rooker/Release1/raw';
outputName= [inputFile,'_raw'];
outputFile= [outputDir, '/', outputName];
% Enter processed output fileName without .mat
L0Dir   = '../../../../Kelp_data/Summer2025/Rooker/Release1/L0';
L0Name  = [inputFile,'_L0'];
% Enter time when instrument was in air for pressure offset
atmTime = [datenum('03-Jul-2024 14:00:00'), datenum('03-Jul-2024 18:00:00')];
depTime = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00')];
% Enter path to save figures
figDir = [outputDir,'/../figures/'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-offset (UTC->EDT) tos = -4 hours
tos = 0;
%
% Generate and save raw data
disp('Generating raw data')
%loadAQD(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime);

% Generate and save L0
disp('Generating L0 data')
%L0_AQD(outputFile, L0Dir, L0Name);

% Rotate to Principas axes
disp('Applying PCA rotation')
pca_function(L0Dir, L0Name)



% Create Plots
disp('Generating figures')
L0_plots(L0Dir, L0Name, figDir, inputFile)


fprintf('\n====================\n\nDone!\n\n====================\n')

