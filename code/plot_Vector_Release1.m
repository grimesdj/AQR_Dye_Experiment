clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/Release1/raw';
inputFile = 'KELP1_Vector';
headingOffset = 331.4;% based on heading from AquodoppHR_KELP1
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/Release1/raw';
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/Release1/L0';
L0Name  = [inputFile,'_L0'];
% Enter path to save figures
figDir = [inputDir,'/../figures/'];
%
% Enter time-period for estimating the atmospheric pressure offset and deployment
atmosphTime = [datetime('03-Jul-2024 17:30:00'), datetime('03-Jul-2024 18:10:00')];
deployTime  = [datetime('03-Jul-2024 18:30:00'), datetime('03-Jul-2024 22:30:00')];
% load([outputDir,'/',outputName,'.mat'])
load([L0Dir,'/',L0Name,'.mat'])
%
% dye release time...
info = get_release_info(1);
%
figure,
ax1 = subplot(2,1,1);
plot(b1,'.')
ax2 = subplot(2,1,2);
plot(pressure,'-')
linkaxes([ax1, ax2],'x')
set(ax2,'xlim',[3.015e5 3.06e5],'ylim',[-0.01 10])


figure,
ax1 = subplot(2,1,1);
plot(a1,'.')
ax2 = subplot(2,1,2);
plot(b1,'.')
linkaxes([ax1, ax2],'x')
set(ax2,'xlim',[3.015e5 3.06e5])
