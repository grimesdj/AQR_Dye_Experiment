clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = 'C:/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/data/Release1/raw';
inputFile = 'KELP1_Vector';
headingOffset = 331.4;% based on heading from AquodoppHR_KELP1
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = 'C:/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release1/raw';
outputName= [inputFile,'_raw'];
outputFile= [outputDir, '/', outputName];
% Enter processed output fileName without .mat
L0Dir   = 'C:/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release1/L0';
L0Name  = [inputFile,'_L0'];
% Enter path to save figures
figDir = [inputDir,'/../figures'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-period for estimating the atmospheric pressure offset and deployment
atmTime = [datenum('03-Jul-2024 17:30:00'), datenum('03-Jul-2024 18:10:00')];
depTime  = [datenum('03-Jul-2024 18:30:00'), datenum('03-Jul-2024 22:30:00')];
%
% time offset if necessary
tos = 0;
%
% returns structure A with all vector data
%if ~exist([outputDir,'/',outputName,'.mat'],'file')
 load_VECTOR_data_function(inputDir, inputFile, fileName, outputFile, tos, depTime, atmTime);
    
%else
%     A = load([outputDir,'/',outputName,'.mat']);
%     pressure = A.Pressure;
%     dt = A.Seconds(2)-A.Seconds(1);
% %end
L0_Vector(outputFile, L0Dir, L0Name);



