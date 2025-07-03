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
% Enter processed output fileName without .mat
L0Dir   = 'C:/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release1/L0';
L0Name  = [inputFile,'_L0'];
% Enter path to save figures
figDir = [inputDir,'/../figures'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-period for estimating the atmospheric pressure offset and deployment
atmosphTime = [datetime('03-Jul-2024 17:30:00'), datetime('03-Jul-2024 18:10:00')];
deployTime  = [datetime('03-Jul-2024 18:30:00'), datetime('03-Jul-2024 22:30:00')];
%
% time offset if necessary
tos = 0;
%
% returns structure A with all vector data
if ~exist([outputDir,'/',outputName,'.mat'],'file')
    A = load_VECTOR_data_function(inputDir, inputFile, fileName, tos);
    save([outputDir,'/',outputName,'.mat'],'-struct','A')
else
    A = load([outputDir,'/',outputName,'.mat']);
    pressure = A.pressure;
    dt = A.seconds(2)-A.seconds(1);
end

A = L0_Vector(A, atmosphTime, deployTime, inputFile, figDir)
L0 = A.L0;

save([L0Dir,'/',L0Name,'.mat'],'-struct','L0')

