clear all
close all
% full path
releaseNumber  = 3;
%ctdSerialNumber=234870;
rootDir  = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release%d/',releaseNumber);
rawDir    = [rootDir,filesep,'raw'];
L0Dir    = [rootDir,filesep,'L0'];
fileRoot = '*.rsk';
files    = dir([rawDir,filesep,fileRoot]);
nf       = length(files);
for ii=1:nf
fin = [rawDir,filesep,files(ii).name];
SN  = split(files(ii).name,'_');
fout= [L0Dir,filesep,SN{1},'_L0.mat'];
% open the file
rsk = RSKopen(fin);
%
rsk = RSKreaddata(rsk);
%  mat = RSK2MAT(rsk);
serialNum = rsk.instruments.serialID;
% $$$ time      = rsk.data.tstamp;
% $$$ pres      = rsk.data.values;
time = rsk.data.tstamp;
data = rsk.data.values;
%
%
cond     = data(:,1);
temp     = data(:,2);
pres_abs = data(:,3);
dye_raw  = data(:,4);
%
%
figure, plot(time,pres_abs,'*')
title({'zoom into region before deployment'; 'press \CR (return/enter) when finished'})
pause()
title({'select times bounding when instrument was in the air'})
inAir = ginput(2);
%
%
iAir  = find(time>=inAir(1,1) & time<=inAir(2,1));
pres0 = nanmean(pres_abs(iAir))
pres  = pres_abs-pres0;
%
%
% $$$ % add the gps data to file
% $$$ Release1GPS = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release1/release1_ctdf_cast_gps.csv';
% $$$ format = '%*s %s %f %f';
% $$$ fid    = fopen(Release1GPS);
% $$$ data   = textscan(fid,format,'HeaderLines',1,'delimiter',{','});
% $$$ time_gps = datenum(data{1})-3/24;% for some reason these are 3-hours off
% $$$ lat      = data{2};
% $$$ lon      = data{3};
% $$$ % save(fout,'-append')
%
% break record into individual casts
% 1) first smooth pressure record to find zero-up crossings
tau = 3;% smooth timescale in seconds
dt  = 1/8;% sample rate Hz
N   = tau/dt;
f   = hanning(N); f = f./sum(f);
pres_avg = conv(pres,f,'same');
% $$$ %
% $$$ figure,
% $$$ xline(datetime(time_gps,'convertfrom','datenum'),'--r'), hold on
% $$$ plot(datetime(time,'convertfrom','datenum'),pres_avg)
% $$$ %
% remove data from surface
oowFlag = (pres<0.02);
%
cond(oowFlag)=nan;
temp(oowFlag)=nan;
pres(oowFlag)=nan;
dye_raw(oowFlag)=nan;
%
% $$$ save(fout,'time','cond','temp','pres','dye_raw','time_gps','lat','lon')
save(fout,'time','cond','temp','pres','dye_raw')
%
close(gcf)
end