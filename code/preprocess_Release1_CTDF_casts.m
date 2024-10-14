clear all
close all
% full path
releaseNumber  = 1;
ATM_Time = [datenum('03-Jul-2024 20:47:27'), datenum('03-Jul-2024 20:59:48')]
ctdSerialNumber=234870;
rootDir  = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release%d/',releaseNumber);
rawDir    = [rootDir,filesep,'raw'];
L0Dir    = [rootDir,filesep,'L0'];
fileRoot = [num2str(ctdSerialNumber),'*.rsk'];
%
%
files    = dir([rawDir,filesep,fileRoot]);
fin = [rawDir,filesep,files(1).name];
SN  = split(files(1).name,'_');
fout= [L0Dir,filesep,SN{1},'_L0.mat'];
SN  = str2num(SN{1});
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
dye      = convert_fluorometer_counts_to_ppb(SN,dye_raw,temp);
%
%
if ~exist('ATM_Time','var')
figure, plot(time,pres_abs)
title({'zoom into region before deployment'; 'press \CR (return/enter) when finished'})
pause()
title({'select times bounding when instrument was in the air'})
inAir = ginput(2);
    datestr(inAir(:,1))
else
    inAir = ATM_Time(1,:)';
end
iAir  = find(time>=inAir(1,1) & time<=inAir(2,1));
pres0 = mean(pres_abs(iAir))
pres  = pres_abs-pres0;
%
%
% add the gps data to file
Release1GPS = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release1/release1_ctdf_cast_gps.csv';
format = '%*s %s %f %f';
fid    = fopen(Release1GPS);
data   = textscan(fid,format,'HeaderLines',1,'delimiter',{','});
gps_time = datenum(data{1});
lat      = data{2};
lon      = data{3};
%
gps_loc  = lat+sqrt(-1)*lon;
ctd_loc  = interp1(gps_time,gps_loc,time);
%
gps    = struct([]);
gps(1).file = Release1GPS;
gps(1).time = time;
gps(1).latitude = real(ctd_loc);
gps(1).longitude= imag(ctd_loc);
rawFiles        = strvcat(fin,Release1GPS);
%
% interpolate downcasts to uniform grid
[time_grid,latitude_grid,longitude_grid,pres_grid,dye_grid,temp_grid,cond_grid] = interpolate_CTDF_downcasts(time,pres,temp,cond,dye,gps);
%
% archive raw and gridded data
save(fout,'time','cond','temp','pres','dye','rawFiles','gps','time_grid','pres_grid','temp_grid','cond_grid','dye_grid','latitude_grid','longitude_grid','pres0')
%
% save gridded data to netcdf files
CTD = struct('Time',time_grid,'Latitude',latitude_grid,'Longitude',longitude_grid,'Pressure',pres_grid,'Dye',dye_grid,'Temperature',temp_grid,'Conductivity',cond_grid)
ncfile = [fout,'.nc'];
if exist(ncfile,'file')
    eval(['!rm ',ncfile])
end
struct2nc(CTD,ncfile,'NETCDF4');
%
%
close(gcf)
