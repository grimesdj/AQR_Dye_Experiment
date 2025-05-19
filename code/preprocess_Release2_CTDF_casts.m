clear all
close all
% full path
releaseNumber  = 2;
%ctdSerialNumber=234870;
rootDir  = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release%d/',releaseNumber);
rawDir    = [rootDir,filesep,'raw'];
L0Dir     = [rootDir,filesep,'L0'];
fileRoot = '2348*.rsk';
files    = dir([rawDir,filesep,fileRoot]);
nf       = length(files);
ATM_Time = [datenum('08-Jul-2024 18:04:08'), datenum('08-Jul-2024 18:09:07');...
            datenum('08-Jul-2024 20:08:21'), datenum('08-Jul-2024 20:44:48')];
%
%
for ii=1:nf
fin = [rawDir,filesep,files(ii).name];
SN  = split(files(ii).name,'_');
fout= [L0Dir,filesep,SN{1},'_L0'];
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
    figure, plot(time,pres_abs,'*')
    title({'zoom into region before deployment'; 'press \CR (return/enter) when finished'})
    pause()
    title({'select times bounding when instrument was in the air'})
    inAir = ginput(2);
    datestr(inAir(:,1))
else
    inAir = ATM_Time(ii,:)';
end
iAir  = find(time>=inAir(1,1) & time<=inAir(2,1));
pres0 = nanmean(pres_abs(iAir))
pres  = pres_abs-pres0;
%
% get the GPS data
[gps,rawFiles] = get_CTDF_gps_data(releaseNumber,SN);
rawFiles = strvcat(fin,char(rawFiles));
%
% interpolate downcasts to uniform grid
[time_grid,latitude_grid,longitude_grid,pres_grid,dye_grid,temp_grid,cond_grid] = interpolate_CTDF_downcasts(time,pres,temp,cond,dye,gps);
%
%
% save gridded data to netcdf files
CTD = struct('Time',time_grid,'Latitude',latitude_grid,'Longitude',longitude_grid,'Pressure',pres_grid,'Dye',dye_grid,'Temperature',temp_grid,'Conductivity',cond_grid)
%
% archive raw and gridded data
save(fout,'time','cond','temp','pres','dye_raw','dye','rawFiles','gps','time_grid','pres_grid','temp_grid','cond_grid','dye_grid','latitude_grid','longitude_grid','pres0')
%
ncfile = [fout,'.nc'];
if exist(ncfile,'file')
    eval(['!rm ',ncfile])
end
struct2nc(CTD,ncfile,'NETCDF4');
%
%
close(gcf)
end
% $$$ % 2) get continuous record of GPS based on dye release notes on which GPS to use
% $$$ %    see ./info/notes/IMG_1925.heic
% $$$ gps_tmp = struct([]);
% $$$ vars    = fields(gps);
% $$$ files   = {};
% $$$ if ismember('234870',SN)
% $$$     gpsIDs = {'GPSUSER06','GPSUSER05'};
% $$$     tchng  = datenum('07/08/2024 20:08:00');
% $$$     for jj=1:length(gpsIDs)
% $$$         for kk=1:length(gps)
% $$$             gpsFile = gps(kk).file
% $$$             gpsSplit = split(gpsFile,{'/','_'});
% $$$             if ismember(gpsIDs(jj),gpsSplit)
% $$$                 switch jj
% $$$                   case 1
% $$$                     inds = (gps(jj).time<=tchng);
% $$$                     files = {files; gps(kk).file};
% $$$                   case 2
% $$$                     inds = (gps.time>=tchng);
% $$$                     files = {files; gps(kk).file};
% $$$                 end
% $$$                 for ll = 1:length(vars)
% $$$                     var = vars(ll);
% $$$                     if ll==1
% $$$                         gps_tmp.file=files;
% $$$                         continue
% $$$                     end
% $$$                     eval(['gps_tmp.',var,'=cat(1,gps_tmp(1).',var,',gps(kk).',var,');'])
% $$$                 end
% $$$             end
% $$$         end
% $$$     end
% $$$     %
% $$$     [~,srt] = sort(gps_tmp.time);
% $$$     for ll = 2:length(vars)
% $$$         var = vars(ll);
% $$$         eval(['gps_tmp.',var,'=gps_tmp.',var,'(srt);'])
% $$$     end
% $$$     return
% $$$ end
