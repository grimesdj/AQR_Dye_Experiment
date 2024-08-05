clear all
close all
% full path
releaseNumber  = 2;
%ctdSerialNumber=234870;
rootDir  = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release%d/',releaseNumber);
rawDir    = [rootDir,filesep,'raw'];
L0Dir     = [rootDir,filesep,'L0'];
fileRoot = '*.rsk';
files    = dir([rawDir,filesep,fileRoot]);
nf       = length(files);
%
% get gps info
gpsInfoFile = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/info/release%d_ctdf_cast_gps/gps_ctd_serial_number_info.csv',releaseNumber);
gpsFID = fopen(gpsInfoFile,'r');
gpsDAT = textscan(gpsFID,'%s %d %s','Headerlines',1,'delimiter',',');
%
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
% Pair CTDF and GPS data
rawFiles = {fin};
isPair = ismember(gpsDAT{2},str2num(SN{1}));
indPair= find(isPair);
% for each pair, get the GPS data file and read time/lat/lon/etc.
gps = struct([]);
for jj=1:sum(isPair)
    string  = char(gpsDAT{3}(indPair(jj)));
    gpsRoot = sprintf([rawDir,filesep,'%s','*.TXT'],string);
    gpsFiles= dir(gpsRoot);
    gpsPath = [gpsFiles.folder,filesep,gpsFiles.name];
    rawFiles{jj+1}=gpsPath;
    % get date from filename
    gpsStr = split(gpsFiles.name,{'_','.'});
    gpsDate= gpsStr{3}();
    gpsYear= str2num(gpsDate(1:4)); gpsMon = str2num(gpsDate(5:6)); gpsDay = str2num(gpsDate(7:8));
    gpsTime= gpsStr{4}();
    gpsHH  = str2num(gpsTime(1:2)); gpsMM = str2num(gpsTime(3:4)); gpsSS = str2num(gpsTime(5:6));
    %
    [gps_hour, gps_minute, gps_second, velocity, numsats, height, lat, lon, quality, geoidsep, pdop, hdop, vdop] = readGPS(gpsPath);
    N_gps = length(gps_hour);
    time_gps = datenum([repmat(gpsYear,N_gps,1),repmat(gpsMon,N_gps,1),repmat(gpsDay,N_gps,1),gps_hour,gps_minute,gps_second]);
    % add a day if necessary
    dTime = diff(time_gps);
    indDayChng = find(dTime<-86398/86400)+1;
    for kk = 1:length(indDayChng)
        time_gps(indDayChng(kk):end)=time_gps(indDayChng(kk):end)+1;
    end
    gps(jj).file=gpsPath;
    gps(jj).time=time_gps;
    gps(jj).velocity=velocity;
    gps(jj).latitude=lat;
    gps(jj).longitude=lon;
    gps(jj).height = height;
    gps(jj).quality=quality;
    gps(jj).geoidsep=geoidsep;
    gps(jj).pdop=pdop;
    gps(jj).hdop=hdop;
    gps(jj).vdop=vdop;
end
%
%
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
%
%
% break record into individual casts
% 1) first smooth pressure record to find zero-up crossings
tau = 5;% smooth timescale in seconds
dt  = 1/8;% sample rate Hz
N   = tau/dt;
f   = hanning(N); f = f./sum(f);
pres_avg = conv(pres,f,'same');
dpres_avg = gradient(pres_avg);
% locate the 'casts'
cast          = pres_avg>=0.5;
downcast      = find(cast & dpres_avg>0);
upcast        = find(cast & dpres_avg<0);
cast          = find(cast);
%
[startINDs,endINDs] = Segment(cast,(1/dt));
% $$$ [~, endINDs]   = Segment(upcast  ,(tau/dt));
%
cast_start     = cast(startINDs);
cast_end       = cast(endINDs);
Ncasts         = length(cast_start);
%
% now create a uniform pressure grid with 10cm resolution to the maximum depth (15 m)
pres_grid = [0:0.1:15]';
Npres     = length(pres_grid);
time_grid = nan(1,Ncasts);
latitude_grid  = nan(1,Ncasts);
longitude_grid = nan(1,Ncasts);
temp_grid = nan(Npres,Ncasts);
dye_grid  = nan(Npres,Ncasts);
cond_grid = nan(Npres,Ncasts);
for jj = 1:Ncasts
    this_cast = find(downcast>=cast_start(jj) & downcast<=cast_end(jj));
    [~,inds]  = unique(pres(downcast(this_cast)));
    time_grid(1,jj) = time(downcast(this_cast(inds(1))));
    latitude_grid(1,jj)  = interp1(gps(1).time,gps(1).latitude,time_grid(1,jj));
    longitude_grid(1,jj) = interp1(gps(1).time,gps(1).longitude,time_grid(1,jj));
    temp_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),temp(downcast(this_cast(inds))),pres_grid);
    cond_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),cond(downcast(this_cast(inds))),pres_grid);
    dye_grid(:,jj)  = interp1(pres(downcast(this_cast(inds))),dye_raw(downcast(this_cast(inds))),pres_grid);
end
%
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
save(fout,'time','cond','temp','pres','dye_raw','rawFiles','gps','time_grid','pres_grid','temp_grid','cond_grid','dye_grid','latitude_grid','longitude_grid')
%
% save to netcdf files

%
close(gcf)
end