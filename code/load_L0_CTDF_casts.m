function CTD = load_L0_CTDF_casts(releaseNumber)


% full path
% releaseNumber    = 3;
switch releaseNumber
  case 1
    releaseLatitude  = 34.4693000;
    releaseLongitude = -120.1277833;
    releaseTime      = [datenum('03-Jul-18:38:00'), datenum('03-Jul-19:47:00')];
    ctdSerialNumber  = 234870;
  case 2
    releaseLatitude  = 34.4693000;
    releaseLongitude = -120.1277833;
    releaseTime      = [datenum('08-Jul-18:41:00'), datenum('08-Jul-20:32:00')];
    ctdSerialNumber  = 2348;
  case 3
    releaseLatitude  = 34.4693000;
    releaseLongitude = -120.1277833;
    releaseTime      = [datenum('11-Jul-18:28:00'), datenum('11-Jul-20:55:00')];
    ctdSerialNumber  = 2348;
end
rootDir  = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/'
dataDir  = sprintf('%s/data/Release%d/',rootDir,releaseNumber);
L0Dir    = [dataDir,filesep,'L0',];
fileRoot = sprintf('%d*.mat',ctdSerialNumber);
files    = dir([L0Dir,filesep,fileRoot]);
%
CTD = struct([]);
for jj  = 1:length(files)
    SN  = split(files(jj).name,'_');
    fin= [L0Dir,filesep,SN{1},'_L0.mat'];
    data = load(fin);
    if jj==1
        CTD = data;
    else
        vars = fields(CTD);
        for kk=1:length(vars)
            if ismember(vars{kk},{'pres','cond','dye','dye_raw','temp','time'})
                eval([ 'CTD(1).',vars{kk},'=cat(1,CTD(1).',vars{kk},',data.',vars{kk},');'])
            elseif ismember(vars{kk},{'rawFiles'})
                eval([ 'CTD(1).',vars{kk},'=strvcat(CTD(1).',vars{kk},',data.',vars{kk},');'])
            elseif ismember(vars{kk},{'pres_grid'})
                continue
            else
                eval([ 'CTD(1).',vars{kk},'=cat(2,CTD(1).',vars{kk},',data.',vars{kk},');'])
            end
        end
    end
end
