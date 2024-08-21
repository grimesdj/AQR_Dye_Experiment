%
whoseIMAGES   = 'centeno_images/latest_images/';%'grimes_images';%;
releaseNumber = 1;
rootDir = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote';
infoDir = [rootDir,filesep,'info'];
dataDir = sprintf([rootDir,filesep,'data',filesep,'Release%d'],releaseNumber);
imageRoot = sprintf('release%d_ctdf_cast_gps',releaseNumber);
releaseDir= [infoDir,filesep,imageRoot,filesep,whoseIMAGES];
inPath    = [releaseDir,filesep];
ext       = {'*.HEIC','*.JPG'};
extensions= cellfun(@(x)dir(fullfile(inPath,x)),ext,'UniformOutput',false);
outPath   = sprintf([dataDir,filesep,imageRoot,'_image_file_times.csv']);
%
% files = dir(inPath,extensions);
files = vertcat(extensions{:});
nf    = length(files);
castTimes = char();
for jj = 1:nf
    castTimes(jj,:) = files(jj).date;
    % the first image in grimes' batch had different time-zone than others
    if releaseNumber==1 & jj==1 & strcmp(whoseIMAGES,'grimes_images')
        castTimes(jj,1:2)='04';
        castTimes(jj,13:14)='00';
    end
    %
    % all other images in grimes' batch are off by 3-hours
    if strcmp(whoseIMAGES,'grimes_images')
        imageTime = datenum(castTimes(jj,:))-3/24;
        castTimes(jj,:) = datestr(imageTime);
    end
    %
    % all of centeno's batch is off by a few hours, too. Maybe the RBR was in wrong time-zone?
    if strcmp(whoseIMAGES,'centeno_images/latest_images/')
        imageTime = datenum(castTimes(jj,:))+4/24;
        castTimes(jj,:) = datestr(imageTime);
    end
end
castTimes = cellstr(strcat('''',castTimes,''''));
imageNames= {files.name}';
header    = {'Image names', 'GPS-cast times'};
data      = cell2table([imageNames,castTimes],'variablenames',header);
if ~exist(outPath,'file')
    writetable(data,outPath)
else
    question = '\n csv-file: %s already exists, overwrite? (1=yes,0=no) \n';
    response = input(question);
    switch response
      case 1
        writetable(data,outPath)
      case 0
        writetable(data,outPath,'WriteMode','append')
        fprintf('appending image names/times to existing file!\n')
    end
end

