%
whoseIMAGES   = 'grimes_images';
releaseNumber = 1;
rootDir = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote';
infoDir = [rootDir,filesep,'info'];
dataDir = sprintf([rootDir,filesep,'data',filesep,'Release%d'],releaseNumber);
imageRoot = sprintf('release%d_ctdf_cast_gps',releaseNumber);
releaseDir= [infoDir,filesep,imageRoot,filesep,whoseIMAGES];
inPath    = [releaseDir,filesep,'*.HEIC'];
outPath   = sprintf([dataDir,filesep,imageRoot,'.csv']);
%
files = dir(inPath);
nf    = length(files);
castTimes = char();
for jj = 1:nf
    castTimes(jj,:) = files(jj).date;
    if releaseNumber==1 & jj==1
        castTimes(jj,1:2)='04';
        castTimes(jj,13:14)='00';
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

