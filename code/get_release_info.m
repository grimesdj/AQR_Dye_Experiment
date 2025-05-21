function info = get_release_info(release_num);
%
% USAGE: info = get_release_info(release_num);
%
% input release number and return:
% info = 
% struct with fields:
%    concentration: 0.0850     [kg-dye/kg-soln]
%    released_mass: 29         [kg-soln]
%       start_time: 7.3944e+05 [UTC-Julian Date]
%         end_time: 7.3944e+05 [UTC-Julian Date]


releaseFile = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/info/dye_mixing_cals_and_releases/dye_release_times.csv';

% Release Number, Dye Concentration [kg-dye/kg-soln], Weight Released (kg), StartTime (UTC), EndTime (UTC)
format = '%d %f %f %s %s';
fid    = fopen(releaseFile);
data   = textscan(fid,format,'delimiter',',','headerlines',1);

info.concentration = data{2}(release_num);
info.released_mass = data{3}(release_num);
info.start_time    = datenum(data{4}{release_num});
info.end_time      = datenum(data{5}{release_num});