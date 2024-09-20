function dye = convert_fluorometer_counts_to_ppb(SN0,dye0,temp0,dye_cal_file);
%
% usage: dye = convert(SN0,dye0,temp0,dye_cal_file);
%
% convert from fluorescent counts to ppb using calibration data in dye_cal_file.
% equation is: dye = slope*(dye0-offset),
% with an exponential correction to the slope: slope = slope*exp(0.027*(temp0-temp_cal))

% load calibration data
format = '%f %f %f %f';
fid    = fopen(dye_cal_file);
header = fgetl(fid);
header_vars = split(header,',');

cal_data = textscan(fid,format,'delimiter',',');
% parse variables
SN     = cal_data{1};
offset = cal_data{2};
slope  = cal_data{3};
temp   = cal_data{4};

ID     = ismember(SN,SN0);

% convert based on slope/offset
dye  = slope(ID).*(dye0-offset(ID));

% correct for thermal effects
if ~isempty(temp0)
    dye  = dye.*exp(0.027*(temp0-temp(ID)));
end

% zero-out negative dye
dye(dye<=0)=0;

% note: I'm not correcting for turbidity!