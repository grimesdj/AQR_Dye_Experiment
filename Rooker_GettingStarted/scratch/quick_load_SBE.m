function [SBE, flg] = quick_load_SBE(fin)

%fin = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', 'FullExperiment', 'raw', 'SBE56s', 'SN_13546.csv');

%% Open file

fid = fopen(fin);
lines = 0;
check = 0;
flg = 0;

while ~feof(fid)
    line = fgetl(fid);

    % check for .hex file
    if contains(line, '.hex')
        fprintf('Cannot open %s\n', line)
        break
    end

    if contains(line, '"Date","Time","Temperature"')
       check = 1;
        break
    end

    lines = lines + 1;

end
fclose(fid);

% some alternate file types
if ~check
    fprintf('file is invalid...\n')
    flg = 1;
    SBE = [];
    return
end


SBE_in = readtable(fin, "NumHeaderLines", lines, 'ExpectedNumVariables', 3, 'Delimiter', ',');

DTime = SBE_in.Date + SBE_in.Time;
SBE.Time  = datenum(DTime);
SBE.Temp  = SBE_in.Temperature;
