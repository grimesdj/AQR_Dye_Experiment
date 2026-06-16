function data = quick_load_rbr(fin)
    % 
    % quick_load_rbr: loads RBR pressure and temp data
    % 
    % USAGE: data = quick_load_rbr(fin)
    % 
    % INPUTS: 
    %             fin = string, file path to .rsk file
    % 
    % OUTPUTS:
    %             data = struct, contains Time, Temp, and Pres (disabled)
    % 
    % written by Derek Grimes and Jason Rooker June 2026


%fin = '../../../Kelp_data/data/FullExperiment/raw/RBR_7_26_24/SN_0953/SN_0953.rsk';
rsk = RSKopen(fin);
rsk = RSKreaddata(rsk);
rbr_time = rsk.data.tstamp;
rbr_data = rsk.data.values;

switch rsk.instruments.model
  case 'RBRsolo'
    if ~strcmp(rsk.channels.longName,'Temperature')
        disp('RBRsolo is not a thermistor')
        keyboard
    end
    rbr_temp = rbr_data;
  case 'RBRduet'
    disp('duet')
    rbr_temp = rbr_data(:,1);
    rbr_pres = rbr_data(:,2);                
end

data.Temp = rbr_temp;
%data.Pres = rbr_pres;
data.Time = rbr_time;

end