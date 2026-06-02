% Need to download "RSKtools for Matlab" from: https://rbr-global.com/products/software/

fin = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/RBR_7_26_24/SN_0953/SN_0953.rsk';
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
