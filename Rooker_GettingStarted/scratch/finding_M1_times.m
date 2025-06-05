



files = dir("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\2024_PROCESSED_DATA\M1\L0\ADCP\ADCP_M1_*.mat");

TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");


ADCP = struct();
M1 = struct();
for i = 1:length(files)
    ind = sprintf('M1_%d', i);
    if ~contains(files(i).name, 'config')  &&  ~contains(files(i).name, 'min')
        ADCP.(ind) = load([files(i).folder,filesep,files(i).name],'Time');
        for  releaseNum = 1:3
            Num = sprintf('release%d', releaseNum);
            time.(Num) = datenum(TRange.StartTime_UTC_(releaseNum));
            if contains(string(datenum(ADCP.(ind).Time)), string(time.(Num)))
                M1.(releaseNum) = load([files(i).folder,filesep,files(i).name]);
            end
            

        end
        
    
    end
    
end

