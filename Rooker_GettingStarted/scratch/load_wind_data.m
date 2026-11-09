%% wind data

clear all
close all

%% load data
files = dir('../../../../Kelp_data/data/wind/*_wind_data_nearshore.txt');

for j = 1:length(files)
    fname = fullfile(files(j).folder, files(j).name);
    fid = fopen(fname);
    
    hdr1 = fgetl(fid);
    hdr2 = fgetl(fid);
    hdr1 = split(hdr1);
    hdr2 = split(hdr2);
    
    format = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
    A = textscan(fid, format);
    fclose(fid);
    
    for i = 1:length(A)
        varname = matlab.lang.makeValidName(hdr1{i});
        varunit = matlab.lang.makeValidName(hdr2{i});
        var = [varname '_' varunit];
        Data.(var) = A{i};
        Data.(var)(A{i} == 99 | A{i} == 999 | A{i} == 9999) = NaN;
        
        
    end


    % build date time
    str = [string(Data.x_YY_x_yr) + '/' + string(Data.MM_mo) + '/' + string(Data.DD_dy) + ' ' + string(Data.hh_hr) + ':' + string(Data.mm_mn)];
    Data.Time = datetime(str, 'InputFormat','yyyy/MM/dd HH:mm');
    dt = seconds(Data.Time(2)-Data.Time(1));
    fs = 1/dt;
        
    %% fix time gaps
    fields = fieldnames(Data);
    for f = 1:length(fields)
        field = fields{f};
        if ~contains(field, 'Time')
            x = Data.(field);
            t = Data.Time;
            [t_fix, x_fix] = correct_burst(t, x, fs);
            Data.(field) = x_fix;
            Data.Time = datetime(t_fix, 'convertfrom', 'datenum');
        end
    end
    
    %% fliter direction
    X = Data.WSPD_m_s .* -cosd(Data.WDIR_degT);
    Y = Data.WSPD_m_s .* -sind(Data.WDIR_degT);
    
    %% low pass
    wpass = 1/(3*3600);
    
    [y, fsd, idx] = hamming_filter(X, wpass, fs);
    Data.North = y;
    
    [y, fsd, idx] = hamming_filter(Y, wpass, fs);
    Data.East = y;
    
    Data.Time_ds = Data.Time(idx);    
    
    %% Summary Figures?
    
    figure
    ax1 = subplot(4, 1, 1);
    plot(Data.Time, Data.WSPD_m_s, 'k', 'LineWidth', 1)
    hold on
    plot(Data.Time, Data.GST_m_s, 'r--', 'LineWidth', 0.5)
    ylabel({'Wind Speed', '[m/s]'})
    
    ax2 = subplot(4, 1, 2);
    plot(Data.Time, Data.WDIR_degT, 'k', 'LineWidth', 1)
    ylabel({'Wind Dir', '[Deg]'})
    yline(180, 'k--')
    ylim([0 360])
    
    ax3 = subplot(4, 1, 3);
    plot(Data.Time, Data.ATMP_degC, 'k', 'LineWidth', 1)
    ylabel({'Air Temp', '[C]'})
    
    ax4 = subplot(4, 1, 4);
    plot(Data.Time, Data.PRES_hPa, 'k', 'LineWidth', 1)
    ylabel({'Atm Pres', '[hPa]'})
    
    axs = [ax1 ax2 ax3 ax4];
    
    set(axs, 'FontSize', 16)
    linkaxes(axs, 'x')
    xticklabels([ax1 ax2 ax3], [])
    sgtitle(sprintf('%s', files(j).name(1:4)))
    
    %% Save
    fprintf('saving %s...\n', files(j).name)
    save(fname(1:end-3) + "mat", '-struct', 'Data')

end