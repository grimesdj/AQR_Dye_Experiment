clear all
close all
%
% 1) load the calibration info (serial numbers, times, thermistors)
infoRoot = '~/OneDriveUNCW/KELP-vingilote/info/dye_mixing_cals_and_releases/';
fluoRoot = '~/OneDriveUNCW/KELP-vingilote/data/dye_calibrations/';
thrmRoot = '~/OneDriveUNCW/KELP-vingilote/data/FullExperiment/SBE56s/';
figRoot  = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/figures/calibrations/';
%
calibration_times_file = [infoRoot,'dye_bucket_times_vs_serial_number.csv'];
thermistor_info_file   = [infoRoot,'dye_bucket_thermistor_serial_numbers.csv'];
fluorometer_info_file  = [infoRoot,'instrument_serial_number_type.csv'];
%
% 1a) load the calibration & thermistor info:
calib_format = ['%d ',repmat('%s ',1,15),' %*s'];
fid          = fopen(calibration_times_file,'r');
calib_header = fgetl(fid);
calib_ppb    = textscan(calib_header,['%*s ',repmat('%f ',1,15),'%*s'],'delimiter',',');
calib_ppb    = cell2mat(calib_ppb(:));
calib_info   = textscan(fid, calib_format,'delimiter',',');
bucket_SN    = calib_info{1};
bucket_str   = calib_info(2:16);
%
% important note: times were logged as HH:MM:SS UTC,
% we started on 06/30 15:30:00 and finished at 07/01 04:45
Ninst = size(bucket_SN,1);
Nbckt = size(bucket_str,2);
bucket_times = nan(Ninst,Nbckt);
cal_start_day_str = '30 June 2024 ';
cal_start_time    = datenum([cal_start_day_str,'15:30:00']);
for i = 1:Ninst
    for j = 1:Nbckt
        bucket_times(i,j) = datenum([cal_start_day_str, char(bucket_str{j}(i,:))]);
        if bucket_times(i,j)<cal_start_time
            bucket_times(i,j) = bucket_times(i,j)+1;
        end
    end
end
%
therm_format = '%f %d';
fid          = fopen(thermistor_info_file,'r');
therm_info   = textscan(fid, therm_format,'delimiter',',','HeaderLines',1);
therm_ppb    = therm_info{1};
therm_SN     = therm_info{2};
therm_time   = {};
therm_temp   = {};
for j = 1:Nbckt
    thrmFile = dir(sprintf('%sSN_%d.csv',thrmRoot,therm_SN(1)));
    fid      = fopen([thrmFile.folder,filesep,thrmFile.name]);
    hdrflag  = 0;
    while ~hdrflag
        line = fgetl(fid);
        strings = split(line,{',','"'});
        if ismember('Date',strings);
            hdrflag = 1;
        end
    end
    therm_data = textscan(fid,'%q %q %q','delimiter',',');
    Nthrm = size(therm_data{1},1);
    therm_time{j} = datenum([char(therm_data{1}),repmat('-',Nthrm,1),char(therm_data{2})],'yyyy-mm-dd-HH:MM:SS');
    therm_temp{j} = cellfun(@str2num,therm_data{3});
end
%
% 2) for each serial number: determine the instrument type, load the calibration data,...
% 2a) get list of SN / instrument type
fid = fopen(fluorometer_info_file);
fluoro_info = textscan(fid,'%d %s','delimiter',',','HeaderLines',1);
fluoro_SN   = fluoro_info{1};
fluoro_type = fluoro_info{2};
%
% 2b) get list of cal files
cal_files   = dir(fluoRoot);
cal_file_names = {cal_files.name}';
%
% 2c) loop over instruments, pair instrument SN to it's cal file, load, and segment
for i = 1:Ninst
    % get current SN
    SN = bucket_SN(i);
    iFL     = find(ismember(fluoro_SN,SN));
    FL_type = fluoro_type{iFL};
    % now load the fluorometer data
    switch FL_type
      case {'ECO','ECO-BIO','ECO-T'}
        format = '%s %s %f %f %f %f %f %f %f';% !!!!stopped here!!!!
        if all(strcmp(FL_type,'ECO-T'))
            format = cat(2,format,' %f');
        end
        fin = sprintf([fluoRoot,filesep,'%d.raw'],SN);
        fid = fopen(fin,'r');
        str = fgetl(fid);
        N   = str2num(str(1:strfind(str,'records')-1));
        dat = textscan(fid,format,'delimiter','\t');
        mmddyy = dat{1}(1:N);
        HHMMSS = dat{2}(1:N);
        time = datenum( strcat(mmddyy, repmat('-',N,1), HHMMSS), 'mm/dd/yy-HH:MM:SS');
        NTU  = dat{4}(1:N);
        CHL  = dat{6}(1:N);
        RWT  = dat{8}(1:N);
        flag = (RWT>4096 | CHL>4096 | NTU>4096);
        NTU(flag)=nan; CHL(flag)=nan; RWT(flag)=nan;
        % despike raw data
        RWT = despike_fluoro(RWT,1);
        % how long was instrument in bucket? [decimal day]
        dt   = 30/86400;
      case 'PME'
        format = '%f %f %f %f %d';
        fin = sprintf([fluoRoot,filesep,'PME_%d.txt'],SN);
        fid = fopen(fin,'r');
        dat = textscan(fid,format,'delimiter',',','headerlines',3);
        SS  = dat{1};
        time= datenum('Jan 01 1970')+SS/86400;
        RWT = dat{4};
        dt  = 300/86400;
        if SN==339378
            time = time + (datenum('30-Jun-2024 22:45:00') - datenum('28-Jun-2024 05:45:00'));
        end
      case 'YSI'
        format = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
        fin = sprintf([fluoRoot,filesep,'SN_%04d.csv'],SN);
        fid = fopen(fin,'r');
        dat = textscan(fid,format,'delimiter',',','headerlines',10);
        N   = size(dat{1},1);
        time = datenum( strcat(dat{2},repmat('-',N,1),dat{1}) );
        RWT  = dat{15};
        dt   = 30/86400;
      case {'RBR-C-F','RBR-V-F','RBR-CTD-F'}
        RBRfile = dir(sprintf([fluoRoot,filesep,'%d_*.rsk'],SN));
        fin     = [RBRfile.folder,filesep,RBRfile.name];
        % open the file, then read data
        rsk = RSKopen(fin);
        rsk = RSKreaddata(rsk);
        time = rsk.data.tstamp;
        data = rsk.data.values;            
            if (strcmp(FL_type,'RBR-V-F') | strcmp(FL_type,'RBR-C-F'))
                RWT  = data(:,1);
            elseif (strcmp(FL_type,'RBR-CTD-F'))
                RWT  = data(:,4);
            end
        RWT = despike_fluoro(RWT,1);
        dt  = 30/86400;
    end
    % 3) segment into 30 second or 5 minute bucket time-series, apply T-correction dye = dye_raw .* exp(0.027*(temp_avg - Tcal));
    for j = 1:Nbckt
        bucket_start_time = bucket_times(i,j);
        fl_in_bucket      = find(time>=bucket_start_time & ...
                                 time<=bucket_start_time+dt);
        thrm_in_bucket    = find(therm_time{j}>=bucket_start_time & ...
                                 therm_time{j}<=bucket_start_time+dt);
        bktRWTavg(j) = nanmean(RWT(fl_in_bucket));
        bktRWTstd(j) = nanstd (RWT(fl_in_bucket));
        bktTEMP(j)   = nanmean(therm_temp{j}(thrm_in_bucket));
    end
    logRWTavg(i,:)=bktRWTavg;
    logRWTstd(i,:)=bktRWTavg;
    logTEMP(i,:)  =bktTEMP;
    % 4) least squares fit regression (slope/intercept) and save in a fluorometer calibration file.
    % there appears to be a measurement threshold of 1ppb for some
    % instruments... lets call the "dark counts" the average of these
    % measurements
    valid_dark = ~isnan(bktRWTavg') & calib_ppb>=0.5 & calib_ppb<=2;
    valid = ~isnan(bktRWTavg') & calib_ppb>=0.5;
    darkcounts = nanmedian(bktRWTavg(valid_dark));%nansum(bktRWTavg(valid_dark)'.*calib_ppb(valid_dark))./nansum(calib_ppb(valid_dark));% the average values are around 1ppb
    X     = bktRWTavg(valid)'-darkcounts;
    L2FIT = X\(calib_ppb(valid)-1);% stoped here... making this work is tricky!
    darkcounts = darkcounts-1/L2FIT;
    X     = bktRWTavg(:)-darkcounts;
    % weighted least squares regression... want good quality at low-mid concentration
% $$$     valid = ~isnan(bktRWTavg') & calib_ppb>=0.5;
% $$$     X     = [ones(sum(valid),1), bktRWTavg(valid)'];
% $$$     Y     = calib_ppb(valid);
% $$$     W     = diag(1./(1+(Y)/100));
% $$$     WL2FIT= ((W*X)'*(W*X))\((W*X)'*(W*Y));
% $$$     darkcounts = -WL2FIT(1)/WL2FIT(2);
% $$$     X     = bktRWTavg(:)-darkcounts;
% $$$     L2FIT = WL2FIT(2);
    % 5) plot results, log slope/offset and temp_avg of buckets
    xm = 2.5;
    ym = 2.5;
    pw = 10;
    ph = 10*(range([min(0,nanmin(bktRWTavg)) max(500,nanmax(bktRWTavg))])/range(calib_ppb));
    ppos = [xm ym pw ph];
    ps   = [2*xm+pw, 2*ym+ph];
    fig  = figure('units','centimeters');
    pos  = get(fig,'position');
    pos(3:4)=ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
    ax   = axes('units','centimeters','position',ppos);
    plot(calib_ppb,bktRWTavg(:),'*r',calib_ppb,X*L2FIT,'or',calib_ppb,calib_ppb,'--k')
    err = 100*abs((X*L2FIT)-calib_ppb)./calib_ppb;
    ERR = mean(err(5:end));
    title(sprintf(' SN: %d    TYPE: %s ',SN,FL_type),'interpreter','latex')
    text(10,450,{sprintf('Avg. \\%%-Error: %2.2f',ERR);['offset: ',num2str(darkcounts)];['slope: ',num2str(L2FIT)]},'interpreter','latex')
    xlabel('bucket concentration [ppb]','interpreter','latex')
    ylabel('sensor [ppb]','interpreter','latex')
    h   = legend('Raw','LS-fit');
    set(h,'location','southeast','interpreter','latex')
    set(ax,'xlim',[calib_ppb(1) calib_ppb(end)],'ylim',[min(0,nanmin(bktRWTavg)) max(500,nanmax(bktRWTavg))],'ticklabelinterpreter','latex','tickdir','out')
    figname = sprintf('%sSN%d_%s_calibration_curve.png',figRoot,SN,FL_type);
    exportgraphics(fig,figname)
    logOFFSET(i) = darkcounts;
    logSLOPE(i)  = L2FIT;
    logERR(i)    = ERR;
end
%
fout = [fluoRoot, 'fluorometer_calibration_coefficients.csv'];
tmp  = cell(Ninst+1,4);
tmp(1,:) = {'Serial Number', 'Offset', 'Slope', 'Temperature'};
tmp(2:Ninst+1,:) = mat2cell([double(bucket_SN(:)), logOFFSET(:), logSLOPE(:), nanmean(logTEMP,2)],ones(Ninst,1),ones(1,4));
writecell(tmp,fout)
%
%
% now use the entire calibration set to estimate the actual bucket concentrations:
X = logRWTavg - logOFFSET';
D = X.*logSLOPE';
%
% now get mean, and remove points outside 2-std
Davg = nanmean(D,1);
Dstd = nanstd(D,[],1);
valid = D>=ones(Ninst,1)*(Davg-Dstd) & D<=ones(Ninst,1)*(Davg+Dstd);
valid(find(sum(valid,2)<10),:) = 0;
%
%
calib_ppb_avg = round(1000*nansum(D.*valid,1)./sum(valid,1))'/1000;
flag = double(valid); flag(flag==0)=nan;
% 7) re-plot results, log slope/offset and temp_avg of buckets
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 10*(range([min(0,nanmin(bktRWTavg)) max(500,nanmax(bktRWTavg))])/range(calib_ppb));
ppos = [xm ym pw ph];
ps   = [2*xm+pw, 2*ym+ph];
fig  = figure('units','centimeters');
pos  = get(fig,'position');
pos(3:4)=ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
ax   = axes('units','centimeters','position',ppos);
plot(flag.*calib_ppb',(flag.*D),'.r',flag.*calib_ppb_avg',(flag.*D),'.k',calib_ppb,calib_ppb,'--k','markersize',8)
err = abs((flag.*D)'-calib_ppb_avg);
ERR = nanmean(nanmean(err(4:end),2));
title(' All Instruments ','interpreter','latex')
text(10,450,{sprintf('Avg. Error: %2.2f ppb',ERR);},'interpreter','latex')
xlabel('bucket concentration [ppb]','interpreter','latex')
ylabel('sensor [ppb]','interpreter','latex')
h   = legend('Target Bucket [ppb]','Estimated Bucket [ppb]');
set(h,'location','southeast','interpreter','latex')
set(ax,'xlim',[calib_ppb(1) calib_ppb(end)],'ylim',[min(0,nanmin(bktRWTavg)) max(500,nanmax(bktRWTavg))],'ticklabelinterpreter','latex','tickdir','out')
figname = sprintf('%sAll_Instrument_calibration_curve.png',figRoot);
exportgraphics(fig,figname)
%
%
%
% now re-do the fits
valid_cal  = calib_ppb_avg>=0 & (calib_ppb_avg-calib_ppb)./(0.1+calib_ppb)<1;
N          = sum(valid_cal);
redoOFFSET = nan(Ninst,1);
redoSLOPE  = nan(Ninst,1);
for i = 1:Ninst
    X     = [ones(N,1), logRWTavg(i,find(valid_cal))'];
    L2FIT = X\calib_ppb_avg(valid_cal);
    redoOFFSET(i) = -L2FIT(1)/L2FIT(2);
    redoSLOPE(i)  = L2FIT(2);
end
%
%
%
% now reconstruct the entire calibration 
X = logRWTavg - redoOFFSET;
D2 = X.*redoSLOPE;
%
% now get mean, and remove points outside 2-std
D2avg = nanmean(D2,1);
D2std = nanstd(D2,[],1);
valid = D2>=ones(Ninst,1)*(D2avg-D2std) & D2<=ones(Ninst,1)*(D2avg+D2std);
valid(find(sum(valid,2)<10),:) = 0;
flag = double(valid); flag(flag==0)=nan;
% 7) re-plot results, log slope/offset and temp_avg of buckets
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 10*(range([min(0,nanmin(bktRWTavg)) max(500,nanmax(bktRWTavg))])/range(calib_ppb));
ppos = [xm ym pw ph];
ps   = [2*xm+pw, 2*ym+ph];
fig  = figure('units','centimeters');
pos  = get(fig,'position');
pos(3:4)=ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
ax   = axes('units','centimeters','position',ppos);
plot(flag.*calib_ppb_avg',(flag.*D),'.r',flag.*calib_ppb_avg',(flag.*D2),'.k',calib_ppb,calib_ppb,'--k','markersize',8)
err2 = abs((flag.*D2)'-calib_ppb_avg);
ERR2 = nanmean(nanmean(err2(4:end),2));
title(' All Instruments ','interpreter','latex')
text(10,450,{sprintf('Avg. Error: %2.2f ppb',ERR2);},'interpreter','latex')
xlabel('bucket concentration [ppb]','interpreter','latex')
ylabel('sensor [ppb]','interpreter','latex')
h   = legend('Preliminary Fit','Secondary Fit');
set(h,'location','southeast','interpreter','latex')
set(ax,'xlim',[calib_ppb(1) calib_ppb(end)],'ylim',[min(0,nanmin(bktRWTavg)) max(500,nanmax(bktRWTavg))],'ticklabelinterpreter','latex','tickdir','out')
figname = sprintf('%sAll_Instrument_calibration_curve_v2.png',figRoot);
exportgraphics(fig,figname)
%
%
%
A = struct('calibration_data_directory',fluoRoot,'calibration_data_filenames',char(cal_file_names),'calibration_start_date_UTC',cal_start_day_str,'calibration_serial_number_and_times_file',calibration_times_file,'calibration_bucket_target_concentration_ppb',calib_ppb,'calibration_bucket_serial_number',bucket_SN,'calibration_bucket_times',bucket_times,'calibration_bucket_estimated_concentration_ppb',calib_ppb_avg,'fluorometer_serial_number_and_type_filename',fluorometer_info_file,'fluorometer_measured_average_concentrations',logRWTavg,'fluorometer_measured_standard_deviation_concentration',logRWTstd,'thermistor_info_filename',thermistor_info_file,'thermistor_serial_number',therm_SN,'thermistor_bucket_target_concentration_ppb',therm_ppb,'thermistor_measured_average_temperature',logTEMP,'estimated_offset',logOFFSET,'estimated_slope',logSLOPE)
%
cd(fluoRoot)
ncfile = ['fluorometer_calibration.nc'];
struct2nc(A,ncfile,'NETCDF4')