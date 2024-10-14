clear all
close all
%
% 0) define directories/file names & start/end times
timeIn  = datenum('07/03/2024 00:00:00');
timeOut = datenum('07/25/2024 00:00:00');
time    = timeIn:1/86400:timeOut;
mooring_info_file = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/info/mooring_designs_and_seial_numbers/mooring_design_and_serial_numbers.csv';
mooring_loc_file  = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/info/mooring_designs_and_seial_numbers/mooring_lat_lon.csv';
dye_cal_file      = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/dye_calibrations/fluorometer_calibration_coefficients.csv';
data_dir          = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/';
archive_dir       = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/2024_PROCESSED_DATA/';
eco_dir = [data_dir,'ECOtriplet',filesep]; 
sbe_dir = [data_dir,'SBE56s',filesep];
rbr_dir = [data_dir,'RBR',filesep];
pme_root= [data_dir,'PME_'];
ysi_dir = [data_dir,'YSI',filesep];
%
% 1) get mooring/instrument/depth info
fid = fopen(mooring_info_file);
mooring_vars = fgetl(fid);
mooring_vars = split(mooring_vars,',');
mooring_info = textscan(fid,'%s %f %s %f','delimiter',',');
mooring_name = deblank(mooring_info{1});
instrument_mab = mooring_info{2};
instrument_type= mooring_info{3};
instrument_SN  = mooring_info{4};
%
%
% 2a) load the mooring deployment location/time/depth
fid = fopen(mooring_loc_file);
mooring_loc_vars = fgetl(fid);
mooring_loc_vars = split(mooring_loc_vars,',');
mooring_loc_data = textscan(fid,'%s %f %f %s %f','delimiter',',');
%
% 2b) split into each mooring,
%    create a mooring structure
moorings = unique(mooring_name);
% loop over moorings
fig0=figure;
fig1=figure;
for ii = 8:length(moorings)
    out = struct([]);
    ind_loc = ismember(mooring_loc_data{1},moorings{ii});
    ind_info= ismember(mooring_name,moorings{ii});
    %
    % log lat/lon/time/depth
    out(1).Name     = moorings{ii};
    out(1).Latitude = mooring_loc_data{2}(ind_loc');
    out(1).Longitude = mooring_loc_data{3}(ind_loc');
    out(1).Time_deploy = char(mooring_loc_data{4}(ind_loc',:));
    out(1).Depth_deploy = mooring_loc_data{5}(ind_loc');
%
% 3) load and interpolate data
% 3c) switch based on instrument type,
%     create instrument specific structure,
%     save to /L0/instrument_type/
    temp = [];
    temp_depth = [];
    dye = [];
    dye_depth = [];
    pres = [];
    pres_depth = [];
    salt = [];
    salt_depth = [];
    inds = find(ind_info);
    fdom = [];
    fdom_depth = [];
    turb = [];
    turb_depth = [];
    for jj = 1:sum(ind_info)
        type = char(instrument_type(inds(jj)));
        SN   = instrument_SN(inds(jj));
        mab  = instrument_mab(inds(jj));
        switch type
          case 'SBE'
            sbeFileStr = sprintf(['%s',filesep,'SBE56s',filesep,'SN_%d.csv'],data_dir,SN);
            sbeFile = dir(sbeFileStr);
            fid      = fopen([sbeFile.folder,filesep,sbeFile.name]);
            if fid==-1
                disp(['missing data file: ', sbeFileStr])
                continue
            end
            hdrflag  = 0;
            while ~hdrflag
                line = fgetl(fid);
                strings = split(line,{',','"'});
                if ismember('Date',strings);
                    hdrflag = 1;
                end
            end
            % read data
            sbe_data = textscan(fid,'%q %q %q','delimiter',',');
            Nsbe = size(sbe_data{1},1);
            sbe_time = datenum([char(sbe_data{1}),repmat('-',Nsbe,1),char(sbe_data{2})],'yyyy-mm-dd-HH:MM:SS');
            sbe_temp = cellfun(@str2num,sbe_data{3});
            % put info into structure array
            sbe = struct('Time',sbe_time,'Temperature',sbe_temp,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            sbe_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'SBE_56',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(sbe_file,'-struct','sbe')
            ncfile = [sbe_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(sbe,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            temp_interp = interp1(sbe_time,sbe_temp,time);
            temp = cat(1,temp,temp_interp);
            temp_depth = cat(1,temp_depth,mab);
            %
            figure(fig0),clf
            plot(datetime(sbe_time,'convertfrom','datenum'),sbe_temp,'-r')
            ylabel('$T$ [$^\circ$ C]','interpreter','latex')
            title(sprintf('SBE-SN%04d',SN),'interpreter','latex')
            set(gca,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            exportgraphics(fig0,[sbe_file,'.pdf'])
          case {'RBR','SOLO','DUET'}
            rbrFileStr = sprintf([data_dir,filesep,'RBR',filesep,'SN_%04d*.rsk'],SN);
            rbrFile = dir(rbrFileStr);
            fin     = [rbrFile.folder,filesep,rbrFile.name];
            % open the file, then read data
            try rsk = RSKopen(fin);
            catch
                disp(['missing data file: ', rbrFileStr])
                continue
            end
            rsk = RSKreaddata(rsk);
            rbr_time = rsk.data.tstamp;
            rbr_data = rsk.data.values;
            rbr_temp = [];
            rbr_pres = [];
            switch rsk.instruments.model
              case 'RBRsolo'
                if ~strcmp(rsk.channels.longName,'Temperature')
                    disp('RBRsolo is not a thermistor')
                    keyboard
                end
                rbr_temp = rbr_data;
              case 'RBRduet'
                disp('duet')
                return
                rbr_temp = rbr_data(:,1);
                rbr_pres = rbr_data(:,2);                
            end
            % put info into structure array
            rbr = struct('Time',rbr_time,'Temperature',rbr_temp,'Pressure',rbr_pres,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            rbr_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'RBR_Temperature',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(rbr_file,'-struct','rbr')
            ncfile = [rbr_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(rbr,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            temp_interp = interp1(rbr_time,rbr_temp,time);
            temp = cat(1,temp,temp_interp);
            temp_depth = cat(1,temp_depth,mab);
            presFlag   = 0;
            if ~isempty(rbr_pres)
                % interpolate to time vector and prep to save in mooring struct
                pres_interp = interp1(rbr_time,rbr_pres,time);
                pres = cat(1,pres,pres_interp);
                pres_depth = cat(1,pres_depth,mab);
                presFlag   = 1;
            end
            figure(fig0),clf
            ax1 = subplot(1+presFlag,1,1);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_temp,'-r')
            ylabel(ax1,'$T$ [$^\circ$ C]','interpreter','latex')
            title(sprintf('RBR-SN%04d',SN),'interpreter','latex')
            set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            if presFlag
                ax2 = subplot(1+presFlag,1,1);
                plot(datetime(rbr_time,'convertfrom','datenum'),rbr_temp,'-r')
                ylabel(ax2,'$p$ [dbar]','interpreter','latex')
                set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            end
            exportgraphics(fig0,[rbr_file,'.pdf'])
          case {'ECO','ECO-BIO','ECO-T'}
            format = '%s %s %f %f %f %f %f %f %f';% 
            if all(strcmp(type,'ECO-T'))
                format = cat(2,format,' %f');
            end
            fin = sprintf([eco_dir,filesep,'%d.raw'],SN);
            fid = fopen(fin,'r');
            str = fgetl(fid);
            N   = str2num(str(1:strfind(str,'records')-1));
            dat = textscan(fid,format,'delimiter','\t');
            mmddyy = dat{1}(1:N);
            HHMMSS = dat{2}(1:N);
            fl_time = datenum( strcat(mmddyy, repmat('-',N,1), HHMMSS), 'mm/dd/yy-HH:MM:SS');
            NTU  = dat{4}(1:N);
            CHL  = dat{6}(1:N);
            RWT  = dat{8}(1:N);
            flag = (RWT>4096 | CHL>4096 | NTU>4096);
            NTU(flag)=nan; CHL(flag)=nan; RWT(flag)=nan;
            % despike raw data
            % RWT = despike_fluoro(RWT,1);
            % convert to ppb (currently neglecting temperature effects)
            RWT = convert_fluorometer_counts_to_ppb(SN,RWT,[],dye_cal_file);
            %
            % put info into structure array
            eco = struct('Time',fl_time,'Dye',RWT,'Chlorophyl',CHL,'Turbidity',NTU,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            eco_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'EcoTriplet',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(eco_file,'-struct','eco')
            ncfile = [eco_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(eco,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            dye_interp = interp1(fl_time,RWT,time);
            dye = cat(1,dye,dye_interp);
            dye_depth = cat(1,dye_depth,mab);
            %
            figure(fig0),clf
            ax1 = subplot(3,1,1);
            plot(datetime(fl_time,'convertfrom','datenum'),RWT,'-m')
            ylabel(ax1,'$D$ [ppb]','interpreter','latex')
            title(sprintf('ECO-SN%04d',SN),'interpreter','latex')
            set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax2 = subplot(3,1,2);
            plot(datetime(fl_time,'convertfrom','datenum'),CHL,'-r')
            ylabel(ax2,'Chl-a [counts]','interpreter','latex')
            set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax3 = subplot(3,1,3);
            plot(datetime(fl_time,'convertfrom','datenum'),NTU,'-r')
            ylabel(ax3,'Turbidity [counts]','interpreter','latex')
            set(ax3,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            exportgraphics(fig0,[eco_file,'.pdf'])
          case 'PME'
            format   = '%f %f %f %f %d';
            pme_dir  = sprintf([pme_root,'%d',filesep],SN);
            pme_files= dir([pme_dir,'*.txt']);
            nf       = length(pme_files);
            pme_time = [];
            pme_temp = [];
            pme_dye  = [];
            for kk = 1:nf
                fin = [pme_dir, pme_files(kk).name];
                fid = fopen(fin,'r');
                dat = textscan(fid,format,'delimiter',',','headerlines',3);
                SS  = dat{1};
                tmp = datenum('Jan 01 1970')+SS/86400;
                pme_time = cat(1,pme_time,tmp);
                tmp1 = dat{3};% temperature
                tmp2 = dat{4};% dye concentration
                RWT  = convert_fluorometer_counts_to_ppb(SN,tmp2,tmp1,dye_cal_file);
                pme_dye  = cat(1,pme_dye,RWT);
                pme_temp = cat(1,pme_temp,tmp1);
            end
            if SN==339378
                % based on max in lagged correlation between PME-339378 and SBE56-13592
                pme_time = pme_time + 2.7040;
                % based on the difference between calibration bucket times and dye signal
                % pme_time = pme_time + (datenum('30-Jun-2024 22:45:00') - datenum('28-Jun-2024 05:45:00'));
            end
            % put info into structure array
            pme = struct('Time',pme_time,'Dye',pme_dye,'Temperature',pme_temp,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            pme_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'PME',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(pme_file,'-struct','pme')
            ncfile = [pme_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(pme,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            dye_interp = interp1(pme_time,pme_dye,time);
            dye = cat(1,dye,dye_interp);
            dye_depth = cat(1,dye_depth,mab);
            %
            temp_interp = interp1(pme_time,pme_temp,time);
            temp = cat(1,temp,temp_interp);
            temp_depth = cat(1,temp_depth,mab);
            %
            figure(fig0),clf
            ax1 = subplot(2,1,1);
            plot(datetime(pme_time,'convertfrom','datenum'),pme_dye,'-m')
            ylabel(ax1,'$D$ [ppb]','interpreter','latex')
            title(sprintf('PME-SN%04d',SN),'interpreter','latex')
            set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax2 = subplot(2,1,2);
            plot(datetime(pme_time,'convertfrom','datenum'),pme_temp,'-r')
            ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
            set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            exportgraphics(fig0,[pme_file,'.pdf'])
          case {'RBR-C-F','RBR-V-F'}
            rbrFileStr = sprintf([data_dir,filesep,'RBR',filesep,'SN_%04d*.rsk'],SN);
            rbrFile = dir(rbrFileStr);
            fin     = [rbrFile.folder,filesep,rbrFile.name];
            % open the file, then read data
            try rsk = RSKopen(fin);
            catch
                disp(['missing data file: ', rbrFileStr])
                continue
            end
            rsk = RSKreaddata(rsk);
            rbr_time = rsk.data.tstamp;
            rbr_data = rsk.data.values;
            RWT  = rbr_data(:,1);
            rbr_dye = convert_fluorometer_counts_to_ppb(SN,RWT,[],dye_cal_file);
            %
            if strcmp(type,'RBR-C-F')
                rbr_fdom = rbr_data(:,2);
                rbr_turb = rbr_data(:,3);
            else
                rbr_fdom = [];
                rbr_turb = [];
            end
            %
            % put info into structure array
            rbr = struct('Time',rbr_time,'Dye',rbr_dye,'CDOM',rbr_fdom,'Turbidity',rbr_turb,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            rbr_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'RBR_Fluorometer',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(rbr_file,'-struct','rbr')
            ncfile = [rbr_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(rbr,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            dye_interp = interp1(rbr_time,rbr_dye,time);
            dye = cat(1,dye,dye_interp);
            dye_depth = cat(1,dye_depth,mab);
            %
            turbFlag = 0;
            if ~isempty(rbr_fdom)
                % interpolate to time vector and prep to save in mooring struct
                fdom_interp = interp1(rbr_time,rbr_fdom,time);
                fdom = cat(1,fdom,fdom_interp);
                fdom_depth = cat(1,fdom_depth,mab);
                %
                turb_interp = interp1(rbr_time,rbr_turb,time);
                turb = cat(1,turb,turb_interp);
                turb_depth = cat(1,turb_depth,mab);
                turbFlag = 1;
            end
            figure(fig0),clf
            ax1 = subplot(1+2*turbFlag,1,1);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_dye,'-m')
            ylabel(ax1,'$D$ [ppb]','interpreter','latex')
            title(sprintf('RBR-SN%04d',SN),'interpreter','latex')
            set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            if turbFlag
            ax2 = subplot(1+2*turbFlag,1,2);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_fdom,'-r')
            ylabel(ax2,'fDOM [?]','interpreter','latex')
            set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax3 = subplot(3,1,3);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_turb,'-r')
            ylabel(ax3,'Turbidity [FTU]','interpreter','latex')
            set(ax3,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            end
            exportgraphics(fig0,[rbr_file,'.pdf'])
          case 'RBR-CTD'
            rbrFileStr = sprintf([data_dir,filesep,'RBR',filesep,'21%04d*.rsk'],SN);
            rbrFile = dir(rbrFileStr);
            fin     = [rbrFile.folder,filesep,rbrFile.name];
            % open the file, then read data
            try rsk = RSKopen(fin);
            catch
                disp(['missing data file: ', rbrFileStr])
                continue
            end
            rsk = RSKreaddata(rsk);
            rbr_time = rsk.data.tstamp;
            rbr_data = rsk.data.values;
            rbr_cond = rbr_data(:,1);
            rbr_temp = rbr_data(:,2);
            rbr_pres = rbr_data(:,3);
            % need to offset pressure,
            rbr_pres_offset = mean(rbr_pres(rbr_pres<10.5 & rbr_pres>9.5));
            rbr_salt = gsw_SP_from_C(rbr_cond,rbr_temp,rbr_pres-rbr_pres_offset);
            % put info into structure array
            rbr = struct('Time',rbr_time,'Temperature',rbr_temp,'Pressure',rbr_pres-rbr_pres_offset,'PressureOffset',rbr_pres_offset,'Salinity',rbr_salt,'Conductivity',rbr_cond,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            rbr_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'RBR_CTD',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(rbr_file,'-struct','rbr')
            ncfile = [rbr_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(rbr,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            temp_interp = interp1(rbr_time,rbr_temp,time);
            temp = cat(1,temp,temp_interp);
            temp_depth = cat(1,temp_depth,mab);
            % interpolate to time vector and prep to save in mooring struct
            pres_interp = interp1(rbr_time,rbr_pres-rbr_pres_offset,time);
            pres = cat(1,pres,pres_interp);
            pres_depth = cat(1,pres_depth,mab);
            % interpolate to time vector and prep to save in mooring struct
            salt_interp = interp1(rbr_time,rbr_salt,time);
            salt = cat(1,salt,salt_interp);
            salt_depth = cat(1,salt_depth,mab);
            %
            figure(fig0),clf
            ax1 = subplot(3,1,1);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_temp,'-r')
            ylabel(ax1,'$T$ [$^\circ$C]','interpreter','latex')
            title(sprintf('RBR-SN%04d',SN),'interpreter','latex')
            set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax2 = subplot(3,1,2);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_pres,'-k')
            ylabel(ax2,'Pressure [dbar]','interpreter','latex')
            set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax3 = subplot(3,1,3);
            plot(datetime(rbr_time,'convertfrom','datenum'),rbr_salt,'-b')
            ylabel(ax3,'Salinity [PSU]','interpreter','latex')
            set(ax3,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            exportgraphics(fig0,[rbr_file,'.pdf'])
          case 'YSI'
            format = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
            fin = sprintf([ysi_dir,filesep,'SN_%04d.csv'],SN);
            fid = fopen(fin,'r');
            if fid==-1
                disp(['missing data file: ', fin])
                continue
            end
            ysi_data = textscan(fid,format,'delimiter',',','headerlines',11);
            Nysi     = size(ysi_data{1},1);
            ysi_time = flipud(datenum( strcat(ysi_data{2},repmat('-',Nysi,1),ysi_data{1}) ));
            ysi_cond = flipud(ysi_data{6});
            ysi_salt = flipud(ysi_data{9});
            ysi_pres = flipud(ysi_data{13});            
            ysi_RWT  = flipud(ysi_data{15});
            ysi_temp = flipud(ysi_data{18});
            ysi_dye = convert_fluorometer_counts_to_ppb(SN,ysi_RWT,ysi_temp,dye_cal_file);
            %
            ysi_pres_offset = mean(ysi_pres<0.5 & ysi_pres>0.1);
            %
            % put info into structure array
            ysi = struct('Time',ysi_time,'Dye',ysi_dye,'Temperature',ysi_temp,'Pressure',ysi_pres-ysi_pres_offset,'PressureOffset',ysi_pres_offset,'Salinity',ysi_salt,'Conductivity',ysi_cond,'Latitude',out.Latitude,'Longitude',out.Longitude,'Time_deploy',out.Time_deploy,'Depth_deploy',out.Depth_deploy,'mab',mab);
            ysi_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L0',filesep,'YSI',filesep,num2str(SN,'%04d')];
            % save raw data to .nc and .mat in mooring dir
            save(ysi_file,'-struct','ysi')
            ncfile = [ysi_file,'.nc'];
            if exist(ncfile)
                eval(['!rm ',ncfile])
            end
            struct2nc(ysi,ncfile,'NETCDF4')
            %
            % interpolate to time vector and prep to save in mooring struct
            dye_interp = interp1(ysi_time,ysi_dye,time);
            dye = cat(1,dye,dye_interp);
            dye_depth = cat(1,dye_depth,mab);
            % interpolate to time vector and prep to save in mooring struct
            temp_interp = interp1(ysi_time,ysi_temp,time);
            temp = cat(1,temp,temp_interp);
            temp_depth = cat(1,temp_depth,mab);
            % interpolate to time vector and prep to save in mooring struct
            pres_interp = interp1(ysi_time,ysi_pres-ysi_pres_offset,time);
            pres = cat(1,pres,pres_interp);
            pres_depth = cat(1,pres_depth,mab);
            % interpolate to time vector and prep to save in mooring struct
            salt_interp = interp1(ysi_time,ysi_salt,time);
            salt = cat(1,salt,salt_interp);
            salt_depth = cat(1,salt_depth,mab);
            %
            figure(fig0),clf
            ax1 = subplot(4,1,1);
            plot(datetime(ysi_time,'convertfrom','datenum'),ysi_dye,'-r')
            ylabel(ax1,'$D$ [ppb]','interpreter','latex')
            title(sprintf('YSI-SN%04d',SN),'interpreter','latex')
            set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax2 = subplot(4,1,2);
            plot(datetime(ysi_time,'convertfrom','datenum'),ysi_temp,'-r')
            ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
            set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax3 = subplot(4,1,3);
            plot(datetime(ysi_time,'convertfrom','datenum'),ysi_pres,'-k')
            ylabel(ax3,'Pressure [dbar]','interpreter','latex')
            set(ax3,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            %
            ax4 = subplot(4,1,4);
            plot(datetime(ysi_time,'convertfrom','datenum'),ysi_salt,'-b')
            ylabel(ax4,'Salinity [PSU]','interpreter','latex')
            set(ax4,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
            exportgraphics(fig0,[ysi_file,'.pdf'])
        end
    end
    % 3b) interpolate to mooring time vector,
    %     log instrument depth, and other notes
    out(1).Time        = time;
    out(1).Temperature = temp;
    out(1).Temperature_mab = temp_depth;
    out(1).Dye         = dye;
    out(1).Dye_mab = dye_depth;
    out(1).CDOM         = fdom;
    out(1).CDOM_mab = fdom_depth;
    out(1).Turbidity         = turb;
    out(1).Turbidity_mab = turb_depth;
    out(1).Pressure   = pres;
    out(1).Pressure_mab = pres_depth;
    out(1).Salinity   = salt;
    out(1).Salinity_mab = salt_depth;
    out(1).Instrument_SN= instrument_SN(inds);
    %
    %
    % 4) archive to /L1/instrument_type/
    out_file = [archive_dir,filesep,deblank(moorings{ii}),filesep,'L1',filesep,'mooring_',deblank(moorings{ii})];
    % save raw data to .nc and .mat in mooring dir
    save(out_file,'-struct','out')
    ncfile = [out_file,'.nc'];
    if exist(ncfile)
        eval(['!rm ',ncfile])
    end
    struct2nc(out,ncfile,'NETCDF4')    
    %
    presFlag = size(pres,2)==length(time);
    saltFlag = size(salt,2)==length(time);
    %
    figure(fig1),clf
    ax1 = subplot(2+presFlag+saltFlag,1,1);
    plot(datetime(time,'convertfrom','datenum'),(dye))
    ylabel(ax1,'$D$ [ppb]','interpreter','latex')
    title(sprintf('Mooring-%s',moorings{ii}),'interpreter','latex')
    set(ax1,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
    %
    ax2 = subplot(2+presFlag+saltFlag,1,2);
    plot(datetime(time,'convertfrom','datenum'),temp)
    ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
    set(ax2,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
    %
    if presFlag
        ax3 = subplot(2+presFlag+saltFlag,1,3);
        ylabel(ax3,'Pressure [dbar]','interpreter','latex')
        plot(datetime(time,'convertfrom','datenum'),pres)
        set(ax3,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')        
    end
    %
    if saltFlag
        ax4 = subplot(2+presFlag+saltFlag,1,3+presFlag);
        plot(datetime(time,'convertfrom','datenum'),salt)
        ylabel(ax4,'Salinity [PSU]','interpreter','latex')
        set(ax4,'tickdir','out','xlim',datetime([min(time) max(time)],'convertfrom','datenum'),'ticklabelinterpreter','latex')
    end        
    exportgraphics(fig1,[out_file,'.pdf'])
end
%
%
%
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//RBR/SN_7222*.rsk
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//RBR/SN_1176*.rsk
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//RBR/SN_1150*.rsk
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//RBR/SN_1932*.rsk
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//RBR/SN_1939*.rsk
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/YSI//SN_5124.csv
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/YSI//SN_5120.csv
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/YSI//SN_5118.csv
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/YSI//SN_5127.csv
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//SBE56s/SN_13532.csv
% $$$ missing data file: /Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw//SBE56s/SN_13646.csv
% $$$ >> 