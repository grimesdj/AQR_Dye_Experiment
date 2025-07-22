clear all
close all
%
% 1) need a list of serial numbers for instruments (see: info/notes/KELP-Dye release notes)
% 2) need the directory list for raw .rbr files
% 3) need a table of depths and deployment times (i
root_dir         = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/';
data_dir         = [root_dir,'/data/Release1/raw/'];
archive_dir      = [root_dir,'/data/Release%d/L0/'];
instruments_file = [root_dir,'/info/dye_mixing_cals_and_releases/dye_release_instruments.csv'];
%
release_times    = readtable([root_dir,'info/dye_mixing_cals_and_releases/dye_release_times.csv']);
%               Type    , SN    , Elevation (mab)  , notes
format   = '%s %d %f %s';
fid      = fopen(instruments_file);
info     = textscan(fid,format,'delimiter',',','headerlines',1);
TYPEs    = deblank(info{1});
SNs      = info{2};
MABs     = info{3};
notes    = info{4};
% 4) need a directory to save the L0 files, but this will be release specific... so hold off
% L0_dir = [root_dir,'data/Release%d/L0/'];
% 5) loop and load data into structure, estimate density (use CTD and assume
% salt is constant)
% 6) estimate N^2 over each vertical interval
% 7) calculate the mean and standard deviation of temperature during
% each release for appropriate instrument?
% 8) make some plots
idx = ismember(TYPEs,{'DUET','RBR-CTD'});
N   = sum(idx);
%
R1  = struct('Time',[],'Temperature',[],'Pressure',[],'PressureOffset',[],'Salinity',[],'Conductivity',[],'Depth_deploy',[],'mab',[]);
R23 = struct('Time',[],'Temperature',[],'Pressure',[],'PressureOffset',[],'Salinity',[],'Conductivity',[],'Depth_deploy',[],'mab',[]);
R2  = struct('Time',[],'Temperature',[],'Pressure',[],'PressureOffset',[],'Salinity',[],'Conductivity',[],'Depth_deploy',[],'mab',[]);
R3  = struct('Time',[],'Temperature',[],'Pressure',[],'PressureOffset',[],'Salinity',[],'Conductivity',[],'Depth_deploy',[],'mab',[]);
for ii = 1:N% i put the adcps at the end, so this shouldn't barf
    % get current SN
    SN   = SNs(ii);
    type = TYPEs{ii};
    mab  = MABs(ii);
    %
    % load/process... below code is from "process_moored_data.m
    switch type
      case {'RBR','SOLO','DUET'}
        rbrFileStr = sprintf([data_dir,filesep,'SN_%04d*.rsk'],SN);
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
        switch rsk.instruments.model(1:7)
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
        % pressure offset!
        rbr_pres_offset = mean(rbr_pres(rbr_pres<10.5 & rbr_pres>9.5),'omitnan');
        % put info into structure array
        rbr = struct('Time',rbr_time,'Temperature',rbr_temp,'Pressure',rbr_pres-rbr_pres_offset,'PressureOffset',rbr_pres_offset,'Depth_deploy',mean(rbr_pres(rbr_pres>10.5)-rbr_pres_offset),'mab',mab);
        rbr_file = [data_dir,filesep,'..',filesep,'L0',filesep,num2str(SN,'21%04d')];
        % save raw data to .nc and .mat in mooring dir
        save(rbr_file,'-struct','rbr')
      case 'RBR-CTD'
        rbrFileStr = sprintf([data_dir,filesep,'%06d*.rsk'],SN);
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
        rbr = struct('Time',rbr_time,'Temperature',rbr_temp,'Pressure',rbr_pres-rbr_pres_offset,'PressureOffset',rbr_pres_offset,'Salinity',rbr_salt,'Conductivity',rbr_cond,'Depth_deploy',mean(rbr_pres(rbr_pres>10.5)-rbr_pres_offset),'mab',mab);
        rbr_file = [data_dir,filesep,'..',filesep,'L0',filesep,num2str(SN,'%04d')];
        % save raw data to .nc and .mat in mooring dir
        save(rbr_file,'-struct','rbr')
    end
    %
    % now loop over rbr_fields to extract each release's data...
    % first get Release1,then combined Release2-Release3, then individual Release2 & Release3,
    Nt  = size(rbr.Time,1);
    iR1 = find(rbr.Time>=datenum(release_times.StartTime_UTC_(1)) & rbr.Time<=datenum(release_times.EndTime_UTC_(1)));
    iR23 = find(rbr.Time>=datenum(release_times.StartTime_UTC_(2)) & rbr.Time<=datenum(release_times.EndTime_UTC_(3)));
    iR2 = find(rbr.Time>=datenum(release_times.StartTime_UTC_(2)) & rbr.Time<=datenum(release_times.EndTime_UTC_(2)));
    iR3 = find(rbr.Time>=datenum(release_times.StartTime_UTC_(3)) & rbr.Time<=datenum(release_times.EndTime_UTC_(3)));
    % now loop over field names
    fields = fieldnames(rbr);
    fields = fields(~ismember(fields,'Time'));
    if isempty(R1.Time)
        R1(1).Time = rbr.Time(iR1);
        R2(1).Time = rbr.Time(iR2);
        R23(1).Time = rbr.Time(iR23);                
        R3(1).Time = rbr.Time(iR3);
    end
    for jj = 1:length(fields)
        field = fields{jj};
        if length(rbr.(field))==Nt
            tmp        = interp1(rbr.Time,rbr.(field),R1.Time);
            R1.(field) = cat(2,R1.(field),tmp);
            tmp        = interp1(rbr.Time,rbr.(field),R2.Time);
            R2.(field) = cat(2,R2.(field),tmp);
            tmp        = interp1(rbr.Time,rbr.(field),R23.Time);
            R23.(field) = cat(2,R23.(field),tmp);
            tmp        = interp1(rbr.Time,rbr.(field),R3.Time);
            R3.(field) = cat(2,R3.(field),tmp);
        else
            R1.(field) = cat(2,R1.(field),rbr.(field));
            R2.(field) = cat(2,R2.(field),rbr.(field));            
            R23.(field) = cat(2,R23.(field),rbr.(field));
            R3.(field) = cat(2,R3.(field),rbr.(field));
        end
    end
end
%
%
% Now estimate some stats and make some plots
g = 9.81;
[Nt,NT]    = size(R1.Temperature);
R1.Density = gsw_rho(R1.Salinity*ones(1,NT),R1.Temperature,R1.Pressure);
R1.N2      = -g/mean(R1.Density,'all','omitnan').*( diff(R1.Density,1,2)./diff(R1.mab)  );
%
[Nt,NT]    = size(R2.Temperature);
R2.Density = gsw_rho(R2.Salinity*ones(1,NT),R2.Temperature,R2.Pressure);
R2.N2      = -g/mean(R2.Density,'all','omitnan').*( diff(R2.Density,1,2)./diff(R2.mab)  );
%
[Nt,NT]     = size(R23.Temperature);
R23.Density = gsw_rho(R23.Salinity*ones(1,NT),R23.Temperature,R23.Pressure);
R23.N2      = -g/mean(R23.Density,'all','omitnan').*( diff(R23.Density,1,2)./diff(R23.mab)  );
%
[Nt,NT]    = size(R3.Temperature);
R3.Density = gsw_rho(R3.Salinity*ones(1,NT),R3.Temperature,R3.Pressure);
R3.N2      = -g/mean(R3.Density,'all','omitnan').*( diff(R3.Density,1,2)./diff(R3.mab)  );
%
%
% make some figures
clrs = flipud(cmocean('thermal',length(R1.mab)));
%
% pannels: Pres, Temp, N
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 4;
ag = 0.5;
ppos1 = [xm ym pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ppos3 = [xm ym+2*ph+2*ag pw ph];
cbpos = [xm+0.75*pw ym+2.9*ph+2*ag 0.25*pw 0.1*ph];
ps    = [1.5*xm+pw, 1.5*ym+2*ag+3*ph];
%
fig = figure('units','centimeters','color','w');
pos = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'paperposition',[0 0 ps],'papersize',ps)
colororder(clrs)
%
ax1      = axes('units','centimeters','position',ppos1);
N        = sqrt(R1.N2);
N(imag(N)~=0)=nan;
plot(ax1,datetime(R1.Time,'convertFrom','datenum'),N,'linewidth',1.5)
ylabel(ax1,'$N$ [s$^{-1}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex',...
   'xlim',datetime([release_times.StartTime_UTC_(1); release_times.EndTime_UTC_(1)]))
colororder(clrs)
grid on
%
ax2      = axes('units','centimeters','position',ppos2);
plot(ax2,datetime(R1.Time,'convertFrom','datenum'),R1.Temperature,'linewidth',1.5)
ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
        'xlim',datetime([release_times.StartTime_UTC_(1); release_times.EndTime_UTC_(1)]))
grid on
%
ax2      = axes('units','centimeters','position',ppos3);
plot(ax2,datetime(R1.Time,'convertFrom','datenum'),R1.Pressure(:,end)+R1.mab(end),'k','linewidth',1.5)
ylabel(ax2,'Depth [m]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
      'xlim',datetime([release_times.StartTime_UTC_(1); release_times.EndTime_UTC_(1)]))
grid on
%
cb       = axes('units','centimeters','position',cbpos);
imagesc([1:length(R1.mab)],0,reshape(clrs,1,length(R1.mab),3));
set(cb,'ytick',[],'xtick',[1:length(R1.mab)],'xticklabel',num2str(R1.mab','%1.2f'))
ylabel(cb,'[m]','rotation',0)
%
figname = [root_dir,'/data/Release1/figures/R1_Release_Stratification.pdf'];
exportgraphics(fig,figname)
%
%
fig = figure('units','centimeters','color','w');
pos = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'paperposition',[0 0 ps],'papersize',ps)
colororder(clrs)
%
ax1      = axes('units','centimeters','position',ppos1);
N        = sqrt(R2.N2);
N(imag(N)~=0)=nan;
plot(ax1,datetime(R2.Time,'convertFrom','datenum'),N,'linewidth',1.5)
ylabel(ax1,'$N$ [s$^{-1}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex',...
   'xlim',datetime([release_times.StartTime_UTC_(2); release_times.EndTime_UTC_(2)]))
colororder(clrs)
grid on
%
ax2      = axes('units','centimeters','position',ppos2);
plot(ax2,datetime(R2.Time,'convertFrom','datenum'),R2.Temperature,'linewidth',1.5)
ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
        'xlim',datetime([release_times.StartTime_UTC_(2); release_times.EndTime_UTC_(2)]))
grid on
%
ax2      = axes('units','centimeters','position',ppos3);
plot(ax2,datetime(R2.Time,'convertFrom','datenum'),R2.Pressure(:,end)+R2.mab(end),'k','linewidth',1.5)
ylabel(ax2,'Depth [m]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
      'xlim',datetime([release_times.StartTime_UTC_(2); release_times.EndTime_UTC_(2)]))
grid on
%
cb       = axes('units','centimeters','position',cbpos);
imagesc([1:length(R2.mab)],0,reshape(clrs,1,length(R2.mab),3));
set(cb,'ytick',[],'xtick',[1:length(R2.mab)],'xticklabel',num2str(R2.mab','%1.2f'))
ylabel(cb,'[m]','rotation',0)
figname = [root_dir,'/data/Release2/figures/R2_Release_Stratification.pdf'];
exportgraphics(fig,figname)
%
%
%
fig = figure('units','centimeters','color','w');
pos = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'paperposition',[0 0 ps],'papersize',ps)
colororder(clrs)
%
ax1      = axes('units','centimeters','position',ppos1);
N        = sqrt(R3.N2);
N(imag(N)~=0)=nan;
plot(ax1,datetime(R3.Time,'convertFrom','datenum'),N,'linewidth',1.5)
ylabel(ax1,'$N$ [s$^{-1}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex',...
   'xlim',datetime([release_times.StartTime_UTC_(3); release_times.EndTime_UTC_(3)]))
colororder(clrs)
grid on
%
ax2      = axes('units','centimeters','position',ppos2);
plot(ax2,datetime(R3.Time,'convertFrom','datenum'),R3.Temperature,'linewidth',1.5)
ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
        'xlim',datetime([release_times.StartTime_UTC_(3); release_times.EndTime_UTC_(3)]))
grid on
%
ax2      = axes('units','centimeters','position',ppos3);
plot(ax2,datetime(R3.Time,'convertFrom','datenum'),R3.Pressure(:,end)+R3.mab(end),'k','linewidth',1.5)
ylabel(ax2,'Depth [m]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
      'xlim',datetime([release_times.StartTime_UTC_(3); release_times.EndTime_UTC_(3)]))
grid on
%
cb       = axes('units','centimeters','position',cbpos);
imagesc([1:length(R3.mab)],0,reshape(clrs,1,length(R3.mab),3));
set(cb,'ytick',[],'xtick',[1:length(R3.mab)],'xticklabel',num2str(R3.mab','%1.2f'))
ylabel(cb,'[m]','rotation',0)
figname = [root_dir,'/data/Release3/figures/R3_Release_Stratification.pdf'];
exportgraphics(fig,figname)
%
%
%
%
fig = figure('units','centimeters','color','w');
pos = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'paperposition',[0 0 ps],'papersize',ps)
colororder(clrs)
%
ax1      = axes('units','centimeters','position',ppos1);
N        = sqrt(R23.N2);
N(imag(N)~=0)=nan;
plot(ax1,datetime(R23.Time,'convertFrom','datenum'),N,'linewidth',1.5)
ylabel(ax1,'$N$ [s$^{-1}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex',...
   'xlim',datetime([release_times.StartTime_UTC_(2); release_times.EndTime_UTC_(3)]))
colororder(clrs)
grid on
%
ax2      = axes('units','centimeters','position',ppos2);
plot(ax2,datetime(R23.Time,'convertFrom','datenum'),R23.Temperature,'linewidth',1.5)
ylabel(ax2,'$T$ [$^\circ$C]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
        'xlim',datetime([release_times.StartTime_UTC_(2); release_times.EndTime_UTC_(3)]))
grid on
%
ax3      = axes('units','centimeters','position',ppos3);
plot(ax3,datetime(R23.Time,'convertFrom','datenum'),R23.Pressure(:,end)+R23.mab(end),'k','linewidth',1.5)
ylabel(ax3,'Depth [m]','interpreter','latex')
set(ax3,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
      'xlim',datetime([release_times.StartTime_UTC_(2); release_times.EndTime_UTC_(3)]))
grid on
%
cb       = axes('units','centimeters','position',cbpos);
imagesc([1:length(R23.mab)],0,reshape(clrs,1,length(R23.mab),3));
set(cb,'ytick',[],'xtick',[1:length(R23.mab)],'xticklabel',num2str(R23.mab','%1.2f'))
ylabel(cb,'[m]','rotation',0)
figname = [root_dir,'/data/Release2/figures/R23_Release_Stratification.pdf'];
exportgraphics(fig,figname)
%
%
%
%
% Compute release temperature PDF to compare with the observed
% define bins and load temperature time series
dbin  = 0.25;
Tbins = [12:dbin:20];
%
% need to extract temperature at level of dye relese!!!
T1    = R1.Temperature(:,3);
T2    = R2.Temperature(:,3);
T3    = R3.Temperature(:,3);
%
% generate histograms
pdf1    = histcounts(T1,[Tbins-dbin/2, Tbins(end)+dbin/2],'normalization','probability');
pdf2    = histcounts(T2,[Tbins-dbin/2, Tbins(end)+dbin/2],'normalization','probability');
pdf3    = histcounts(T3,[Tbins-dbin/2, Tbins(end)+dbin/2],'normalization','probability');
%
R1.Temperature_bins = Tbins;
R1.Temperature_pdf  = pdf1;
%
R2.Temperature_bins = Tbins;
R2.Temperature_pdf  = pdf2;
%
R3.Temperature_bins = Tbins;
R3.Temperature_pdf  = pdf3;
%
pw = 9;
ph = 4;
ppos= [xm ym pw ph];
ps = [2*xm+pw 2*ym+ph];
fig = figure('units','centimeters');
pos = get(fig,'position'); pos(3:4) = ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos);
plot(Tbins, pdf1,'-k',Tbins,pdf2,'-r',Tbins,pdf3,'-b','linewidth',2)
ll = legend('R1','R2','R3');
ylabel('$\mathrm{p}(T)$ [none]','interpreter','latex')
xlabel('temperature [C]','interpreter','latex')
figname = [root_dir,filesep,'figures',filesep,'Release_temp_pdf.png'];
exportgraphics(fig,figname)
%
% mean temperature of dye during release
T1dye = sum(Tbins.*pdf1)./sum(pdf1);
T2dye = sum(Tbins.*pdf2)./sum(pdf2);
T3dye = sum(Tbins.*pdf3)./sum(pdf3);
%
% mean stratification at release depth
dT1   = diff(R1.Temperature,1,2);
dTdz1 = dT1./diff(R1.mab,1,2);
dTdz1 = mean(dTdz1(2:3),2);
dT2   = diff(R2.Temperature,1,2);
dTdz2 = dT2./diff(R2.mab,1,2);
dTdz2 = mean(dTdz2(2:3),2);
dT3   = diff(R3.Temperature,1,2);
dTdz3 = dT3./diff(R3.mab,1,2);
dTdz3 = mean(dTdz3(2:3),2);
%
% estimate the standard plume (vertical) thickness using mean stratification
sigT1 = sqrt( sum( (Tbins-T1dye).^2.*pdf1, 'omitnan') ./ sum(pdf1,'omitnan') );
sigT2 = sqrt( sum( (Tbins-T2dye).^2.*pdf2, 'omitnan') ./ sum(pdf2,'omitnan') );
sigT3 = sqrt( sum( (Tbins-T3dye).^2.*pdf3, 'omitnan') ./ sum(pdf3,'omitnan') );
sigP1 = sigT1/dTdz1
sigP2 = sigT2/dTdz2
sigP3 = sigT3/dTdz3
%
fig = figure('units','centimeters');
pos = get(fig,'position'); pos(3:4) = ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos);
plot((Tbins-T1dye)/dTdz1, pdf1*dTdz1,'-k',(Tbins-T2dye)/dTdz2,pdf2*dTdz2,'-r',(Tbins-T3dye)/dTdz3,pdf3*dTdz3,'-b','linewidth',2)
ll = legend('R1','R2','R3');
ylabel('$\mathrm{p}(z'')$ [none]','interpreter','latex')
xlabel('z'' [m]','interpreter','latex')
set(ax1,'ticklabelinterpreter','latex','xlim',[-15 15])
figname = [root_dir,filesep,'figures',filesep,'Release_thickness_pdf.png'];
exportgraphics(gcf,figname)
%
% from: estimate_bulk_CTDF_stats_by_release.m
% $$$ dTdp_T1avg = -0.4797;
% $$$ dTdp_T2avg = -1.3773;
% $$$ dTdp_T3avg = -0.0468;
% $$$ sigP1 = 2.1914
% $$$ sigP2 = 0.8097
% $$$ sigP3 = 6.6642
archive_dir      = [root_dir,'/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat'];
save(archive_dir)