clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/Release2/raw';
inputFile = 'KELP2_Vector';
headingOffset = 181.31;% based on heading from AquodoppHR_KELP2
fileName  = [inputDir,'/',inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/Release2/raw';
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/Release2/L0';
L0Name  = [inputFile,'_L0'];
% Enter path to save figures
figDir = [inputDir,'/../figures/'];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-period for estimating the atmospheric pressure offset and deployment
atmosphTime = [datetime('08-Jul-2024 16:00:00'), datetime('08-Jul-2024 16:30:00')];
deployTime  = [datetime('08-Jul-2024 18:30:00'), datetime('11-Jul-2024 19:30:00')];
%
% time offset if necessary
tos = 0;
%
% returns structure A with all aquadopp data
load_VECTOR_data
save([outputDir,'/',outputName,'.mat'],'-struct','A')
%
date = datetime(A.sensor.date,'convertFrom','datenum');
time = datetime(A.sensor.date(1)+A.seconds/86400,'convertFrom','datenum');
%
if ~exist('atmosphTime','var')
% plot some stuff
disp('pick 2 points bounding when out of water for ATM pressure offset')
plot(pressure)
l = ginput(2);
atmosphTime = [time(l(1,1)) time(l(2,1))]
end
if ~exist('deployTime','var')
% now trim the data to when it was in the water
disp('pick start/end points of deployment')
plot(pressure)
l = ginput(2);
deployTime = [time(l(1,1)) time(l(2,1))]
end
iATM = find(time>=atmosphTime(1) & time<=atmosphTime(2));
iDEP = find(time>=deployTime(1) & time<=deployTime(2));
%    
A.pressureOffset = mean(pressure(iATM));
%
vel_valid = iDEP;
vars  = {'seconds','pressure','a1','a2','a3','v1','v2','v3','c1','c2','c3','SNR1','SNR2','SNR3','b1','b2','b3','east','north','up'};
for jj = 1:length(vars)
    eval(['A.',vars{jj},' = A.',vars{jj},'(vel_valid);'])
end
nsamples = length(vel_valid);
time     = time(vel_valid);
A.time   = datenum(time);
%
sen_valid = find(date>=deployTime(1) &...
                  date<=deployTime(2));
vars  = {'date','battery_voltage','sound_speed','heading','pitch','roll','temperature'};
for jj = 1:length(vars)
    eval(['A.',vars{jj},' = A.sensor.',vars{jj},'(sen_valid);'])
end
A.seconds = A.seconds-A.seconds(1);
date = date(sen_valid);
%
A.qcFlag   =  min(A.a1,min(A.a2,A.a3))>60 & min(A.c1,min(A.c2,A.c3))>50 & min(A.SNR1,min(A.SNR2,A.SNR3))>15;
%
% make a few quick plots
fig0 = figure;
ax1 = subplot(2,1,1);
plot(date,A.temperature)
ylabel(ax1,'$T$ [$^\circ$]','interpreter','latex')
set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
ax2 = subplot(2,1,2);
plot(time,A.pressure)
ylabel(ax2,'$P$ [m]','interpreter','latex')
xlabel(ax2,'time','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out','xlim',get(ax1,'xlim'))
figName = [figDir,'/',inputFile,'_temp_pres.png'];
exportgraphics(fig0,figName)
%
%
% use acceleration and jolt to filter bad data
utot = sqrt( A.v1.^2 + A.v2.^2 + A.v3.^2);
du   = gradient(utot)/dt;
ddu  = gradient(du)/dt;
%
r0   =  8*nanstd(utot);
r1   =  5*nanstd(du);
r2   =  4*nanstd(ddu);
%
valid = (utot/r0).^2 + (du/r1).^2 + (ddu/r2).^2 <= 1;
A.qcFlag = A.qcFlag & valid;
%
% need to apply this mask before conducting analysis
% make a few quick plots
fig0 = figure;
ax1 = subplot(3,1,1);
plot(time,A.east.*valid.*A.qcFlag,'.')
ylabel(ax1,'$v_\mathrm{\scriptscriptstyle{E}}$ [m/s]','interpreter','latex')
set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
ax2 = subplot(3,1,2);
plot(time,A.north.*valid.*A.qcFlag,'.')
ylabel(ax2,'$v_\mathrm{\scriptscriptstyle{N}}$ [m/s]','interpreter','latex')
xlabel(ax2,'time','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out','xlim',get(ax1,'xlim'))
ax3 = subplot(3,1,3);
plot(time,A.north.*valid.*A.qcFlag,'.')
ylabel(ax3,'$v_\mathrm{\scriptscriptstyle{U}}$ [m/s]','interpreter','latex')
xlabel(ax3,'time','interpreter','latex')
set(ax3,'ticklabelinterpreter','latex','tickdir','out','xlim',get(ax1,'xlim'))
figName = [figDir,'/',inputFile,'_velocity_ENU.png'];
exportgraphics(fig0,figName)
%
save([L0Dir,'/',L0Name,'.mat'],'-struct','A')
%
% add the config info to the structure A to quick save as netcdf4
fieldNames = fields(A.config);
originalFields = fields(A);
%
for j = 1:length(fieldNames)
 A.(fieldNames{j}) = A.config.(fieldNames{j});
end
A = orderfields(A,cat(1,fieldNames,originalFields));
ncfile = [L0Dir,'/',L0Name,'.nc'];
if exist(ncfile,'file')
    system(['rm ',ncfile]);
end
struct2nc(A,ncfile,'NETCDF4');
%
% $$$ %
% $$$ % plot some stuff
% $$$ disp('pick 2 points bounding when out of water for ATM pressure offset')
% $$$ plot(pressure)
% $$$ l = ginput(2);
% $$$ A.pressureOffset = mean(pressure(round(l(1,1)):round(l(2,1))));
% $$$ %
% $$$ % now trim the data to when it was in the water
% $$$ disp('pick start/end points of deployment')
% $$$ l = ginput(2);
% $$$ valid = [round(l(1,1)):round(l(2,1))]';
% $$$ vars  = {'seconds','pressure','a1','a2','a3','v1','v2','v3','c1','c2','c3','SNR1','SNR2','SNR3'};
% $$$ for jj = 1:length(vars)
% $$$     eval(['A.',vars{jj},' = A.',vars{jj},'(valid,:);'])
% $$$ end
% $$$ nsamples = length(valid);
% $$$ %
% $$$ date_valid = find(A.sensor.date>=A.sensor.date(1)+valid(1)/A.config.sample_rate/86400 &...
% $$$                   A.sensor.date<=A.sensor.date(1)+valid(end)/A.config.sample_rate/86400);
% $$$ vars  = {'date','battery_voltage','sound_speed','heading','pitch','roll','temperature'};
% $$$ for jj = 1:length(vars)
% $$$     eval(['A.sensor.',vars{jj},' = A.sensor.',vars{jj},'(date_valid);'])
% $$$ end
% $$$ A.seconds = A.seconds-A.seconds(1);
% $$$ %
% $$$ A.qcFlag   =  min(A.a1,min(A.a2,A.a3))>130 & min(A.c1,min(A.c2,A.c3))>75 & min(A.SNR1,min(A.SNR2,A.SNR3))>15;
% $$$ %
% $$$ date = datetime(A.sensor.date,'convertFrom','datenum');
% $$$ time = datetime(A.sensor.date(1)+A.seconds/86400,'convertFrom','datenum');
% $$$ %
% $$$ % make a few quick plots
% $$$ fig0 = figure;
% $$$ ax1 = subplot(2,1,1);
% $$$ plot(date,A.sensor.temperature)
% $$$ ylabel(ax1,'$T$ [$^\circ$]','interpreter','latex')
% $$$ set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
% $$$ ax2 = subplot(2,1,2);
% $$$ plot(time,A.pressure)
% $$$ ylabel(ax2,'$P$ [m]','interpreter','latex')
% $$$ xlabel(ax2,'time','interpreter','latex')
% $$$ set(ax2,'ticklabelinterpreter','latex','tickdir','out','xlim',get(ax1,'xlim'))
% $$$ figName = [figDir,'/',inputFile,'_temp_pres.png'];
% $$$ exportgraphics(fig0,figName)
% $$$ %
% $$$ %
% $$$ % $$$ % quick convolution running mean filter
% $$$ % $$$ np1 = round(0.3/A.config.binSize);% 10 cm vertical 
% $$$ % $$$ np2 = 31;% 5min for 1Hz data
% $$$ % $$$ f1 = hamming(np1);f1 = f1./sum(f1);
% $$$ % $$$ f2 = hamming(np2);f2 = f2./sum(f2);
% $$$ % $$$ % do a nan-mean filter, keep track of normalization
% $$$ % $$$ on = conv2(f1,f2,A.qcFlag','same');
% $$$ % $$$ %
% $$$ % $$$ A1 = conv2(f1,f2,(A.a1.*A.qcFlag)','same')./on;
% $$$ % $$$ A2 = conv2(f1,f2,(A.a2.*A.qcFlag)','same')./on;
% $$$ % $$$ A3 = conv2(f1,f2,(A.a3.*A.qcFlag)','same')./on;
% $$$ % $$$ %
% $$$ % $$$ fig1 = figure;
% $$$ % $$$ ax1 = subplot(3,1,1);
% $$$ % $$$ imagesc(time,A.dbins',A1.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'X')
% $$$ % $$$ set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ title(ax1,'Amplitude')
% $$$ % $$$ %
% $$$ % $$$ ax2 = subplot(3,1,2);
% $$$ % $$$ imagesc(time,A.dbins',A2.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'Y')
% $$$ % $$$ ylabel('mab','interpreter','latex')
% $$$ % $$$ set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ %
% $$$ % $$$ ax3 = subplot(3,1,3);
% $$$ % $$$ imagesc(time,A.dbins',A3.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'Z')
% $$$ % $$$ set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ xlabel('time [s]','interpreter','latex')
% $$$ % $$$ figName = [figDir,'/',inputFile,'_amplitude.png'];
% $$$ % $$$ exportgraphics(fig1,figName)
% $$$ % $$$ %
% $$$ % $$$ %
% $$$ % $$$ C1 = conv2(f1,f2,(A.c1.*A.qcFlag)','same')./on;
% $$$ % $$$ C2 = conv2(f1,f2,(A.c2.*A.qcFlag)','same')./on;
% $$$ % $$$ C3 = conv2(f1,f2,(A.c3.*A.qcFlag)','same')./on;
% $$$ % $$$ %
% $$$ % $$$ fig1 = figure;
% $$$ % $$$ ax1 = subplot(3,1,1);
% $$$ % $$$ imagesc(time,A.dbins',C1.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'X')
% $$$ % $$$ set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ title(ax1,'Correlation')
% $$$ % $$$ %
% $$$ % $$$ ax2 = subplot(3,1,2);
% $$$ % $$$ imagesc(time,A.dbins',C2.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'Y')
% $$$ % $$$ ylabel('mab','interpreter','latex')
% $$$ % $$$ set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ %
% $$$ % $$$ ax3 = subplot(3,1,3);
% $$$ % $$$ imagesc(time,A.dbins',C3.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'Z')
% $$$ % $$$ set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ xlabel('time [s]','interpreter','latex')
% $$$ % $$$ figName = [figDir,'/',inputFile,'_correlation.png'];
% $$$ % $$$ exportgraphics(fig1,figName)
% $$$ % $$$ %
% $$$ % $$$ %
% $$$ % $$$ % Now plot currents
% $$$ % $$$ V1 = conv2(f1,f2,(A.v1.*A.qcFlag)','same')./on;
% $$$ % $$$ V2 = conv2(f1,f2,(A.v2.*A.qcFlag)','same')./on;
% $$$ % $$$ V3 = conv2(f1,f2,(A.v3.*A.qcFlag)','same')./on;
% $$$ % $$$ %
% $$$ % $$$ fig2 = figure;
% $$$ % $$$ ax1 = subplot(3,1,1);
% $$$ % $$$ imagesc(time,A.dbins',V1.*qcFlag0'),caxis([-0.25 0.25]),colormap(cmocean('balance')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'X')
% $$$ % $$$ %
% $$$ % $$$ set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ ax2 = subplot(3,1,2);
% $$$ % $$$ imagesc(time,A.dbins',V2.*qcFlag0'),caxis([-0.25 0.25]),colormap(cmocean('balance')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'Y')
% $$$ % $$$ set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ ylabel('mab','interpreter','latex')
% $$$ % $$$ %
% $$$ % $$$ ax3 = subplot(3,1,3);
% $$$ % $$$ imagesc(time,A.dbins',V3.*qcFlag0'),caxis([-0.125 0.125]),colormap(cmocean('balance')),colorbar
% $$$ % $$$ text(time(1),ylims(2),'Z')
% $$$ % $$$ set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
% $$$ % $$$ xlabel('time [s]','interpreter','latex')
% $$$ % $$$ figName = [figDir,'/',inputFile,'_velocity.png'];
% $$$ % $$$ exportgraphics(fig2,figName)
% $$$ % $$$ %
% $$$ % $$$ % get the time-averaged current.
% $$$ % $$$ % Note, this is not in depth normalized (sigma) coordinates.
% $$$ % $$$ % first, nan any values that don't pass QC.
% $$$ % $$$ V1(~A.qcFlag')=nan;
% $$$ % $$$ V2(~A.qcFlag')=nan;
% $$$ % $$$ V3(~A.qcFlag')=nan;
% $$$ % $$$ flag = sum(A.qcFlag,1) > 0.50*nsamples;
% $$$ % $$$ % Uz = time averegd; U = depth & time averaged
% $$$ % $$$ Uz = nanmean(V1,2); U = nanmean(Uz); Uz(~flag)=nan;
% $$$ % $$$ Vz = nanmean(V2,2); V = nanmean(Vz); Vz(~flag)=nan;
% $$$ % $$$ Wz = nanmean(V3,2); W = nanmean(Wz); Wz(~flag)=nan;
% $$$ % $$$ %
% $$$ % $$$ fig3 = figure;
% $$$ % $$$ plot(Uz,A.dbins,'.-r',Vz,A.dbins,'.-b',Wz,A.dbins,'.-k')
% $$$ % $$$ hh = legend('X','Y','Z');
% $$$ % $$$ set(hh,'fontsize',9)
% $$$ % $$$ xlabel('m/s','interpreter','latex')
% $$$ % $$$ ylabel('mab','interpreter','latex')
% $$$ % $$$ set(gca,'ticklabelinterpreter','latex','tickdir','out')
% $$$ % $$$ figName = [figDir,'/',inputFile,'_mean_velocity.png'];
% $$$ % $$$ exportgraphics(fig3,figName)
% $$$ % $$$ %
% $$$ % $$$ A.VelX = V1;
% $$$ % $$$ A.VelY= V2;
% $$$ % $$$ A.VelZ   = V3;
% $$$ % $$$ %
% $$$ % $$$ A.VelXAvg = Uz;
% $$$ % $$$ A.VelYAvg= Vz;
% $$$ % $$$ A.VelZAvg   = Wz;
% $$$ % $$$ %
% $$$ % $$$ save([L0Dir,'/',L0Name,'.mat'],'-struct','A')
% $$$ % $$$ %
% $$$ % $$$ %
% $$$ % $$$ %
% $$$ % $$$ v1bar = nansum( A.v1.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ % $$$ v2bar = nansum( A.v2.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ % $$$ v3bar = nansum( A.v3.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ % $$$ % rotate 60 degrees:
% $$$ % $$$ vrot = [cosd(60) sind(60);-sind(60) cosd(60)]*[v1bar';v2bar'];
% $$$ % $$$ v1rot = vrot(1,:)';
% $$$ % $$$ v2rot = vrot(2,:)';
% $$$ % $$$ %
% $$$ % $$$ %
% $$$ % $$$ dv    = 0.01;
% $$$ % $$$ vbins = [-1:dv:1];
% $$$ % $$$ v1H   = hist(v1rot,vbins);
% $$$ % $$$ v2H   = hist(v2rot,vbins);
% $$$ % $$$ v3H   = hist(v3bar,vbins);
% $$$ % $$$ %
% $$$ % $$$ figure, plot(v1rot,v2rot,'.')
% $$$ % $$$ hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
% $$$ % $$$ xlabel('$\bar{v}_X$','interpreter','latex')
% $$$ % $$$ ylabel('$\bar{v}_Y$','interpreter','latex')
% $$$ % $$$ figName = [figDir,'/',inputFile,'_depth_avg_XYvel_checker_pattern.png'];
% $$$ % $$$ exportgraphics(gcf,figName)
% $$$ % $$$ %
% $$$ % $$$ %
% $$$ % $$$ figure,
% $$$ % $$$ subplot(3,1,1)
% $$$ % $$$ bar(vbins, v1H./length(v1bar))
% $$$ % $$$ ylabel('$p(\bar{v}_{X''})$','interpreter','latex')
% $$$ % $$$ subplot(3,1,2)
% $$$ % $$$ bar(vbins, v2H./length(v1bar))
% $$$ % $$$ ylabel('$p(\bar{v}''_{Y''})$','interpreter','latex')
% $$$ % $$$ subplot(3,1,3)
% $$$ % $$$ bar(vbins, v3H./length(v1bar))
% $$$ % $$$ ylabel('$p(\bar{v}_Z)$','interpreter','latex')
% $$$ % $$$ xlabel('velocity')
% $$$ % $$$ figName = [figDir,'/',inputFile,'_depth_avg_vel_XYrot60deg_histograms.png'];
% $$$ % $$$ exportgraphics(gcf,figName)
