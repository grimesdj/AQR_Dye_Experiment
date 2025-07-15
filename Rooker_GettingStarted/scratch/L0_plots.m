


function L0_plots(A, inputFile, figDir)

% make a few quick plots
fig0 = figure;
ax1 = subplot(2,1,1);
plot(A.time,A.temperature)
ylabel(ax1,'$T$ [$^\circ$]','interpreter','latex')
set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
ax2 = subplot(2,1,2);
plot(A.time,A.pressure)
ylabel(ax2,'$P$ [m]','interpreter','latex')
xlabel(ax2,'time [s]','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out')
figName = [figDir,'/',inputFile,'_temperature_pressure.png'];
exportgraphics(fig0,figName)
%
%

fig1 = figure;
ax1 = subplot(3,1,1);
imagesc(A.time,A.dbins',A.A1.*A.qcFlag'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(A.time(1),A.ylims(2),'X')
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
title(ax1,'Amplitude')
%
ax2 = subplot(3,1,2);
imagesc(A.time,A.dbins',A.A2.*A.qcFlag'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(A.time(1),A.ylims(2),'Y')
ylabel('mab','interpreter','latex')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
%
ax3 = subplot(3,1,3);
imagesc(A.time,A.dbins',A.A3.*A.qcFlag'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(A.time(1),A.ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir,'/',inputFile,'_amplitude.png'];
exportgraphics(fig1,figName)
%

fig1 = figure;
ax1 = subplot(3,1,1);
imagesc(A.time,A.dbins',A.C1.*A.qcFlag'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(A.time(1),A.ylims(2),'X')
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
title(ax1,'Correlation')
%
ax2 = subplot(3,1,2);
imagesc(A.time,A.dbins',A.C2.*A.qcFlag'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(A.time(1),A.ylims(2),'Y')
ylabel('mab','interpreter','latex')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
%
ax3 = subplot(3,1,3);
imagesc(A.time,A.dbins',A.C3.*A.qcFlag'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(A.time(1),A.ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir,'/',inputFile,'_correlation.png'];
exportgraphics(fig1,figName)
%
%

fig2 = figure;
ax1 = subplot(3,1,1);
imagesc(A.time,A.dbins',A.PCA_X),caxis([-0.25 0.25]),colormap(cmocean('balance')),colorbar
text(A.time(1),A.ylims(2),'X')
%
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
ax2 = subplot(3,1,2);
imagesc(A.time,A.dbins',A.PCA_Y),caxis([-0.25 0.25]),colormap(cmocean('balance')),colorbar
text(A.time(1),A.ylims(2),'Y')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
ylabel('mab','interpreter','latex')
%
ax3 = subplot(3,1,3);
imagesc(A.time,A.dbins',A.VelZ.*A.qcFlag'),caxis([-0.125 0.125]),colormap(cmocean('balance')),colorbar
text(A.time(1),A.ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',A.ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir,'/',inputFile,'_velocity.png'];
exportgraphics(fig2,figName)
%

%
fig3 = figure;
plot(A.VelXAvg,A.dbins,'.-r',A.VelYAvg,A.dbins,'.-b',A.VelZAvg,A.dbins,'.-k')
hh = legend('X','Y','Z');
set(hh,'fontsize',9)
xlabel('m/s','interpreter','latex')
ylabel('mab','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
figName = [figDir,'/',inputFile,'_mean_velocity.png'];
exportgraphics(fig3,figName)

v1bar = nansum( A.v1.*A.qcFlag ,2)./sum(A.qcFlag,2);
v2bar = nansum( A.v2.*A.qcFlag ,2)./sum(A.qcFlag,2);
v3bar = nansum( A.v3.*A.qcFlag ,2)./sum(A.qcFlag,2);
%rotate 60 degrees:
vrot = [cosd(60) sind(60);-sind(60) cosd(60)]*[v1bar';v2bar'];
v1rot = vrot(1,:)';
v2rot = vrot(2,:)';
%
%

dv    = 0.01;
vbins = [-1:dv:1];
v1H   = hist(v1rot,vbins);
v2H   = hist(v2rot,vbins);
v3H   = hist(v3bar,vbins);
%
figure, plot(v1rot,v2rot,'.')
hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
xlabel('$\bar{v}_X$','interpreter','latex')
ylabel('$\bar{v}_Y$','interpreter','latex')
figName = [figDir,'/',inputFile,'_depth_avg_XYvel_checker_pattern.png'];
exportgraphics(gcf,figName)
%
%
figure,
subplot(3,1,1)
bar(vbins, v1H./length(v1bar))
ylabel('$p(\bar{v}_{X''})$','interpreter','latex')
subplot(3,1,2)
bar(vbins, v2H./length(v1bar))
ylabel('$p(\bar{v}''_{Y''})$','interpreter','latex')
subplot(3,1,3)
bar(vbins, v3H./length(v1bar))
ylabel('$p(\bar{v}_Z)$','interpreter','latex')
xlabel('velocity')
figName = [figDir,'/',inputFile,'_depth_avg_vel_XYrot60deg_histograms.png'];
exportgraphics(gcf,figName)
%
%
%
%
b1bar = nansum( A.b1.*A.qcFlag ,2)./sum(A.qcFlag,2);
b2bar = nansum( A.b2.*A.qcFlag ,2)./sum(A.qcFlag,2);
b3bar = nansum( A.b3.*A.qcFlag ,2)./sum(A.qcFlag,2);
figure, plot(b1bar,b2bar,'.')
hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
xlabel('$\bar{v}_X$','interpreter','latex')
ylabel('$\bar{v}_Y$','interpreter','latex')
figName = [figDir,'/',inputFile,'_depth_avg_Beam1Beam2_checker_pattern.png'];
exportgraphics(gcf,figName)




end