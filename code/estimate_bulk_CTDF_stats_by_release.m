clear all
close all
%
rootDir  = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/'
%
% 0) load CTDF data from each release
releaseNumber = 1;
CTD1 = load_L0_CTDF_casts(releaseNumber);
CTD2 = load_L0_CTDF_casts(2);
CTD3 = load_L0_CTDF_casts(3);
%
N1 = size(CTD1.dye_grid,2);
D1 = CTD1.dye_grid(:); T1 = CTD1.temp_grid(:); P1 = repmat(CTD1.pres_grid,N1,1); 
P1dye = sum( D1(D1>=1).*P1(D1>=1), 'omitnan')./sum(D1(D1>=1),'omitnan');
T1dye = sum( D1(D1>=1).*T1(D1>=1), 'omitnan')./sum(D1(D1>=1),'omitnan');
%
N2 = size(CTD2.dye_grid,2);
D2 = CTD2.dye_grid(:); T2 = CTD2.temp_grid(:); P2 = repmat(CTD2.pres_grid,N2,1); 
P2dye = sum( D2(D2>=1).*P2(D2>=1), 'omitnan')./sum(D2(D2>=1),'omitnan');
T2dye = sum( D2(D2>=1).*T2(D2>=1), 'omitnan')./sum(D2(D2>=1),'omitnan');
%
N3 = size(CTD3.dye_grid,2);
D3 = CTD3.dye_grid(:); T3 = CTD3.temp_grid(:); P3 = repmat(CTD3.pres_grid,N3,1); 
P3dye = sum( D3(D3>=1).*P3(D3>=1), 'omitnan')./sum(D3(D3>=1),'omitnan');
T3dye = sum( D3(D3>=1).*T3(D3>=1), 'omitnan')./sum(D3(D3>=1),'omitnan');
%
[~,srt1] = sort(D1,'ascend');
[~,srt2] = sort(D2,'ascend');
[~,srt3] = sort(D3,'ascend');
%
clims     = [0 10];
[cm,conc] = make_DYE_colormap(clims,1);
%
xm = 2.5;
ym = 2.5;
pw = 4;
ph = 9;
ag = 0.5;
ppos1 = [xm ym pw ph];
ppos2 = [xm+pw+ag ym pw ph];
ppos3 = [xm+2*pw+2*ag ym pw ph];
cbpos = [xm+pw-1.5*ag ym+0.5*ag ag ph/3];
ps    = [1.5*xm+3*pw+2*ag 2*ym+ph]
fig = figure('units','centimeters');
pos = get(fig,'position'); pos(3:4) = ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos1);
scatter(ax1,T1(srt1),P1(srt1),15,log2(D1(srt1)),'filled'), set(ax1,'ydir','reverse','ticklabelinterpreter','latex','xlim',[14 20],'ylim',[0 15],'tickdir','out'), colormap(ax1,cm), caxis(ax1,clims)
grid on
xline(ax1,T1dye,'--k','linewidth',1)
yline(ax1,P1dye,'--k','linewidth',1)
ylabel(ax1,' depth [m] ','interpreter','latex')
title('Release 1','interpreter','latex')
%
ax2 = axes('units','centimeters','position',ppos2);
scatter(ax2,T2(srt2),P2(srt2),15,log2(D2(srt2)),'filled'), set(ax2,'ydir','reverse','ticklabelinterpreter','latex','yticklabel',[],'xlim',[14 20],'ylim',[0 15],'tickdir','out'), colormap(ax2,cm), caxis(ax2,clims)
grid on
xline(ax2,T2dye,'--k','linewidth',1)
yline(ax2,P2dye,'--k','linewidth',1)
xlabel(ax2,' temperature [C] ','interpreter','latex')
title('Release 2','interpreter','latex')
%
ax3 = axes('units','centimeters','position',ppos3);
scatter(ax3,T3(srt3),P3(srt3),15,log2(D3(srt3)),'filled'), set(ax3,'ydir','reverse','ticklabelinterpreter','latex','yticklabel',[],'xlim',[14 20],'ylim',[0 15],'tickdir','out'), colormap(ax3,cm), caxis(ax3,clims)
grid on
xline(ax3,T3dye,'--k','linewidth',1)
yline(ax3,P3dye,'--k','linewidth',1)
title('Release 3','interpreter','latex')
%
cb0=axes('units','centimeters','position',cbpos);
imagesc(0,conc,reshape(cm,256,1,3))
ticks = get(cb0,'ytick');
set(cb0,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex','xaxislocation','top','ydir','normal','xtick',[])
xlabel(cb0,{'Dye';' [ppb]'},'interpreter','latex')
figdir  = [rootDir,filesep,'figures',filesep];
if ~exist(figdir,'dir')
    eval(['!mkdir -p ',figdir])
end
figname = [figdir,'CTDF_depth_temp_scatter_dye_colormap.png'];
exportgraphics(gcf,figname)
%
%
% plot 2D histograms
nbins = [];,
Tbins = [12:0.25:20];
Pbins = [0:0.125:15];
[T1out, P1out, M1, S1] = bin_stats2D(T1,P1,D1,nbins,Tbins,Pbins);
[T2out, P2out, M2, S2] = bin_stats2D(T2,P2,D2,nbins,Tbins,Pbins);
[T3out, P3out, M3, S3] = bin_stats2D(T3,P3,D3,nbins,Tbins,Pbins);
%
Q1 = sum(M1(:),'omitnan');
m1 = sum(M1,1,'omitnan')/Q1;
%
Q2 = sum(M2(:),'omitnan');
m2 = sum(M2,1,'omitnan')/Q2;
%
Q3 = sum(M3(:),'omitnan');
m3 = sum(M3,1,'omitnan')/Q3;
%
pw = 9;
ph = 4;
ppos= [xm ym pw ph];
ps = [2*xm+pw 2*ym+ph];
fig = figure('units','centimeters');
pos = get(fig,'position'); pos(3:4) = ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos);
plot(Tbins, m1,'-k',Tbins,m2,'-r',Tbins,m3,'-b','linewidth',2)
ll = legend('R1','R2','R3');
ylabel('$\mathrm{p}(T)$ [none]','interpreter','latex')
xlabel('temperature [C]','interpreter','latex')
figname = [figdir,'CTDF_temp_pdf.png'];
exportgraphics(gcf,figname)
%
%
% now estimate dT/dp at the mean temperature
dT1 = gradient(T1);
dP1 = gradient(P1);dP1(dP1<0)=nan;
% average within 0.1 C
dTdp1= dT1./(sign(dP1).*max(abs(dP1),0.01));
dTdp_T1avg= mean( dTdp1(T1>=T1dye-0.1 & T1<=T1dye+0.1),'omitnan');
% 
dT2 = gradient(T2);
dP2 = gradient(P2);dP2(dP2<0)=nan;
dTdp2= dT2./(sign(dP2).*max(abs(dP2),0.01));
dTdp_T2avg= mean( dTdp2(T2>=T2dye-0.1 & T2<=T2dye+0.1),'omitnan');
%
dT3 = gradient(T3);
dP3 = gradient(P3);dP3(dP3<0)=nan;
dTdp3= dT3./(sign(dP3).*max(abs(dP3),0.01));
dTdp_T3avg= mean( dTdp3(T3>=T3dye-0.1 & T3<=T3dye+0.1),'omitnan');
%
Tlow = Tbins - 0.5*gradient(Tbins);
Thigh= Tbins + 0.5*gradient(Tbins);
for jj=1:length(Tbins)
    in1 = find(T1>=Tlow(jj) & T1<Thigh(jj));
    in2 = find(T2>=Tlow(jj) & T2<Thigh(jj));
    in3 = find(T3>=Tlow(jj) & T3<Thigh(jj));    
    strat1(jj) = mean(dTdp1(in1),'omitnan');
    strat2(jj) = mean(dTdp2(in2),'omitnan');
    strat3(jj) = mean(dTdp3(in3),'omitnan');    
end
%
% where was release relative to maximum stratification?
% recreate the time/space mean temperature profile
dp = min(abs((Tbins(2)-Tbins(1))./min(strat1)),abs((Tbins(2)-Tbins(1))./min(strat2)));
dp = min(dp,abs((Tbins(2)-Tbins(1))./min(strat3)));
dp = min(dp, 0.01);
%
fig = figure('units','centimeters');
pos = get(fig,'position'); pos(3:4) = ps; set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos);
plot(-(Tbins-T1dye)/dTdp_T1avg, m1*-dTdp_T1avg,'-k',-(Tbins-T2dye)/dTdp_T2avg,m2*-dTdp_T2avg,'-r',-(Tbins-T3dye)/dTdp_T3avg,m3*-dTdp_T3avg,'-b','linewidth',2)
ll = legend('R1','R2','R3');
ylabel('$\mathrm{p}(z'')$ [none]','interpreter','latex')
xlabel('z'' [m]','interpreter','latex')
set(ax1,'ticklabelinterpreter','latex','xlim',[-15 15])
figname = [figdir,'CTDF_thickness_pdf.png'];
exportgraphics(gcf,figname)
%
% estimate standard thickness
sigT1 = sqrt( sum( (Tbins-T1dye).^2.*m1, 'omitnan') ./ sum(m1,'omitnan') );
sigT2 = sqrt( sum( (Tbins-T2dye).^2.*m2, 'omitnan') ./ sum(m2,'omitnan') );
sigT3 = sqrt( sum( (Tbins-T3dye).^2.*m3, 'omitnan') ./ sum(m3,'omitnan') );
sigP1 = -sigT1/dTdp_T1avg
sigP2 = -sigT2/dTdp_T2avg
sigP3 = -sigT3/dTdp_T3avg
%
% 1) for each release, estimate the sigma2_x, sigma2_y, sigma_xy,
D1sumz = CTD1.dye_grid; D1sumz(P1<1 | D1sumz(:)<1) = nan;
D1sumz = sum(D1sumz,1,'omitnan');
lat1  = CTD1.latitude_grid;
lon1  = CTD1.longitude_grid;
%
D2sumz = CTD2.dye_grid; D2sumz(P2<1 | D2sumz(:)<1) = nan;
D2sumz = sum(D2sumz,1,'omitnan');
lat2  = CTD2.latitude_grid;
lon2  = CTD2.longitude_grid;
%
D3sumz = CTD3.dye_grid; D3sumz(P3<1 | D3sumz(:)<1) = nan;
D3sumz = sum(D3sumz,1,'omitnan');
lat3  = CTD3.latitude_grid;
lon3  = CTD3.longitude_grid;
%
angle = -70;
[x1,y1] = lltoxy_AQR(lat1,lon1,angle);
[x2,y2] = lltoxy_AQR(lat2,lon2,angle);
[x3,y3] = lltoxy_AQR(lat3,lon3,angle);
%
mX1 = sum( D1sumz.*x1,'omitnan')./sum(D1sumz,'omitnan');
mX2 = sum( D2sumz.*x2,'omitnan')./sum(D2sumz,'omitnan');
mX3 = sum( D3sumz.*x3,'omitnan')./sum(D3sumz,'omitnan');
%
mY1 = sum( D1sumz.*y1,'omitnan')./sum(D1sumz,'omitnan');
mY2 = sum( D2sumz.*y2,'omitnan')./sum(D2sumz,'omitnan');
mY3 = sum( D3sumz.*y3,'omitnan')./sum(D3sumz,'omitnan');
%
sX1 = sum( D1sumz.*(x1-mX1).^2,'omitnan')./sum(D1sumz,'omitnan');
sX2 = sum( D2sumz.*(x2-mX2).^2,'omitnan')./sum(D2sumz,'omitnan');
sX3 = sum( D3sumz.*(x3-mX3).^2,'omitnan')./sum(D3sumz,'omitnan');
%
sXY1 = sum( D1sumz.*(x1-mX1).*(y1-mY1),'omitnan')./sum(D1sumz,'omitnan');
sXY2 = sum( D2sumz.*(x2-mX1).*(y2-mY2),'omitnan')./sum(D2sumz,'omitnan');
sXY3 = sum( D3sumz.*(x3-mX3).*(y3-mY3),'omitnan')./sum(D3sumz,'omitnan');
%
sY1 = sum( D1sumz.*(y1-mY1).^2,'omitnan')./sum(D1sumz,'omitnan');
sY2 = sum( D2sumz.*(y2-mY2).^2,'omitnan')./sum(D2sumz,'omitnan');
sY3 = sum( D3sumz.*(y3-mY3).^2,'omitnan')./sum(D3sumz,'omitnan');
fprintf ('\n absolute plume scales: sigma_x, sigma_xy, sigma_y \n')
[sX1, sXY1, sY1;...
 sX2, sXY2, sY2;...
 sX3, sXY3, sY3]
%
% rotate these...
[V1,S1] = eig([sX1 sXY1; sXY1 sY1]); theta1 = 90-atan2d((V1(1)), (V1(1,2)));
[V2,S2] = eig([sX2 sXY2; sXY2 sY2]); theta2 = 90-atan2d((V2(1)), (V2(2,1)));
[V3,S3] = eig([sX3 sXY3; sXY3 sY3]); theta3 = 90-atan2d((V3(1)), (V3(2,1)));
%
fprintf ('\n absolute plume scales: sigma_x'', sigma_y'', theta \n')
S1 = sqrt(S1);
S2 = sqrt(S2);
S3 = sqrt(S3);
[S1(1), S1(4), theta1;...
 S2(1), S2(4), theta2;...
 S3(1), S3(4), theta3]
%
figure,
plot(-dTdp_T1avg,sqrt(S1(1)*S1(4)),'*k',-dTdp_T2avg,sqrt(S2(1)*S2(4)),'*r',-dTdp_T3avg,sqrt(S3(1)*S3(4)),'*b',...
     -dTdp_T1avg,S1(1),'^k',-dTdp_T2avg,S2(1),'xr',-dTdp_T3avg,S3(1),'xb',...
     -dTdp_T1avg,S1(4),'xk',-dTdp_T2avg,S2(4),'^r',-dTdp_T3avg,S3(4),'^b')
leg = legend( '$\sqrt{\lambda_x\lambda_y}$','$\lambda_x$','$\lambda_y$');
set(leg,'interpreter','latex')
xlabel('$dT/dz$ [C/m]','interpreter','latex')
ylabel('$\lambda_{x,y}$ [m]','interpreter','latex')
%
figure,
plot(-dTdp_T1avg,sigP1,'*k',-dTdp_T2avg,sigP2,'*r',-dTdp_T3avg,sigP3,'*b')
xlabel('$dT/dz$ [C/m]','interpreter','latex')
ylabel('$\lambda_z$ [m]','interpreter','latex')
%
figure,
releaseLatitude  = 34.4693000;
releaseLongitude = -120.1277833;
% plot plume lat/lon
t = [0:0.01:2*pi];
x = cos(t);
y = sin(t);
X1 = mX1 + S1(1)*x*cosd(theta1) - S1(4)*y*sind(theta1);
Y1 = mY1 + S1(1)*x*sind(theta1) + S1(4)*y*cosd(theta1);
X2 = mX2 + S2(1)*x*cosd(theta2) - S2(4)*y*sind(theta2);
Y2 = mY2 + S2(1)*x*sind(theta2) + S2(4)*y*cosd(theta2);
X3 = mX3 + S3(1)*x*cosd(theta3) - S3(4)*y*sind(theta3);
Y3 = mY3 + S3(1)*x*sind(theta3) + S3(4)*y*cosd(theta3);
[Lat1,Lon1] = xytoll_AQR(X1,Y1,angle);
[Lat2,Lon2] = xytoll_AQR(X2,Y2,angle);
[Lat3,Lon3] = xytoll_AQR(X3,Y3,angle);
%
%
% get sensor lat/lon1
fin = '~/git/kelp/info/mooring_lat_lon.csv';
dat = readtable(fin);
Lat = table2array(dat(:,2));
Lon = table2array(dat(:,3));
geoplot(Lat1,Lon1,'k',Lat2,Lon2,'r',Lat3,Lon3,'b'),
hold on, geoplot(Lat,Lon,'ok','markerfacecolor','m','markersize',10)
hold on, geoplot(releaseLatitude,releaseLongitude,'ok',releaseLatitude,releaseLongitude,'xk','markersize',20,'linewidth',4,'markerfacecolor','m')
geobasemap('satellite')
geolimits([34.4667 34.4718],[-120.1324 -120.1246])
% $$$ ax = axes('position',[0.8 0.3 0.03 0.2]);
% $$$ imagesc(0,2.^conc,reshape(cm,256,1,3))
% $$$ ylabel('$\bar{D}$ [ppb]','interpreter','latex')
% $$$ set(ax,'ydir','normal','xtick',[],'ticklabelinterpreter','latex','xcolor','w','ycolor','w','fontsize',16)
figname = [figdir,'Release_plume_scales.png'];
exportgraphics(gcf,figname)
%
%
archive_dir      = [rootDir,'/data/2024_PROCESSED_DATA/VesselCTDFData.mat'];
save(archive_dir)
% 2) for each release, estimate the vertical and thermal width of the plume.
% 3) did we observe the plume relatively soon or late after release start? I think all relatively similar
