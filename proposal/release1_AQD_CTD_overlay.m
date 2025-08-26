clear all
close all
fin = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/Summer2025/Rooker/Release1/L0/KELP1_AquadoppHR_L0.mat';
aqd = load(fin);
%
% do some smoothing
dt_aqd = (aqd.Time(2)-aqd.Time(1))*86400;
dz_aqd = aqd.Config.binSize;
zbins  = aqd.Config.blank+[1:aqd.Config.Nbins]*aqd.Config.binSize;
%
Nt     = round(5*60/dt_aqd); if iseven(Nt), Nt=Nt+1; end
Nz     = round(0.1/dz_aqd);  if iseven(Nz), Nz=Nz+1; end
%
ft     = hamming(Nt); ft = ft./sum(ft);
fz     = hamming(Nz); fz = fz./sum(fz);
%
tmpU    = aqd.Velocity_East';
tmpV    = aqd.Velocity_North';
tmpQC   =~isnan(tmpU);
tmpU(~tmpQC)=0;
tmpV(~tmpQC)=0;
%
% rotate...
rotUV = [cosd(-45) sind(-45); -sind(-45), cosd(-45)]*[tmpU(:)';tmpV(:)'];
rotU  = reshape(rotUV(1,:),size(tmpU));
rotV  = reshape(rotUV(2,:),size(tmpV));
%
on      = conv2(fz,ft,tmpQC,'same');
Uavg    = conv2(fz,ft,rotU,'same')./on;
Vavg    = conv2(fz,ft,rotV,'same')./on;
%
% remove ends
Uavg    = Uavg(Nz:end-Nz,Nt:end-Nt);
Vavg    = Vavg(Nz:end-Nz,Nt:end-Nt);
time    = aqd.Time(Nt:end-Nt);
zbins   = zbins(Nz:end-Nz);
%
%
% make the velocity plot
xm = 2;
ym = 2;
pw = 8;
ph = 3.5;
ag = 0.5;
%
ppos1 = [xm ym pw ph];
ppos2 = [xm ym+ag+ph pw ph];
cbpos1= [xm+pw+ag ym ag ph*0.75];
cbpos2= [xm+pw+ag ym+ag+ph ag ph*0.75];
ps    = [2*xm+pw+2*ag 2*ym+2*ph+ag];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
% x-tick locations/times
xnum = datenum('July-03-2024 19:00:00') + [0:1:3]*3600/86400;
xstr = datestr(xnum,'HH:MM');
%
a1 = axes('units','centimeters','position',ppos1);
imagesc(time,zbins,Vavg*100),caxis(100*[-0.1 0.1]), colormap(cmocean('balance'))
set(a1,'ticklabelinterpreter','latex','ydir','normal','tickdir','out','xtick',xnum,'xticklabel',xstr)

annotation('textbox','units','centimeters','Position',[ppos1(1)+0.05 ppos1(2)+0.875*ppos1(4) ph/2 xm/4],'string','b) Alongshore','interpreter','latex','edgecolor','none','fitboxtotext','on','rotation',0,'fontsize',12)

% $$$ c1 = axes('units','centimeters','position',cbpos1);
% $$$ imagesc(0,100*[-0.1:0.2/255:0.1],reshape(cmocean('balance'),256,1,3))
% $$$ xlabel(c1,'~~~~~[cm/s]','interpreter','latex')
% $$$ set(c1,'ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','top','xtick',[],'tickdir','out')

a2 = axes('units','centimeters','position',ppos2);
imagesc(time,zbins,100*Uavg),caxis(100*[-0.1 0.1]), colormap(cmocean('balance'))
set(a2,'ticklabelinterpreter','latex','ydir','normal','tickdir','out','xticklabel',[],'xtick',xnum)

annotation('textbox','units','centimeters','Position',[ppos2(1)+0.05 ppos2(2)+0.875*ppos2(4) ph/2 xm/4],'string','a) Cross-shore','interpreter','latex','edgecolor','none','fitboxtotext','on','rotation',0,'fontsize',12)

c2 = axes('units','centimeters','position',cbpos2);
imagesc(0,100*[-0.1:0.2/255:0.1],reshape(cmocean('balance'),256,1,3))
xlabel(c2,'~~~~~$u$ [cm/s]','interpreter','latex')
set(c2,'ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','top','xtick',[],'tickdir','out','ydir','normal')

% ylabel([a1 a2],'meters above bottom','interpreter','latex')
annotation('textbox','units','centimeters','Position',[xm*2/3 ym+.5*ph ph/2 xm/4 ],'string','meters above bottom','interpreter','latex','edgecolor','none','fitboxtotext','on','rotation',90)
%
%
%
%
% 1) now load release CTD data
fin = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat';
load(fin,'ADCP1')
Tmap = [13:0.5:15.5];
clrs = cmocean('thermal',length(Tmap));
%
ctdTime = ADCP1.Time';
dt_ctd  = (ctdTime(2)-ctdTime(1))*86400;
ctdTemp = ADCP1.Temperature;
%
Nt     = round(5*60/dt_ctd); if iseven(Nt), Nt=Nt+1; end
ft     = hamming(Nt); ft = ft./sum(ft);
ctdTemp = conv2(ctdTemp,ft,'same');
ctdTemp = ctdTemp(Nt:end-Nt,:)';
ctdTime = ctdTime(Nt:end-Nt);
ctdZ    = ADCP1.mab';
% $$$ tmp = figure;
% $$$ [cont, dum] = contour(ADCP1.Time',ADCP1.mab',ADCP1.Temperature',Tmap);
% $$$ close(tmp);
% $$$ idx = 1;
% $$$ for jj  = 1:length(Tmap)
% $$$ 
% $$$     lvl = cont(1,idx);
% $$$     num = cont(2,idx)
% $$$     tt  = cont(1,idx+[1:num]);
% $$$     zz  = cont(2,idx+[1:num]);
% $$$     idx = idx+num+1;
% $$$     
% $$$     hold(a1,'on'), plot(a1,datetime(tt,'convertfrom','datenum'),zz,'o-','color',clrs(jj,:),'markerfacecolor',clrs(jj,:),'linewidth',2)
% $$$ 
% $$$     
% $$$     hold(a2,'on'), plot(a2,datetime(tt,'convertfrom','datenum'),zz,'o-','color',clrs(jj,:),'markerfacecolor',clrs(jj,:),'linewidth',2)
% $$$ 
% $$$ end


for jj  = 1:length(Tmap)
    hold(a1,'on'), contour(a1,ctdTime,ctdZ,ctdTemp,Tmap(jj)*[1 1],'linecolor',clrs(jj,:),'linewidth',2.5)
    hold(a2,'on'), contour(a2,ctdTime,ctdZ,ctdTemp,Tmap(jj)*[1 1],'linecolor',clrs(jj,:),'linewidth',2.5)
end


c1 = axes('units','centimeters','position',cbpos1);
imagesc(0,Tmap,reshape(clrs,length(Tmap),1,3))
xlabel(c1,'~~~~~[$^\circ$C]','interpreter','latex')
set(c1,'ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','top','xtick',[],'tickdir','out','ydir','normal')

figname = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/figures/Release1/Release_ADCP_CTD_overlay.pdf';
exportgraphics(fig,figname)


% 2) map T(z,t) to z(T,t)
% $$$ Tmap = [13:0.5:15.5];
% $$$ [tt,TT] = meshgrid(time,Tmap');
% $$$ [tt0,zz0] = meshgrid(ADCP1.Time',ADCP1.mab');
% $$$ zmap    = griddata(tt0,ADCP1.Temperature',zz0,tt,TT,'linear');
% $$$ % 3) contour on above axes
% $$$ clrs    = cmocean('thermal',length(Tmap));
% $$$ for jj  = 1:length(Tmap)
% $$$ 
% $$$     axes(a1)
% $$$     hold on, plot(time,zmap(jj,:),'linewidth',2,'color',clrs(jj,:))
% $$$ 
% $$$     axes(a2)
% $$$     hold on, plot(time,zmap(jj,:),'linewidth',2,'color',clrs(jj,:))
% $$$ 
% $$$ end