clear all
% full path
releaseNumber    = 2;
releaseLatitude  = 34.4693000;
releaseLongitude = -120.1277833;
releaseTime      = [datenum('08-Jul-18:41:00'), datenum('08-Jul-20:32:00')];
ctdSerialNumber  = 2348;
rootDir  = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/'
dataDir  = sprintf('%s/data/Release%d/',rootDir,releaseNumber);
L0Dir    = [dataDir,filesep,'L0',];
fileRoot = sprintf('%d*.mat',ctdSerialNumber);
files    = dir([L0Dir,filesep,fileRoot]);
% preallocate
time = [];
pres = [];
dye_grid = [];
temp_grid = [];
latitude_grid = [];
longitude_grid = [];
for jj  = 1:length(files)
    SN  = split(files(jj).name,'_');
    fin= [L0Dir,filesep,SN{1},'_L0.mat'];
    data = load(fin);
    %
    time      = cat(2,time,data.time_grid);
    temp_grid = cat(2,temp_grid,data.temp_grid);
    dye_grid = cat(2,dye_grid,data.dye_grid);
    latitude_grid = cat(2,latitude_grid,data.latitude_grid);
    longitude_grid = cat(2,longitude_grid,data.longitude_grid);    
    if jj==1
        pres_grid   = data.pres_grid;
    end
end
    time_grid = repmat(time,size(temp_grid,1),1);
    pres_grid = repmat(pres_grid,1,size(temp_grid,2));
    dye_grid(pres_grid<1) = 0;
    %
    Tdye = sum(temp_grid(:).*dye_grid(:),'omitnan')./sum(dye_grid(:),'omitnan')
    Pdye = sum(pres_grid(:).*dye_grid(:),'omitnan')./sum(dye_grid(:),'omitnan')    
    %
% $$$     figure, scatter(dye,pres, 10,temp,'filled')
% $$$     figure, plot(temp,dye_raw,'.')
    %
    [dye_sort,srt] = sort(dye_grid(:),'ascend');
    temp_sort = temp_grid(srt);
    pres_sort = pres_grid(srt);
    [cm,conc] = make_DYE_colormap([0 8],1);
    %    
    figure,
    scatter(temp_sort,pres_sort,15,log2(dye_sort),'filled'), set(gca,'ydir','reverse','ticklabelinterpreter','latex'), colormap(cm), caxis([0 8])
    xline(Tdye,'--r','linewidth',2)
    xlabel('Temp [C]','interpreter','latex')
    ylabel('Depth [m]','interpreter','latex')
    cb0=colorbar;
    ticks = get(cb0,'ytick');
    set(cb0,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb0,'Dye [ppb]','interpreter','latex')
    figdir  = [rootDir,filesep,'figures',filesep,sprintf('Release%d',releaseNumber),filesep];
    if ~exist(figdir,'dir')
        eval(['!mkdir -p ',figdir])
    end
    figname = [figdir,sprintf('Release%d_depth_temp_scatter_dye_colormap.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %    
    %
    figure,
    ax1 = subplot(2,1,1);
    scatter(datetime(time_grid(:),'convertfrom','datenum'),pres_grid(:),20,log2(dye_grid(:)),'filled'),
    set(ax1,'ydir','reverse','ticklabelinterpreter','latex'), colormap(ax1,cm), caxis(ax1,[0 8]),
    ylabel('Depth [m]','interpreter','latex')
    cb1=colorbar;
    ticks = get(cb1,'ytick');
    set(cb1,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb1,'Dye [ppb]','interpreter','latex')
    yline(ax1,Pdye,'--r','linewidth',2)       
    ax2 = subplot(2,1,2);
    scatter(datetime(time_grid(:),'convertfrom','datenum'),pres_grid(:),20,temp_grid(:),'filled'),
    ylabel('Depth [m]','interpreter','latex')
    set(ax2,'ydir','reverse','ticklabelinterpreter','latex'),
    colormap(ax2,cmocean('thermal')), caxis(ax2,[12 19])
    cb2 = colorbar; set(cb2,'ticklabelinterpreter','latex')
    ylabel(cb2,'Temp [C]','interpreter','latex')       
    yline(ax2,Pdye,'--r','linewidth',2)
    figname = [figdir,sprintf('Release%d_depth_time_scatter_vs_dye_temp_colormaps.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %
    %
    %
    % convert lat/lon to (cross,along) shore
    latitude = latitude_grid;
    longitude= longitude_grid;
    Np = size(pres_grid,1);
    latitude_grid = ones(Np,1)*latitude;
    longitude_grid = ones(Np,1)*longitude;
    % 
    [x,y] = lltoxy_AQR(latitude,longitude);
    y_grid = ones(Np,1)*y;
    x_grid = ones(Np,1)*x;
    %
    % break into two periods (early, late) to showcase shoaling
    t0 = datenum('08-Jul-2024 18:00:00');
    early = find(time< t0);
    late  = find(time> t0);    
    early_grid = find(time_grid< t0);
    late_grid  = find(time_grid> t0);
    %
    %
    % next we're going to do cross-shore and alongshore transects of early/late. Show those on a geoscatter plot...
    xi = [-150:10:150];
    pi = [0.5:0.2:13];
    yi = [-300:10:50];
    [xxi,yyi,ppi] = meshgrid(xi,yi,pi);
    %
    ybar_early    = nanmean(y(early));
    [lat_early, lon_early] = xytoll_AQR(xi,ybar_early*ones(1,length(xi)));
    xbar_late     = nanmean(x(late));
    [lat_late, lon_late] = xytoll_AQR(xbar_late*ones(1,length(yi)),yi);
    %
    figure,
    geoscatter(latitude,longitude, 60, log2(nanmean(dye_grid,1)),'filled'),
    hold on, geoplot(releaseLatitude,releaseLongitude,'ok',releaseLatitude,releaseLongitude,'xk','markersize',20,'linewidth',4,'markerfacecolor','m')
    geoplot(lat_early,lon_early,'--k',lat_late,lon_late,':k','linewidth',1.5)
    colormap(cm), caxis([0 8])
    geobasemap('satellite')
    geolimits([34.4667 34.4718],[-120.1324 -120.1246])
    ax = axes('position',[0.8 0.3 0.03 0.2]);
    imagesc(0,2.^conc,reshape(cm,256,1,3))
    ylabel('$\bar{D}$ [ppb]','interpreter','latex')
    set(ax,'ydir','normal','xtick',[],'ticklabelinterpreter','latex','xcolor','w','ycolor','w','fontsize',16)
    figname = [figdir,sprintf('Release%d_depth_avg_dye_vs_latitude_longitude.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %
    if ~isempty(early_grid)
        figure, scatter(latitude_grid(early_grid),pres_grid(early_grid),20,log2(dye_grid(early_grid)),'filled')
        %
        %
            tmp_x   = x_grid(:,early);    
    tmp_y   = y_grid(:,early);
    tmp_p   = pres_grid(:,early);
    tmp_dye = dye_grid(:,early); tmp_dye(tmp_dye==0)=1e-3;
    tmp_temp= temp_grid(:,early);
    valid  = ~isnan(tmp_temp) & ~isnan(tmp_dye);
    log2Di = griddata(tmp_x(valid),tmp_y(valid),tmp_p(valid),log2(tmp_dye(valid)), xxi, yyi, ppi);
    Ti     = griddata(tmp_x(valid),tmp_y(valid),tmp_p(valid),tmp_temp(valid), xxi, yyi, ppi);
    %
    %
    figure,
    ax1 = subplot(2,1,1);
    %scatter(x_grid(early_grid),pres_grid(early_grid),20,log2(dye_grid(early_grid)),'filled'),
    i1 = imagesc(xi,pi,squeeze(log2(nanmean(2.^log2Di,1)))');
    set(i1, 'AlphaData', ~isnan(squeeze(nanmean(log2Di,1))'))
    set(ax1,'ydir','reverse','ticklabelinterpreter','latex'), colormap(ax1,cm), caxis(ax1,[0 8]),
    ylabel('Depth [m]','interpreter','latex')
    cb1=colorbar;
    ticks = get(cb1,'ytick');
    set(cb1,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb1,'Dye [ppb]','interpreter','latex')
    yline(ax1,8,'--r','linewidth',2)
    xline(ax1,nanmean(x(late)),'--k','linewidth',2)
    titleStr = sprintf('$\\overline{\\Delta t} = %2.1f~\\mathrm{h}$~~~$\\bar{y}=%3.0f~\\mathrm{m}$',24*(mean(time(early))-releaseTime(1)),mean(y(early)));
    title(ax1,titleStr,'interpreter','latex')
    ax2 = subplot(2,1,2);
% $$$     scatter(x_grid(early_grid),pres_grid(early_grid),20,temp_grid(early_grid),'filled'),
    i2 = imagesc(xi,pi,squeeze(nanmean(Ti,1))');
    set(i2, 'AlphaData', ~isnan(squeeze(nanmean(Ti,1))'))    
    ylabel('Depth [m]','interpreter','latex')
    set(ax2,'ydir','reverse','ticklabelinterpreter','latex'),
    colormap(ax2,cmocean('thermal')), caxis(ax2,[12 19])
    cb2 = colorbar; set(cb2,'ticklabelinterpreter','latex')
    ylabel(cb2,'Temp [C]','interpreter','latex')       
    yline(ax2,8,'--r','linewidth',2)
    xline(ax2,nanmean(x(late)),'--k','linewidth',2)    
    figname = [figdir,sprintf('Release%d_dye_vs_alongshore_before_shoaling.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %
    %

    end
    if ~isempty(late_grid)
        figure, scatter(latitude_grid(late),pres_grid(late),20,log2(dye_grid(late)),'filled')
        %
            %
    tmp_x   = x_grid(:,late);    
    tmp_y   = y_grid(:,late);
    tmp_p   = pres_grid(:,late);
    tmp_dye = dye_grid(:,late); tmp_dye(tmp_dye==0)=1e-3;
    tmp_temp= temp_grid(:,late);
    tmp_time= (time_grid(:,late)-releaseTime(1))*24;
    valid  = ~isnan(tmp_temp) & ~isnan(tmp_dye);
    log2Di = griddata(tmp_x(valid),tmp_y(valid),tmp_p(valid),log2(tmp_dye(valid)), xxi, yyi, ppi);
    Ti     = griddata(tmp_x(valid),tmp_y(valid),tmp_p(valid),tmp_temp(valid), xxi, yyi, ppi);
    ti     = griddata(tmp_x(valid),tmp_y(valid),tmp_p(valid),tmp_time(valid), xxi, yyi, ppi);
% $$$     yi = [-200:10:30];
% $$$     pi = [0.5:0.2:13];
% $$$     [yyi,ppi] = meshgrid(yi,pi);
% $$$     tmp_y   = y_grid(:,late);
% $$$     tmp_p   = pres_grid(:,late);
% $$$     tmp_dye = dye_grid(:,late); tmp_dye(tmp_dye==0)=1e-3;
% $$$     tmp_temp= temp_grid(:,late);
% $$$     valid  = ~isnan(tmp_temp) & ~isnan(tmp_dye);
% $$$     log2Di = griddata(tmp_y(valid),tmp_p(valid),log2(tmp_dye(valid)), yyi, ppi);
% $$$     Ti     = griddata(tmp_y(valid),tmp_p(valid),tmp_temp(valid), yyi, ppi);
    %
    figure,
    ax1 = subplot(2,1,1);
    %scatter(x_grid(early_grid),pres_grid(early_grid),20,log2(dye_grid(early_grid)),'filled'),
    i1 = imagesc(yi,pi,squeeze(log2(nanmean(2.^log2Di,2)))');
    set(i1, 'AlphaData', ~isnan(squeeze(nanmean(log2Di,2))'))
    set(ax1,'ydir','reverse','ticklabelinterpreter','latex'), colormap(ax1,cm), caxis(ax1,[0 8]),
    ylabel('Depth [m]','interpreter','latex')
    cb1=colorbar;
    ticks = get(cb1,'ytick');
    set(cb1,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb1,'Dye [ppb]','interpreter','latex')
    yline(ax1,Pdye,'--r','linewidth',2)
    xline(ax1,nanmean(y(early)),'--k','linewidth',2)        
    titleStr = sprintf('$\\overline{\\Delta t} = %2.1f~\\mathrm{h}$~~~$\\bar{x}=%3.0f~\\mathrm{m}$',24*(mean(time(late))-releaseTime(1)),mean(x(late)));
    title(ax1,titleStr,'interpreter','latex')
    ax2 = subplot(2,1,2);
% $$$     scatter(x_grid(early_grid),pres_grid(early_grid),20,temp_grid(early_grid),'filled'),
    i2 = imagesc(yi,pi,squeeze(nanmean(Ti,2))');
    set(i2, 'AlphaData', ~isnan(squeeze(nanmean(Ti,2))'))    
    ylabel('Depth [m]','interpreter','latex')
    set(ax2,'ydir','reverse','ticklabelinterpreter','latex'),
    colormap(ax2,cmocean('thermal')), caxis(ax2,[12 19])
    cb2 = colorbar; set(cb2,'ticklabelinterpreter','latex')
    ylabel(cb2,'Temp [C]','interpreter','latex')       
    yline(ax2,Pdye,'--r','linewidth',2)
    xline(ax2,nanmean(y(early)),'--k','linewidth',2)            
    figname = [figdir,sprintf('Release%d_dye_vs_crossshore_after_shoaling.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %
    %
    %
    Di     = 2.^log2Di;
    Dsum_z = squeeze(nansum( Di , 3));
    mask_z = (Dsum_z>=1);
    Dbar_z = Dsum_z./squeeze(sum(~isnan(Di),3));
    zbar_D = squeeze(nansum( Di.*ppi , 3))./Dsum_z;
    zbar_D(~mask_z) = nan;
    z2bar_D = squeeze(nansum( Di.*(ppi-zbar_D).^2 , 3))./Dsum_z;
    z2bar_D(~mask_z) = nan;
    figure,
    i1 = imagesc(xi,yi,log2(Dbar_z));
    set(i1, 'AlphaData', ~isnan(Dbar_z))    
    set(gca,'ydir','normal','ticklabelinterpreter','latex')
    colormap(gca,cm), caxis(gca,[0 8]),
    xlabel('$x$~[m]','interpreter','latex')
    ylabel('$y$~[m]','interpreter','latex')
    title('Gridded $\overline{D}$','interpreter','latex')
    cb = colorbar;
    ylabel(cb,'$\bar{D}$ [ppb]','interpreter','latex')
    %
    figure,
    i3 = imagesc(xi,yi,squeeze(ti(:,:,1)));
    set(i3, 'AlphaData', ~isnan(Dbar_z))        
    set(gca,'ydir','normal','ticklabelinterpreter','latex')
    xlabel('$x$~[m]','interpreter','latex')
    ylabel('$y$~[m]','interpreter','latex')
    title('Gridded elapsed time $\Delta t$','interpreter','latex')
    cb = colorbar;
    ylabel(cb,'$\overline{\Delta t}$ [h]','interpreter','latex')
    figure,
    i2 = imagesc(xi,yi,squeeze(zbar_D));
    set(i2, 'AlphaData', ~isnan(Dbar_z))        
    set(gca,'ydir','normal','ticklabelinterpreter','latex')
    xlabel('$x$~[m]','interpreter','latex')
    ylabel('$y$~[m]','interpreter','latex')
    title('First moment $m_1 = m_0^{-1} \int z\,D\,\mathrm{d}z$','interpreter','latex')
    cb = colorbar;
    ylabel(cb,'$\bar{z}$ [m]','interpreter','latex')
    %
    figure,
    i2 = imagesc(xi,yi,sqrt(squeeze(z2bar_D)));
    set(i2, 'AlphaData', ~isnan(Dbar_z))        
    set(gca,'ydir','normal','ticklabelinterpreter','latex')
    xlabel('$x$~[m]','interpreter','latex')
    ylabel('$y$~[m]','interpreter','latex')
    title('Std-thickness $\sigma_z=m_2^{(1/2)}$','interpreter','latex')
    cb = colorbar; caxis([0 1])
    ylabel(cb,'$\sigma_z$ [m]','interpreter','latex')
    %
    %
    % map D(x,y,T)
    xi = [-250:10:100];
    Ti = [12:0.1:19];
    yi = [-250:10:50];
    [xxi,yyi,TTi] = meshgrid(xi,yi,Ti);
    tmp_x   = x_grid;    
    tmp_y   = y_grid;
    tmp_p   = pres_grid;    
    tmp_dye = dye_grid; tmp_dye(tmp_dye==0)=1e-3;
    tmp_temp= temp_grid;
    valid  = ~isnan(tmp_temp) & ~isnan(tmp_dye);
    log2Di = griddata(tmp_x(valid),tmp_y(valid),tmp_temp(valid),log2(tmp_dye(valid)), xxi, yyi, TTi);
    pi     = griddata(tmp_x(valid),tmp_y(valid),tmp_temp(valid),tmp_p(valid), xxi, yyi, TTi);
    %
    figure,
    ax1 = axes;
    %scatter(x_grid(early_grid),pres_grid(early_grid),20,log2(dye_grid(early_grid)),'filled'),
    i1 = imagesc(yi,Ti,squeeze(log2(nanmax(2.^log2Di,[],2)))');
    set(i1, 'AlphaData', ~isnan(squeeze(nanmax(log2Di,[],2))'))
    set(ax1,'ydir','normal','ticklabelinterpreter','latex'), colormap(ax1,cm), caxis(ax1,[0 8]),
    xlabel('$y$ [m]','interpreter','latex')
    ylabel('$T$ [$^\circ$ C]','interpreter','latex')
    cb1=colorbar;
    ticks = get(cb1,'ytick');
    set(cb1,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb1,'Dye [ppb]','interpreter','latex')
    titleStr = sprintf('$\\overline{\\Delta t} = %2.1f~\\mathrm{h}$~~~$\\bar{x}=%3.0f~\\mathrm{m}$',24*(mean(time)-releaseTime(1)),mean(x));
    title(ax1,titleStr,'interpreter','latex')
    figname = [figdir,sprintf('Release%d_dye_vs_crossshore_and_temperature.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %
    % find the center of the dye plume
    figure,
    ax1 = axes;
    %scatter(x_grid(early_grid),pres_grid(early_grid),20,log2(dye_grid(early_grid)),'filled'),
    i1 = imagesc(xi,Ti,squeeze(log2(nanmax(2.^log2Di,[],1)))');
    set(i1, 'AlphaData', ~isnan(squeeze(nanmax(log2Di,[],1))'))
    set(ax1,'ydir','normal','ticklabelinterpreter','latex'), colormap(ax1,cm), caxis(ax1,[0 8]),
    xlabel('$x$ [m]','interpreter','latex')
    ylabel('$T$ [$^\circ$ C]','interpreter','latex')
    cb1=colorbar;
    ticks = get(cb1,'ytick');
    set(cb1,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb1,'Dye [ppb]','interpreter','latex')
    titleStr = sprintf('$\\overline{\\Delta t} = %2.1f~\\mathrm{h}$~~~$\\bar{x}=%3.0f~\\mathrm{m}$',24*(mean(time)-releaseTime(1)),mean(x));
    title(ax1,titleStr,'interpreter','latex')
    figname = [figdir,sprintf('Release%d_dye_vs_alongshore_and_temperature.png',releaseNumber)];
    exportgraphics(gcf,figname)
    %
    end
    %
    %
    %
% $$$     dx = mean(gradient(x_grid(1,early)));
% $$$     dp = mean(gradient(pres_grid(:,1)));    
% $$$     TRI = delaunay(x_grid(:,early)/dx,pres_grid(:,early)/dp);
% $$$ figure, trisurf(TRI,x_grid(:,early),pres_grid(:,early),log2(dye_grid(:,early)),'edgecolor','none')
% $$$     dye_grid(isnan(dye_grid) | dye_grid==0)=1e-3;
    %
    %
    
