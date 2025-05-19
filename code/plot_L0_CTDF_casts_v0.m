clear all
% full path
releaseNumber  = 3;
ctdSerialNumber=234870;
rootDir  = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/'
dataDir  = sprintf('%s/data/Release%d/',rootDir,releaseNumber);
L0Dir    = [dataDir,filesep,'L0',];
fileRoot = sprintf('%d*.mat',ctdSerialNumber);
files    = dir([L0Dir,filesep,fileRoot]);
for jj  = 1:length(files)
    SN  = split(files(jj).name,'_');
    fin= [L0Dir,filesep,SN{1},'_L0.mat'];
    load(fin);
    %
% $$$     figure, scatter(dye,pres, 10,temp,'filled')
% $$$     figure, plot(temp,dye_raw,'.')
    %
    [dye_sort,srt] = sort(dye_raw,'ascend');
    temp_sort = temp(srt);
    pres_sort = pres(srt);
    [cm,conc] = make_DYE_colormap([0 8],1);
    %    
    figure,
    scatter(temp_sort,pres_sort,15,log2(dye_sort),'filled'), set(gca,'ydir','reverse','ticklabelinterpreter','latex'), colormap(cm), caxis([0 8])
    xline(15,'--r','linewidth',2)
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
    figname = [figdir,sprintf('%d_depth_temp_scatter_dye_colormap.png',ctdSerialNumber)];
    exportgraphics(gcf,figname)
    %    
    %
    figure,
    ax1 = subplot(2,1,1);
    scatter(datetime(time,'convertfrom','datenum'),pres,20,log2(dye_raw),'filled'),
    set(ax1,'ydir','reverse','ticklabelinterpreter','latex'), colormap(ax1,cm), caxis(ax1,[0 8]),
    ylabel('Depth [m]','interpreter','latex')
    cb1=colorbar;
    ticks = get(cb1,'ytick');
    set(cb1,'yticklabel', cellfun(@num2str,num2cell(2.^ticks)','uniformoutput',0),'ticklabelinterpreter','latex')
    ylabel(cb1,'Dye [ppb]','interpreter','latex')
    yline(ax1,8,'--r','linewidth',2)       
    ax2 = subplot(2,1,2);
    scatter(datetime(time,'convertfrom','datenum'),pres,20,temp,'filled'),
    ylabel('Depth [m]','interpreter','latex')
    set(ax2,'ydir','reverse','ticklabelinterpreter','latex'),
    colormap(ax2,cmocean('thermal')), caxis(ax2,[12 19])
    cb2 = colorbar; set(cb2,'ticklabelinterpreter','latex')
    ylabel(cb2,'Temp [C]','interpreter','latex')       
    yline(ax2,8,'--r','linewidth',2)
    figname = [figdir,sprintf('%d_depth_time_scatter_vs_dye_temp_colormaps.png',ctdSerialNumber)];
    exportgraphics(gcf,figname)
end