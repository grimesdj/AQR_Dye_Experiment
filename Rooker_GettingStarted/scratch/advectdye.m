%% Advect DYE?

clear all
close all

%% Load
load('../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat')
AQD = load('../../../../Kelp_data/data/Release2/L0/ADCP/AQD_ADCP.mat');
%% Assumtions:
% no diffusion
% horizontal homogenaety or however you spell that
% dye density == water density
% uniform stratification


%% North

for dirs = 1:2
    cards = {'North', 'East'};
    card = cards{dirs};

    %% interpolate Velocitiy
    [t_grid, z_grid] = meshgrid(AQD.Time, AQD.dbins);
    msk = ~isnan(AQD.("Velocity_" + card));
    u = AQD.("Velocity_" + card)(msk);
    t = t_grid(msk);
    z = z_grid(msk);
    F = scatteredInterpolant(t(:), z(:), u(:), "linear", "nearest");
    
    Nt = length(R2.Time);
    T = R2.Temperature(:, 5);
    Tbar = mean(R2.Temperature);
    dTdz = mean(diff(Tbar)./diff(R2.mab));
    dt = mode(diff(R2.Time))*86400;
    A = zeros(1, Nt);
    z = ones(1, Nt);
    c = 1;
% release new dye every 10 mins?
for j = 2:round(600/dt):3600/dt
    fprintf('Dye release at %s\n', datestr(R2.Time(j)))
    A1 = zeros(1, Nt);
    z1 = ones(1, Nt);
    for i = j:Nt
            z1(i) = z1(1) + (T(i) - T(j))*1./dTdz; % vertical position
            u_inst = F(R2.Time(i), z1(i));
            A_inst = u_inst * dt;
            A1(i) = A1(i-1) + A_inst;
    end

    A(c, :) = A1;
    z(c, :) = z1;
    c = c+1;
end

    figure
    plot(R2.Time, z)
    title(sprintf('z vs time %s', card))
    drawnow
    
    figure
    plot(R2.Time, A)
    title(sprintf('A vs time %s', card))
    drawnow

    figure
    plot(A, z, 'r.')
    title(sprintf('A vs z %s', card))
    L.(card) = A;

end
    
   


%% Both
cmap = cmocean('deep');
n = length(L.North);
idx = round(linspace(1, size(cmap,1), n));
colors = cmap(idx,:);
colors = colors;

figure
scatter(L.East,L.North, 10, colors, "filled")
axis equal
drawnow

%% mappp
c = projcrs(32611);
lat = dms2degrees([34 28 8.68]);
lon = dms2degrees([-120 7 40.28]);
[origin_x, origin_y] = projfwd(c, lat, lon);
X = origin_x + L.East;
Y = origin_y + L.North;
[lats, lons] = projinv(c, X, Y);
figure
geoscatter(lats(:), lons(:),  [], z(:), "filled")
geobasemap('satellite')
set(gca, 'FontSize', 18)
clim([min(z, [], 'all'), max(z, [], 'all')])
colormap(flipud(cmap))

ax = gca;
c1 = axes('Position', [ax.Position(1)+ax.Position(3)*0.7, ax.Position(2)+ax.Position(4)*0.2, 0.01, ax.Position(4)/4.25]);
clims = clim(ax);
imagesc(c1, [1 2], clims, (1:size(cmap, 1))')
set(gca, 'YDir', 'normal')
ylabel('Dye MAB')
colormap(flipud(cmap))
c1.YColor = 'white';
set(c1, 'FontSize', 18)
ylim(c1,clims)
set(c1,'XTick',[])
drawnow

[latlim, lonlim] = geolimits(ax);
latPad = 0.10 * diff(latlim);
lonPad = 0.10 * diff(lonlim);
geolimits(ax, ...
    [latlim(1)-latPad, latlim(2)+latPad], ...
    [lonlim(1)-lonPad, lonlim(2)+lonPad]);