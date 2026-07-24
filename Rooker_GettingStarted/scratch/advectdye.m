%% Advect DYE?

clear all
close all

%% Load
load('../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat')
AQD = load('../../../../Kelp_data/data/Release1/L0/ADCP/AQD_ADCP.mat');
%% Assumtions:
% no diffusion
% horizontal homogenaety or however you spell that
% dye density == water density
% uniform stratification


%% North

%% interpolate Velocitiy
[t_grid, z_grid] = meshgrid(AQD.Time, AQD.dbins);
msk = ~isnan(AQD.Velocity_North);
u = AQD.Velocity_North(msk);
t = t_grid(msk);
z = z_grid(msk);
F = scatteredInterpolant(t(:), z(:), u(:), "linear", "nearest");

z = 1; % mab
A = 0; % meters from origin
Nt = length(ADCP1.Time);
T = ADCP1.Temperature(:, 5);
Tbar = mean(ADCP1.Temperature);
dTdz = mean(diff(Tbar)./diff(ADCP1.mab));
T0 = T1dye;
dt = mean(diff(ADCP1.Time))*86400;

for i = 1:Nt
    z(i+1) = z(1) + (T(i) - T0)*1./dTdz; % vertical position
    u_inst = (F(ADCP1.Time(i), z(i+1)));
    A_inst = u_inst * dt;
    A(i+1) = A(i) + A_inst;
end

figure
plot(z)

figure
plot(A)

figure
plot(A, z, 'r.')

%% East

msk = ~isnan(AQD.Velocity_East);
u = AQD.Velocity_East(msk);
t = t_grid(msk);
z = z_grid(msk);
F = scatteredInterpolant(t(:), z(:), u(:), "linear", "nearest");

z = 1; % mab
B = 0; % meters from origin
Nt = length(ADCP1.Time);
T = ADCP1.Temperature(:, 5);
Tbar = mean(ADCP1.Temperature);
dTdz = mean(diff(Tbar)./diff(ADCP1.mab));

for i = 1:Nt
    z(i+1) = z(1) + (T(i) - T0)*1./dTdz; % vertical position
    u_inst = (F(ADCP1.Time(i), z(i+1)));
    B_inst = u_inst * dt;
    B(i+1) = B(i) + B_inst;
end

figure
plot(z)

figure
plot(B)

figure
plot(B, z, 'r.')


%% Both
cmap = cmocean('deep');
n = length(A);
idx = round(linspace(1, size(cmap,1), n));
colors = cmap(idx,:);

figure
scatter(B, A, 10, colors, "filled")
axis equal


%% mappp
c = projcrs(32611);
lat = dms2degrees([34 28 8.68]);
lon = dms2degrees([-120 7 40.28]);
[origin_x, origin_y] = projfwd(c, lat, lon);
X = origin_x + B;
Y = origin_y + A;
[lats, lons] = projinv(c, X, Y);
figure
geoscatter(lats, lons,  [], colors, "filled")
geobasemap('satellite')
set(gca, 'FontSize', 18)

ax = gca;
c1 = axes('Position', [ax.Position(1)+ax.Position(3)*0.8, ax.Position(2)+ax.Position(4)*0.75, 0.01, ax.Position(4)/4.25]);
clim = ylim(c1);
axe = [1;1] * linspace(clim(1), clim(2), length(cmap)) .* max(z)./256;
imagesc(c1, axe')
colormap(cmap)
c1.YColor = 'white';
set(c1, 'FontSize', 18)