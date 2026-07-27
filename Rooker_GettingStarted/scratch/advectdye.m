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
x = AQD.Velocity_North;
wpass = 1/600;
fs = 1/round(mean(diff(AQD.Time))*86400);
 [y, fsd, idx] = hamming_filter(x, wpass, fs);
%% interpolate Velocitiy
[t_grid, z_grid] = meshgrid(AQD.Time, AQD.dbins);
msk = ~isnan(y);
u = y(msk);
t = t_grid(msk);
z = z_grid(msk);
F = scatteredInterpolant(t(:), z(:), u(:), "linear", "nearest");

Time = ADCP1.Time(1:300:end);
Nt = length(Time);

T = ADCP1.Temperature(:, 5);
Tbar = mean(ADCP1.Temperature);
dTdz = mean(diff(Tbar)./diff(ADCP1.mab));
T0 = T1dye;
dt = mean(diff(Time))*86400; % HAVE TO SMOOTH TEMPERATURE TOO
c = 1;
% [Tq, Zq] = meshgrid(ADCP1.Time, AQD.dbins);
% U = F(Tq, Zq);
figure
hold on
for j = 2:13
    fprintf('j == %d\n', j)
    z = 1; % mab
    A = zeros(1, j); % meters from origin
    for i = j:Nt
        z(i) = z(1) + (T(i) - T(j-1))*1./dTdz; % vertical position
        u_inst = F(ADCP1.Time(i), z(i));
        A_inst = u_inst * dt;
        A(i) = A(i-1) + A_inst;
    end

    La(:, c) = A;
    Lz(:, c) = z;
    c = c+1;
    plot(Time, A)
    drawnow
end


figure
plot(z)
drawnow

figure
plot(A)
drawnow

figure
plot(La(:), Lz(:), 'r.')
drawnow
keyboard
%% East

msk = ~isnan(AQD.Velocity_North);
u = AQD.Velocity_East(msk);
t = t_grid(msk);
z = z_grid(msk);
F = scatteredInterpolant(t(:), z(:), u(:), "linear", "nearest");
Nt = length(ADCP1.Time);
z = 1; % mab
B = 0; % meters from origin
T = ADCP1.Temperature(:, 5);
Tbar = mean(ADCP1.Temperature);
dTdz = mean(diff(Tbar)./diff(ADCP1.mab));
T0 = T1dye;
dt = mean(diff(ADCP1.Time))*86400*10;

for i = 2:Nt
    z(i) = z(i-1) + (T(i) - T(i-1))*1./dTdz; % vertical position
    u_inst = F(ADCP1.Time(i), z(i));
    B_inst = u_inst * dt;
    B(i) = B(i-1) + B_inst;
end

figure
plot(z)
drawnow

figure
plot(B)
drawnow

figure
plot(B(:), z(:), 'r.')
drawnow



%% Both
cmap = cmocean('deep');
n = length(A);
idx = round(linspace(1, size(cmap,1), n));
colors = cmap(idx,:);

figure
scatter(A, B, 10, colors, "filled")
axis equal
drawnow

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

drawnow