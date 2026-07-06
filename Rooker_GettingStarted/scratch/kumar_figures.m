%% Figures inspired by Kumar

clear all
close all

%% Load Data

moorings = {'M1', 'M2', 'M3', 'M4'};
for mooring_ID = 1:length(moorings)

% load temp
mooring = moorings{mooring_ID};
fprintf('loading %s data...\n', mooring)
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1');
savestr = mooring + "_10min_gridded.mat";
Temp(mooring_ID) = load(fullfile(fpath, savestr));

% load temp EOF
savestr = mooring + "_EOF.mat";
TEOF(mooring_ID) = load(fullfile(fpath, savestr));

% ADCP
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L1', 'ADCP');
fname = mooring + "_10min_gridded";
if exist(fullfile(fpath), 'dir')
    Vel.(mooring).North = load(fullfile(fpath,fname + "_North.mat"));
    Vel.(mooring).East = load(fullfile(fpath,fname + "_East.mat"));
else
    %
end

% load vel EOF
savestr = mooring + "_EOF_depth_coords.mat";
path = fullfile(fpath, '..', '..', 'L1', 'ADCP');
if exist(fullfile(path, savestr), 'file')
    VEOF(mooring_ID) = load(fullfile(path, savestr));
end
end
%% Make Figures

% need to put velocities on depth grid!

figure("Position", [2250 30 500 1000])
H = arrayfun(@(x) x.dz(1), TEOF);
rowScale = H;
Ncols = 2;
ax = scaled_figure(rowScale,Ncols);


% M1 Velocity EOF mode 1
axes(ax(1, 1));
plot(abs(VEOF(1).EOFs(:, 1)), VEOF(1).dz, 'k-s', 'LineWidth', 2)
grid minor
axis ij
xticklabels([])


% M1 Temp variablilty
axes(ax(1, 2))
phi = TEOF(1).EC(:, 1);
A   = TEOF(1).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(1).dz, 'r-s', 'LineWidth', 2)
grid(gca, 'minor')
axis ij
yticklabels([])
xticklabels([])

% M2 Velocity
axes(ax(2, 1))
plot(abs(VEOF(2).EOFs(:, 1)), VEOF(2).dz, 'k-s', 'LineWidth', 2)
grid(gca, 'minor')
axis ij
xticklabels([])

% M2 Temp variablilty
axes(ax(2, 2))
phi = TEOF(2).EC(:, 1);
A   = TEOF(2).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(2).dz, 'r-s', 'LineWidth', 2)
grid(gca, 'minor')
axis ij
yticklabels([])
xticklabels([])

% M3 Velocity
axes(ax(3, 1))
plot(abs(VEOF(3).EOFs(:, 1)), VEOF(3).dz, 'k-s', 'LineWidth', 2)
grid(gca, 'minor')
axis ij
xticklabels([])

% M3 Temp variablilty
axes(ax(3, 2))
phi = TEOF(3).EC(:, 1);
A   = TEOF(3).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(3).dz, 'r-s', 'LineWidth', 2)
grid(gca, 'minor')
axis ij
xticklabels([])
yticklabels([])

% M4 Vel
 % no data

% M4 Temp variablilty
axes(ax(4, 2))
phi = TEOF(4).EC(:, 1);
A   = TEOF(4).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(4).dz, 'r-s', 'LineWidth', 2)
grid(gca, 'minor')
axis ij
yticklabels([])

for rows = 1:length(H)
    linkaxes(ax(rows, [1 2]), 'y')
end
for cols = [1 2]
    %linkaxes(ax(1:4, cols), 'x')
end

ylabel(ax(:, 1), '$h$ [m]', 'Interpreter','latex')
xlabel(ax(end, 1), '$U_{EOF1}$', 'Interpreter','latex')
xlabel(ax(end, 2), '$\sigma_{T}(^{\circ}C)$', 'Interpreter','latex')



% rowH = H/sum(H);
% 
% lM  = 0.08;
% rM  = 0.03;
% tM  = 0.03;
% bM  = 0.07;
% gap = 0.03;
% 
% pW = 1 - lM - rM - 2*gap;
% colW = pW/3;
% 
% pH = 1 - tM - bM - gap*(length(H)-1);
% 
% y = 1 - tM;
% 
% for i = 1:4
% 
%     h = pH * rowH(i);
%     y = y - h;
% 
%     for j = 1:3
% 
%         x = lM + (j-1) * (colW+gap);
% 
%         ax(i, j) = axes('Position',[x y colW h]);
% 
%         % plot here
%         Y = Temp(i).Temp_grid';
%         Tbar = mean(Y, 1);
%         plot(Tbar, TEOF(i).dz)
%         grid minor
%         axis ij
%     end
%     y = y-gap;
% end



% mean temperature profile
figure
ax = scaled_figure(H(1:4), 1, 0.1, 0.1, 0.1, 0.05, 1);

% M1
axes(ax(1))
Y = Temp(1).Temp_grid';
Tbar = mean(Y, 1);
plot(Tbar, TEOF(1).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(1), 'on')
title('M1')

% M2
axes(ax(2))
Y = Temp(2).Temp_grid';
Tbar = mean(Y, 1);
plot(Tbar, TEOF(2).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(2), 'on')
title('M2')

% M3
axes(ax(3))
Y = Temp(3).Temp_grid';
Tbar = mean(Y, 1);
plot(Tbar, TEOF(3).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(3), 'on')
title('M3')

% M4
axes(ax(4))
Y = Temp(4).Temp_grid';
Tbar = mean(Y, 1);
plot(Tbar, TEOF(4).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(4), 'on')
title('M4')


% std temperature profile
figure
ax = scaled_figure(H(1:4), 1, 0.1, 0.1, 0.1, 0.05, 1);

% M1
axes(ax(1))
Y = Temp(1).Temp_grid';
Tbar = std(Y, [], 1);
plot(Tbar, TEOF(1).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(1), 'on')
title('M1')

% M2
axes(ax(2))
Y = Temp(2).Temp_grid';
Tbar = std(Y, [], 1);
plot(Tbar, TEOF(2).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(2), 'on')
title('M2')

% M3
axes(ax(3))
Y = Temp(3).Temp_grid';
Tbar = std(Y, [], 1);
plot(Tbar, TEOF(3).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(3), 'on')
title('M3')

% M4
axes(ax(4))
Y = Temp(4).Temp_grid';
Tbar = std(Y, [], 1);
plot(Tbar, TEOF(4).dz, 'k-s', 'LineWidth', 1.5)
axis ij
grid(ax(4), 'on')
title('M4')





% mean velocity profile
figure
ax = scaled_figure(H(1:3), 1, 0.1, 0.1, 0.1, 0.05, 1);

% M1
axes(ax(1))
N = Vel.M1.North.Vel_grid;
N = mean(N, 2);
E = Vel.M1.East.Vel_grid;
E = mean(E, 2);
plot(N, VEOF(1).dz, 'k-s', 'LineWidth', 1.5)
hold on
plot(E, VEOF(1).dz, 'r-s', 'LineWidth', 1.5)
axis ij
grid(ax(1), 'minor')
xline(0, 'k--', 'LineWidth', 1)
xlim([-0.05 0.05])
title('M1')

% M2
axes(ax(2))
N = Vel.M2.North.Vel_grid;
N = mean(N, 2);
E = Vel.M2.East.Vel_grid;
E = mean(E, 2);
plot(N, VEOF(2).dz, 'k-s', 'LineWidth', 1.5)
hold on
plot(E, VEOF(2).dz, 'r-s', 'LineWidth', 1.5)
axis ij
grid(ax(2), 'minor')
xline(0, 'k--', 'LineWidth', 1)
xlim([-0.05 0.05])
title('M2')

% M3
axes(ax(3))
N = Vel.M3.North.Vel_grid;
N = mean(N, 2);
E = Vel.M3.East.Vel_grid;
E = mean(E, 2);
plot(N, VEOF(3).dz, 'k-s', 'LineWidth', 1.5)
hold on
plot(E, VEOF(3).dz, 'r-s', 'LineWidth', 1.5)
axis ij
grid(ax(3), 'minor')
xline(0, 'k--', 'LineWidth', 1)
xlim([-0.05 0.05])
title('M3')

% std velocity profile
figure
ax = scaled_figure(H(1:3), 1, 0.1, 0.1, 0.1, 0.05, 1);

% M1
axes(ax(1))
N = Vel.M1.North.Vel_grid;
N = std(N, [], 2);
E = Vel.M1.East.Vel_grid;
E = std(E, [], 2);
plot(N, VEOF(1).dz, 'k-s', 'LineWidth', 1.5)
hold on
plot(E, VEOF(1).dz, 'r-s', 'LineWidth', 1.5)
axis ij
grid(ax(1), 'minor')
xline(0, 'k--', 'LineWidth', 1)
title('M1')

% M2
axes(ax(2))
N = Vel.M2.North.Vel_grid;
N = std(N, [], 2);
E = Vel.M2.East.Vel_grid;
E = std(E, [], 2);
plot(N, VEOF(2).dz, 'k-s', 'LineWidth', 1.5)
hold on
plot(E, VEOF(2).dz, 'r-s', 'LineWidth', 1.5)
axis ij
grid(ax(2), 'minor')
xline(0, 'k--', 'LineWidth', 1)
title('M2')

% M3
axes(ax(3))
N = Vel.M3.North.Vel_grid;
N = std(N, [], 2);
E = Vel.M3.East.Vel_grid;
E = std(E, [], 2);
plot(N, VEOF(3).dz, 'k-s', 'LineWidth', 1.5)
hold on
plot(E, VEOF(3).dz, 'r-s', 'LineWidth', 1.5)
axis ij
grid(ax(3), 'minor')
xline(0, 'k--', 'LineWidth', 1)
title('M3')

return



figure
plot(Tbar, TEOF(1).dz, 'k-s', 'LineWidth', 1.5)
grid on
axis ij
axis square
ylabel('Depth [m]')
xlabel('$^{\circ}C$', 'Interpreter','latex')
set(gca, 'FontSize', 18)



% u vs sig(T)


