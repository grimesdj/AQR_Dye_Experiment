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
fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L0', 'ADCP');
fname = "ADCP_" + mooring + "_L0_10min.mat";
% if exist(fullfile(fpath, fname), 'file')
%     Vel(mooring_ID) = load(fullfile(fpath,fname));
% else
%     Vel(mooring_ID) = Vel(mooring_ID-1);
% end

% load vel EOF
savestr = mooring + "_EOF_depth_coords.mat";
path = fullfile(fpath, '..', '..', 'L1', 'ADCP');
if exist(fullfile(path, savestr), 'file')
    VEOF(mooring_ID) = load(fullfile(path, savestr));
end
end
%% Make Figures

% need to put velocities on depth grid!

figure("Position", [2250 30 750 1000])
H = arrayfun(@(x) x.dz(1), TEOF);
rowScale = H;
Ncols = 3;
ax = scaled_figure(rowScale,Ncols);


% M1 Velocity EOF mode 1
axes(ax(1, 1));
plot(abs(VEOF(1).EOFs(:, 1)), VEOF(1).dz, 'k-s', 'LineWidth', 2)
grid minor
axis ij
xticks([])


% M1 Temp variablilty
axes(ax(1, 3))
phi = TEOF(1).EC(:, 1);
A   = TEOF(1).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(1).dz, 'r-s', 'LineWidth', 2)
grid minor
axis ij
yticks([])
xticks([])

% M2 Velocity
axes(ax(2, 1))
plot(abs(VEOF(2).EOFs(:, 1)), VEOF(2).dz, 'k-s', 'LineWidth', 2)
grid minor
axis ij
xticks([])

% M2 Temp variablilty
axes(ax(2, 3))
phi = TEOF(2).EC(:, 1);
A   = TEOF(2).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(2).dz, 'r-s', 'LineWidth', 2)
grid minor
axis ij
yticks([])
xticks([])

% M3 Velocity
axes(ax(3, 1))
plot(abs(VEOF(3).EOFs(:, 1)), VEOF(3).dz, 'k-s', 'LineWidth', 2)
grid minor
axis ij
xticks([])

% M3 Temp variablilty
axes(ax(3, 3))
phi = TEOF(3).EC(:, 1);
A   = TEOF(3).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(3).dz, 'r-s', 'LineWidth', 2)
grid minor
axis ij
yticks([])
xticks([])

% M4 Vel
 % no data

% M4 Temp variablilty
axes(ax(4, 3))
phi = TEOF(4).EC(:, 1);
A   = TEOF(4).EOFs(:, 1);
S = phi * A';
sig = std(S, [], 1);
plot(sig, TEOF(4).dz, 'r-s', 'LineWidth', 2)
grid minor
axis ij
yticks([])

for rows = 1:length(H)
    for cols = [1 3]

        linkaxes(ax(rows, [1 3]), 'y')
        linkaxes(ax(1:3, cols), 'x')
       
    end
end

ylabel(ax(:, 1), '$h$ [m]', 'Interpreter','latex')
xlabel(ax(end, 1), '$U_{EOF1}$', 'Interpreter','latex')
xlabel(ax(end, 3), '$\sigma_{T}(^{\circ}C)$', 'Interpreter','latex')



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



return
% mean temperature profile
Y = Temp(1).Temp_grid';
Tbar = mean(Y, 1);








figure
plot(Tbar, TEOF(1).dz, 'k-s', 'LineWidth', 1.5)
grid on
axis ij
axis square
ylabel('Depth [m]')
xlabel('$^{\circ}C$', 'Interpreter','latex')
set(gca, 'FontSize', 18)



% u vs sig(T)