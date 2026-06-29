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
savestr = mooring + "_EOF.mat";
path = fullfile(fpath, '..', '..', 'L1', 'ADCP');
if exist(fullfile(path, savestr), 'file')
    VEOF(mooring_ID) = load(fullfile(path, savestr));
end
end
%% Make Figures

H = arrayfun(@(x) x.dz(1), TEOF);
rowH = H/sum(H);

lM  = 0.08;
rM  = 0.03;
tM  = 0.03;
bM  = 0.07;
gap = 0.03;

pW = 1 - lM - rM - 2*gap;
colW = pW/3;

pH = 1 - tM - bM - gap*(length(H)-1);

y = 1 - tM;

for i = 1:4

    h = pH * rowH(i);
    y = y - h;

    for j = 1:3

        x = lM + (j-1) * (colW+gap);

        ax(i, j) = axes('Position',[x y colW h]);

        % plot here
        Y = Temp(i).Temp_grid';
        Tbar = mean(Y, 1);
        plot(Tbar, TEOF(i).dz)
        grid minor
        axis ij
    end
    y = y-gap;
end



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