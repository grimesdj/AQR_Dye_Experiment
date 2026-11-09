%% Wake Scatter
% inspired by Rosman et al 07

clear all
close all

%% load data

fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M1', 'L1', 'ADCP');
fname = 'M1_10min_gridded_PCA.mat';
M1 = load(fullfile(fpath, fname));

fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M2', 'L1', 'ADCP');
fname = 'M2_10min_gridded_PCA.mat';
M2 = load(fullfile(fpath, fname));

fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', 'M3', 'L1', 'ADCP');
fname = 'M3_10min_gridded_PCA.mat';
M3 = load(fullfile(fpath, fname));

x1s = [mean(M1.U_grid(:, :), 1); M1.U_grid(end-1, :); M1.U_grid(2, :)];
x2s = [mean(M2.U_grid(:, :), 1); M2.U_grid(end-1, :); M2.U_grid(2, :)];
x3s = [mean(M3.U_grid(:, 1:end-1), 1); M3.U_grid(end-1, 1:end-1); M3.U_grid(2, 1:end-1)];

y1s = [mean(M1.V_grid(:, :), 1); M1.V_grid(end-1, :); M1.V_grid(2, :)];
y2s = [mean(M2.V_grid(:, :), 1); M2.V_grid(end-1, :); M2.V_grid(2, :)];
y3s = [mean(M3.V_grid(:, 1:end-1), 1); M3.V_grid(end-1, 1:end-1); M3.V_grid(2, 1:end-1)];


name = {'Depth-Averaged', '~1 MAB', 'Near-surface'};

for i = 1:3
   
    x1 = x1s(i,:);
    x2 = x2s(i,:);
    x3 = x3s(i,:);

    y1 = y1s(i,:);
    y2 = y2s(i,:);
    y3 = y3s(i,:);

    %% plots


    %fig(i) = figure('Position', [0 0 500 1000]);
    
    
    %% U M3 M1
    figure('Position', [0 0 1000 500])%subplot(3, 2, 1)
    subplot(1, 2, 1)
    plot(x3, x1, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    xlim([-0.25 0.25])
    ylim([-0.25 0.25])
    xlabel('$U_{M3}$ [m/s]', 'Interpreter','latex')
    ylabel('$U_{M1}$ [m/s]', 'Interpreter','latex')
    set(gca, 'fontsize', 15)
    grid on
    axis square
    
    % X = [x3' x1'];
    % coeff = pca(X);
    % v = coeff(:, 1);
    % m = v(2)/v(1);
    %t = linspace(-0.2, 0.2);
    % y = m .* t;
    
    west = x1 < 0 & x3 < 0;
    WX = [x3(west)', x1(west)'];
    
    east = x1 > 0 & x3 > 0;
    EX = [x3(east)', x1(east)'];

    Wcoeff = pca(WX);
    Ecoeff = pca(EX);

    Wv = Wcoeff(:, 1);
    Ev = Ecoeff(:, 1);

    Wm = Wv(2)/Wv(1);
    Em = Ev(2)/Ev(1);

    Wt = linspace(-0.2, 0);
    Et = linspace(0, 0.2);

    Wy = Wm .* Wt;
    Ey = Em .* Et;

    hold on
    plot(Wt, Wy, 'r--', 'LineWidth', 2)
    plot(Et, Ey, 'r--', 'LineWidth', 2)

    text(0, -0.2, sprintf('m1 = %.2f\nm2 = %.2f', Wm, Em), ...
        'FontSize', 20, ...
        'Color', 'r', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.5)
    



    %% V M3 M1
    %subplot(3, 2, 2)
    subplot(1, 2, 2)
    plot(y3, y1, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    xlim([-0.25 0.25])
    ylim([-0.25 0.25])
    xlabel('$V_{M3}$ [m/s]', 'Interpreter','latex')
    ylabel('$V_{M1}$ [m/s]', 'Interpreter','latex')
    set(gca, 'fontsize', 15)
    grid on
    axis square
    % 
    % X = [y3' y1'];
    % coeff = pca(X);
    % v = coeff(:, 1);
    % m = v(2)/v(1);
    % t = linspace(-0.2, 0.2);
    % y = m .* t;
    % 
    % hold on
    % plot(t, y, 'r--', 'LineWidth', 1.5)


    west = y1 < 0 & y3 < 0;
    WX = [y3(west)', y1(west)'];
    
    east = y1 > 0 & y3 > 0;
    EX = [y3(east)', y1(east)'];

    Wcoeff = pca(WX);
    Ecoeff = pca(EX);

    Wv = Wcoeff(:, 1);
    Ev = Ecoeff(:, 1);

    Wm = Wv(2)/Wv(1);
    Em = Ev(2)/Ev(1);

    Wt = linspace(-0.2, 0);
    Et = linspace(0, 0.2);

    Wy = Wm .* Wt;
    Ey = Em .* Et;

    hold on
    plot(Wt, Wy, 'r--', 'LineWidth', 2)
    plot(Et, Ey, 'r--', 'LineWidth', 2)


    text(0, -0.2, sprintf('m1 = %.2f\nm2 = %.2f', Wm, Em), ...
        'FontSize', 20, ...
        'Color', 'r', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.5)

    sgtitle(name{i})
    fpath = fullfile('..', 'figures');
    exportgraphics(gcf, fullfile(fpath, sprintf('%s_scatter.png', name{i})), 'Resolution',600)
end




return 
sgtitle(fig(1), 'Depth Averaged')
sgtitle(fig(2), '~1 MAB')
sgtitle(fig(3), 'Near Surface')


fpath = fullfile('..', 'figures');
exportgraphics(fig(1), fullfile(fpath, 'Depth_avg_scatter.png'), 'Resolution',600)
exportgraphics(fig(2), fullfile(fpath, 'Near_bottom_scatter.png'), 'Resolution',600)
exportgraphics(fig(3), fullfile(fpath, 'Near_Surface_scatter.png'), 'Resolution',600)
