% Clean 10-minute average velocity processing and mapping script

clear all; close all;

%% File Setup and Color Definitions
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/*.mat');
colors = {[0, 0, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0]};
labels = cell(1, length(files));
h1 = gobjects(1, length(files));
h2 = gobjects(1, length(files));

%% Time Series Plot Setup
figure;
ax1 = subplot(2, 1, 1); hold(ax1, 'on'); ylabel(ax1, 'East Velocity');
ax2 = subplot(2, 1, 2); hold(ax2, 'on'); ylabel(ax2, 'North Velocity'); xlabel(ax2, 'Time');
datetick(ax1, 'x', 'keeplimits'); datetick(ax2, 'x', 'keeplimits');
sgtitle('10-min Avg Velocities at 1.2 MAB')

%% Quiver Plot Setup
figure;
for i = 1:length(files)
    data = load(fullfile(files(i).folder, files(i).name));

    % Determine appropriate bin
    if size(data.Velocity_East, 2) == 1
        bin = 1;
    else
        bin = round((1.237 - data.Config.blank) / data.Config.binSize);
    end

    % Extract and clean data
    u = data.Velocity_East(:, bin);
    v = data.Velocity_North(:, bin);
    t = data.Time;
    dt = double(data.Config.dt);

    % Lowpass filter
    Nf = 600 / dt;
    flt = hamming(Nf); flt = flt / sum(flt);
    nanFlag = ~isnan(u);
    flag_flt = conv(nanFlag, flt, 'same');
    u(~nanFlag) = 0; v(~nanFlag) = 0;
    u_filt = conv(u, flt, 'same') ./ flag_flt;
    v_filt = conv(v, flt, 'same') ./ flag_flt;

    % Plot time series
    figure(1);
    h1(i) = plot(ax1, t, u_filt);
    h2(i) = plot(ax2, t, v_filt);
    labels{i} = data.Config.SN;

    % Downsample to 10-minute intervals
    [~, unique_idx] = unique(minutes(t - t(1)));
    t = t(unique_idx);
    u_filt = u_filt(unique_idx);
    v_filt = v_filt(unique_idx);
    tq = t(1):minutes(10):t(end);
    u_ds = interp1(t, u_filt, tq, 'linear');
    v_ds = interp1(t, v_filt, tq, 'linear');

    % Plot quiver per instrument
    figure(2)
    subplot(length(files), 1, i)
    quiver(datenum(tq), zeros(size(tq)), u_ds, v_ds, 'AutoScale', 'off')
    ylabel('Velocity'); title(data.Config.SN)
    datetick('x', 'keeplimits')

    % Save Global Vars
    U_all{i} = u_ds;
    V_all{i} = v_ds;
    T_all{i} = tq;
    
end
xlabel('Time')
sgtitle('Current Vectors at 1.2 MAB')

% Add legends
figure(1)
legend(ax1, h1, labels, 'Location', 'best')
legend(ax2, h2, labels, 'Location', 'best')
datetick(ax1); datetick(ax2)

% Set colors
for i = 1:length(files)
    h1(i).Color = colors{i};
    h2(i).Color = colors{i};
end
yline(ax1, 0); yline(ax2, 0);

%% Approximate Map Coordinates
lats = [34.4690789, 34.469615, 34.4690706, 34.4690789];
longs = [-120.1278549, -120.130375, -120.1278088, -120.1278549];

% Setup map figure
figure;
geobasemap satellite; hold on
for i = 1:length(lats)
    geoplot(lats(i), longs(i), 'o', 'Color', colors{i})
end
close(3)



[~, sizeLimit] = size(U_all{2});

% Compare inside instruments to outside instruments
AQDvsM1 = [U_all{1}(1 , 1:sizeLimit)-U_all{2}; V_all{1}(1, 1:sizeLimit)-V_all{2}];
M2vsM1  = [U_all{3}-U_all{2}; V_all{3}-V_all{2}];
VECvsM1 = [U_all{4}(1 , 1:sizeLimit)-U_all{2}; V_all{4}(1, 1:sizeLimit)-V_all{2}];

%% Plot difference scatter in-out
figure
scatter(AQDvsM1(1,:), AQDvsM1(2,:), 'b.')
hold on
scatter(M2vsM1(1,:),M2vsM1(2,:) , 'm.')
scatter(VECvsM1(1,:),VECvsM1(2,:) , 'g.')

% Generate and plot trendlines
AQDtrend = polyfit(AQDvsM1(1,:), AQDvsM1(2,:), 1);
yfit = polyval(AQDtrend, AQDvsM1(1,:));
plot(AQDvsM1(1,:), yfit, 'b-', 'LineWidth', 2)

M2trend = polyfit(M2vsM1(1,:), M2vsM1(2,:), 1);
yfit = polyval(M2trend, M2vsM1(1,:));
plot(M2vsM1(1,:), yfit, 'm-', 'LineWidth', 2)

VECtrend = polyfit(VECvsM1(1,:), VECvsM1(2,:), 1);
yfit = polyval(VECtrend, VECvsM1(1,:));
plot(VECvsM1(1,:), yfit, 'g-', 'LineWidth', 2)

% Formatting and such
title("Inside Instruments vs M1 ADCP ( X - M1 )", "FontSize", 18)
legend(labels{[1 3 4]})
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)

%% Compare AQD - VEC
AQDvsVEC = [U_all{1}(1 , 1:sizeLimit)-U_all{4}(1, 1:sizeLimit); V_all{1}(1, 1:sizeLimit)-V_all{4}(1, 1:sizeLimit)];
figure, scatter(AQDvsVEC(1,:), AQDvsVEC(2,:), '.')

% mean and covariance
mu = [mean(AQDvsVEC(1,:)), mean(AQDvsVEC(2,:))];
Sigma = cov(AQDvsVEC(1,:), AQDvsVEC(2,:));

% eigen decomposition
[eigvec,eigval] = eig(Sigma);

% parametric ellipse
theta = linspace(0,2*pi,200);
ellipse = [cos(theta); sin(theta)]';

% scale by sqrt of eigenvalues (std dev)
ellipse = ellipse * sqrt(eigval) * eigvec';

% shift to mean
ellipse = ellipse + mu;
hold on
plot(ellipse(:,1), ellipse(:,2), 'r','LineWidth',2)

% Formatting
title([labels{1} ' vs ' labels{4}], "FontSize", 18)
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)
legend('AQD vs  VEC', '1 \sigma')

%% Compare AQD - M2
AQDvsM2 = [U_all{1}(1 , 1:sizeLimit)-U_all{3}(1, 1:sizeLimit); V_all{1}(1, 1:sizeLimit)-V_all{3}(1, 1:sizeLimit)];
figure, scatter(AQDvsM2(1,:), AQDvsM2(2,:), '.')

% mean and covariance
mu = [mean(AQDvsM2(1,:)), mean(AQDvsM2(2,:))];
Sigma = cov(AQDvsM2(1,:), AQDvsM2(2,:));

% eigen decomposition
[eigvec,eigval] = eig(Sigma);

% parametric ellipse
theta = linspace(0,2*pi,200);
ellipse = [cos(theta); sin(theta)]';

% scale by sqrt of eigenvalues (std dev)
ellipse = ellipse * sqrt(eigval) * eigvec';

% shift to mean
ellipse = ellipse + mu;
hold on
plot(ellipse(:,1), ellipse(:,2), 'r','LineWidth',2)

% Formatting
title([labels{1} ' vs ' labels{3}], "FontSize", 18)
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)
legend('AQD vs M2', '1 \sigma')

%% Compare VEC - M2
VECvsM2 = [U_all{4}(1 , 1:sizeLimit)-U_all{3}(1, 1:sizeLimit); V_all{4}(1, 1:sizeLimit)-V_all{3}(1, 1:sizeLimit)];
figure, scatter(VECvsM2(1,:), VECvsM2(2,:), '.')

% mean and covariance
mu = [mean(VECvsM2(1,:)), mean(VECvsM2(2,:))];
Sigma = cov(VECvsM2(1,:), VECvsM2(2,:));

% eigen decomposition
[eigvec,eigval] = eig(Sigma);

% parametric ellipse
theta = linspace(0,2*pi,200);
ellipse = [cos(theta); sin(theta)]';

% scale by sqrt of eigenvalues (std dev)
ellipse = ellipse * sqrt(eigval) * eigvec';

% shift to mean
ellipse = ellipse + mu;
hold on
plot(ellipse(:,1), ellipse(:,2), 'r','LineWidth',2)

% Formatting
title([labels{4} ' vs ' labels{3}], "FontSize", 18)
ylabel('Latitude Current Difference (m/s)', 'FontSize', 14)
xlabel('Longitude Current Difference (m/s)', 'FontSize', 14)
legend('VEC vs M2', '1 \sigma')

% %% Make MP4 of Vectors on Map with Midpoint Diff Arrow
% v = VideoWriter('../../../../Kelp_data/Summer2025/Rooker/figures/Release2/animations/Release2_Currents.mp4', 'MPEG-4');
% v.FrameRate = 5;
% open(v)
% 
% % Create visible figure
% fig = figure('Visible','on');
% gx = geoaxes(fig);
% grid(gx,'off');
% hold(gx,'on');
% geobasemap(gx,'satellite');
% 
% % --- Create dummy arrows for persistent legend ---
% hArrowsLegend = gobjects(length(colors)+1, 1); % +1 for diff arrow
% for i = 1:length(colors)
%     hArrowsLegend(i) = geoplot(gx, [0 0], [0 0], 'Color', colors{i}, 'LineWidth', 2);
% end
% % Dummy for diff arrow (black)
% lon_diffdot = [];
% lat_diffdot = [];
% 
% hArrowsLegend(end) = geoplot(gx, [0 0], [0 0], 'Color', [0 0 0], 'LineWidth', 2);
% 
% legend(hArrowsLegend, [labels, {'Diff 2-3'}], 'Location', 'northeast');
% 
% % --- Animate arrows and base circles ---
% for idx = 1:length(tq)
%     % Delete previous arrows/markers
%     if exist('hArrows','var'); delete(hArrows); end
%     if exist('hMarkers','var'); delete(hMarkers); end
% 
%     hold(gx,'on');
%     for i = 1:length(lats)
%         if idx > length(U_all{i}); continue; end
% 
%         lat_base = lats(i);
%         lon_base = longs(i);
%         scale = 0.03;
%         lat_tip(i) = lat_base + V_all{i}(idx) * scale;
%         lon_tip(i) = lon_base + U_all{i}(idx) * scale;
% 
%         % Plot arrow
%         hArrows(i) = geoplot(gx, [lat_base, lat_tip(i)], [lon_base, lon_tip(i)], ...
%                              'Color', colors{i}, 'LineWidth', 2, 'DisplayName','');
% 
%         % Open circle at base
%         hMarkers(i) = geoscatter(gx, lat_base, lon_base, 50, ...
%                                  'MarkerEdgeColor', colors{i}, ...
%                                  'MarkerFaceColor','none', 'LineWidth',1.5,'DisplayName','');
%     end
% 
%     % --- Compute difference arrow (2 - 3) ---
%     % Midpoint between instrument 2 and 3
%     lon_base_diff = (longs(2) + longs(3)) / 2;
%     lat_base_diff = (lats(2) + lats(3)) / 2;
% 
%     % Compute difference vector
%     lon_diff = (lon_tip(2) - longs(2)) - (lon_tip(3) - longs(3));
%     lat_diff = (lat_tip(2) - lats(2)) - (lat_tip(3) - lats(3));
%     lon_diffdot = [lon_diffdot lon_diff];
%     lat_diffdot = [lat_diffdot lat_diff];
%     % Tip of difference arrow
%     lon_tip_diff = lon_base_diff + lon_diff;
%     lat_tip_diff = lat_base_diff + lat_diff;
% 
%     % Plot difference arrow (black)
%     hArrows(end+1) = geoplot(gx, [lat_base_diff, lat_tip_diff], ...
%                              [lon_base_diff, lon_tip_diff], ...
%                              'Color', [0 0 0], 'LineWidth',2,'DisplayName','');
% 
%     % Open circle at base of difference arrow
%     hMarkers(end+1) = geoscatter(gx, lat_base_diff, lon_base_diff, 50, ...
%                                  'MarkerEdgeColor',[0 0 0], ...
%                                  'MarkerFaceColor','none','LineWidth',1.5,'DisplayName','');
% 
%     % Set map limits
%     geolimits(gx, [34.4653 34.4736], [-120.133 -120.125]);
% 
%     % Add timestamp title
%     title(gx, datestr(T_all{4}(idx), 'mmm dd, yyyy HH:MM'));
% 
%     drawnow
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% end
% 
% close(v);
% 
% 
% figure, scatter(lon_diffdot, lat_diffdot)
% hold on
% p = polyfit(lon_diffdot, lat_diffdot, 1);
% yfit = polyval(p, lon_diffdot);
% plot(lon_diffdot, yfit, 'r-', 'LineWidth', 2)
% title('M1 vs M2')
% 
% exportgraphics(gcf, '../../../../Kelp_data/Summer2025/Rooker/figures/Release2/difference_scatter.png')




% %% Make MP4 of Vectors on Map
% v = VideoWriter('../../../../Kelp_data/Summer2025/Rooker/figures/Release2/animations/geoplot_animation.mp4', 'MPEG-4');
% v.FrameRate = 4;
% open(v)
% 
% figure;
% gx = geoaxes;
% grid(gx, 'off');       % turn off the grid
% hold(gx, 'on');
% geobasemap satellite
% 
% % --- Create dummy arrows for persistent legend ---
% hArrowsLegend = gobjects(length(colors), 1);
% for i = 1:length(colors)
%     hArrowsLegend(i) = geoplot(gx, [0 0], [0 0], 'Color', colors{i}, 'LineWidth', 2);
% end
% legend(hArrowsLegend, labels, 'Location', 'northeast')
% 
% % --- Animate arrows and base circles ---
% for idx = 1:length(tq)
%     % Delete previous frame arrows and markers
%     if exist('hArrows','var'); delete(hArrows); end
%     if exist('hMarkers','var'); delete(hMarkers); end
% 
%     hold(gx, 'on')
%     for i = 1:length(lats)
%         if idx > length(U_all{i})
%             continue
%         end
%         lat_base = lats(i);
%         lon_base = longs(i);
%         scale = 0.03;
%         lat_tip = lat_base + V_all{i}(idx) * scale;
%         lon_tip = lon_base + U_all{i}(idx) * scale;
% 
%         % Plot arrow
%         hArrows(i) = geoplot(gx, [lat_base, lat_tip], [lon_base, lon_tip],...
%                              'Color', colors{i}, 'LineWidth', 2, 'DisplayName','');
% 
%         % Plot open circle at base (excluded from legend)
%         hMarkers(i) = geoscatter(gx, lat_base, lon_base, 50, ...
%                                  'MarkerEdgeColor', colors{i}, ...
%                                  'MarkerFaceColor','none', 'LineWidth', 1.5, 'DisplayName','');
%     end
% 
%     % Set map limits
%     geolimits(gx, [34.4653 34.4736], [-120.133 -120.125])
% 
%     % Add timestamp title
% 
%     title(gx, datestr(T_all{4}(idx), 'mmm dd, yyyy HH:MM'))
% 
%     % Capture frame
%     drawnow
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% end
% 
% close(v)
