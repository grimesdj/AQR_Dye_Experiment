% Let's make some 10 min averages:
clear all
close all


% Define the files to be loaded
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release2/L0/*.mat');
colors = {[0, 0, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0]};

% Make figure
figure;
ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
datetick(ax1, 'x', 'keeplimits')
datetick(ax2, 'x', 'keeplimits')
sgtitle('10 min avg velocities at 1.2 MAB')
ylabel(ax1, 'East Velocity')
ylabel(ax2,'North Velocity')
xlabel(ax2, 'Time')

figure;
for i = 1:length(files)
    data = load(fullfile(files(i).folder, files(i).name));

    % Choose bin (same as before)
    if size(data.Velocity_East, 2) == 1
        bin = 1;
    else
        bin = round((1.237 - data.Config.blank) / data.Config.binSize);
    end
    %bin = max(1, min(bin, size(data.Velocity_East, 2)));

    % Extract and clean data
    u = data.Velocity_East(:,bin);
    v = data.Velocity_North(:,bin);
    t = data.Time;

    % Filter using same lowpass
    dt = double(data.Config.dt);
    Nf = 600 / dt;
    flt = hamming(Nf); flt = flt/sum(flt);

    nanFlag = ~isnan(u);
    flag_flt = conv(nanFlag, flt, 'same');
    u(~nanFlag) = 0; v(~nanFlag) = 0;

    u_filt = conv(u, flt, 'same') ./ flag_flt;
    v_filt = conv(v, flt, 'same') ./ flag_flt;

    % Plot and store handles and labels
    figure(1);
    hold(ax1, 'on')
    hold(ax2, 'on')
    h1(i) = plot(ax1, data.Time, u_filt);
    h2(i) = plot(ax2, data.Time, v_filt);
    labels{i} = data.Config.SN;

    % Downsample every 10 minutes
    [~, unique_idx] = unique(minutes(t - t(1))); % remove duplicate timestamps
    t = t(unique_idx);
    u_filt = u_filt(unique_idx);
    v_filt = v_filt(unique_idx);

    % Resample to 10-minute intervals
    tq = t(1):minutes(10):t(end);
    u_ds = interp1(t, u_filt, tq, 'linear');
    v_ds = interp1(t, v_filt, tq, 'linear');

    % Quiver plot setup
    figure(2)
    subplot(length(files),1,i)
    quiver(datenum(tq), zeros(size(tq)), u_ds, v_ds, 'AutoScale', 'off')
    ylabel('Velocity')
    title(sprintf('%s', data.Config.SN))
    datetick('x', 'keeplimits')
end
xlabel('Time', 'FontSize', 18)
sgtitle('Current Vectors at 1.2 MAB', 'Fontsize', 25)

% Add legends
figure(1)
legend(ax1, h1, labels, 'Location', 'best')
legend(ax2, h2, labels, 'Location', 'best')
datetick(ax1)
datetick(ax2)

for i = 1:numel(colors)
    h1(i).Color = colors{i};
    h2(i).Color = colors{i};
end

yline(ax1,0)
yline(ax2,0)

%% Map Plots

% These are NOT the correct lat/lons. having a hard time getting them from
% the info files so I approximated them for now
return
% Assign Coords
lats = [34.4690789 34.469615 34.4690706 34.4690789];
longs = [-120.1278549 -120.130375 -120.1278088 -120.1278549];

% Make a map
figure, 
geobasemap satellite
hold on
geoplot(lats(1), longs(1), 'bo')
geoplot(lats(2), longs(2), 'ro')
geoplot(lats(3), longs(3), 'mo')
geoplot(lats(4), longs(4), 'g.')

% Set up video writer
v = VideoWriter('geoplot_animation.mp4', 'MPEG-4');
v.FrameRate = 5;  
open(v)

for idx = 1:length(data.Time)
    for i = 1:length(lats)

% Define arrow scaling (adjust as needed to make arrows visible on map)
scale = 0.05;

% Arrow base coordinates (lat/lon of instrument 1)
lat_base = lats(i);
lon_base = longs(i);

% Compute the tip of the arrow
lat_tip = lat_base + v_ds(idx) * scale; % v = north/south
lon_tip = lon_base + u_ds(idx) * scale; % u = east/west

% Plot the arrow as a line
geoplot([lat_base, lat_tip], [lon_base, lon_tip], 'Color', colors{i}, 'LineWidth', 2)
     
    end
    % Capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end