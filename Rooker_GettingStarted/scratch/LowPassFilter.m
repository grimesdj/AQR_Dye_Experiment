% Let's make some 10 min averages:
clear all
close all


% Define the files to be loaded
files = dir('../../../../Kelp_data/Summer2025/Rooker/Release1/L0/*.mat');

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
xlabel('Time')
sgtitle('Current Vectors at 1.2 MAB')

% Add legends
figure(1)
legend(ax1, h1, labels, 'Location', 'best')
legend(ax2, h2, labels, 'Location', 'best')
datetick(ax1)
datetick(ax2)

h1(1).Color = 'b';
h2(1).Color = 'b';

h1(2).Color = 'r';
h2(2).Color = 'r';

h1(3).Color = 'm';
h2(3).Color = 'm';

h1(4).Color = 'g';
h2(4).Color = 'g';

yline(ax1,0)
yline(ax2,0)
