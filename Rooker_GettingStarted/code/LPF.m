function [T_all, U_all, V_all] = LPF(dataCell, labels, colors)
%   
% USAGE: [T_all, U_all, V_all] = LPF(dataCell, labels, colors)
% 
%   takes time series data and 
%  returns 3 cell arrays, containing 10-min averaged data, u, v, and time
% 
%   dataCell = cell array (# of files x 1) containing data structures
%   labels = {optional} Cell array of labels for graphs
%   colors = {optional} Cell array of RGB vals for graphs
%    
%
%   T_all = interpolated time
%   U_all = East Velocities
%   V_all = North Velocities
%   

%% Defaults
if nargin < 3
    colors = num2cell(lines(length(dataCell)), 2);
    if nargin < 2
        labels = num2cell(1:length(dataCell));
    end
end


%% Time Series Plot Setup
figure;
ax1 = subplot(2, 1, 1); hold(ax1, 'on'); ylabel(ax1, 'East Velocity');
ax2 = subplot(2, 1, 2); hold(ax2, 'on'); ylabel(ax2, 'North Velocity'); xlabel(ax2, 'Time');
datetick(ax1, 'x', 'keeplimits'); datetick(ax2, 'x', 'keeplimits');
sgtitle('10-min Avg Velocities at 1.2 MAB')

%% Quiver Plot Setup
figure;
for i = 1:length(dataCell)

    data = dataCell{i};

    % Extract and clean data
    u = data.Velocity_East(:, data.bin);
    v = data.Velocity_North(:, data.bin);
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
    subplot(length(dataCell), 1, i)
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
for i = 1:length(dataCell)
    h1(i).Color = colors{i};
    h2(i).Color = colors{i};
end
yline(ax1, 0); yline(ax2, 0);

end