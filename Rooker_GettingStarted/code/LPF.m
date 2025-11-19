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
tp = figure;
ax1 = subplot(2, 1, 1); hold(ax1, 'on'); ylabel(ax1, 'East Velocity');
ax2 = subplot(2, 1, 2); hold(ax2, 'on'); ylabel(ax2, 'North Velocity'); xlabel(ax2, 'Time');
datetick(ax1, 'x', 'keeplimits'); datetick(ax2, 'x', 'keeplimits');
sgtitle('10-min Avg Velocities at 1.2 MAB')

%% Quiver Plot Setup
qp = figure;
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
      disp('UNDERCONSTRUCTION')  
        %fixing NaNs with linear interp
        ids = 1:length(data.Velocity_East);
        valid = ~isnan(u);
        u_interp = u;
        v_interp = v;
        u_interp(~valid) = 
        v_interp(~valid) = 
    
    flag_flt = conv(nanFlag, flt, 'same');
    u(~nanFlag) = 0; v(~nanFlag) = 0;
    u_filt = conv(u, flt, 'same') ./ flag_flt;
    v_filt = conv(v, flt, 'same') ./ flag_flt;

    % Plot time series
    figure(1);
    h1(i) = plot(ax1, t, u_filt,'LineWidth', 1);
    h2(i) = plot(ax2, t, v_filt,'LineWidth', 1);
    
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
    
    % Time relative to start, in seconds
    tq_rel = seconds(tq - tq(1));

    % Convert velocity to displacement per 10 min
    dt = 600*100*3;  % scaling
    u_disp = u_ds * dt;
    v_disp = v_ds * dt;

    % Plot quiver
    quiver(tq_rel, zeros(size(tq_rel)), u_disp, v_disp, 'AutoScale', 'off')
    axis equal  % ensure same scale in x and y for correct arrow directions
    yticklabels([])
    xticklabels([])
    title(data.Config.SN)
 
    % Save Global Vars
    U_all{i} = u_ds;
    V_all{i} = v_ds;
    T_all{i} = tq;
    
end
xlabel('Time')
sgtitle('Current Vectors at 1.2 MAB')

quivers = get(qp, 'Children');
linkaxes(quivers, 'x')

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