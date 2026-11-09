%% wind climate

clear all
close all

%% Load
files = dir('../../../../Kelp_data/data/wind/*_wind_data_nearshore.mat');

fname = fullfile(files(1).folder, files(1).name);
A = load(fullfile(files(1).folder, files(1).name));
fields = fieldnames(A);

for j = 2:length(files)
    fname = fullfile(files(j).folder, files(j).name);
    tmp = load(fullfile(files(j).folder, files(j).name));
    for f = 1:length(fields)
        field = fields{f};
        A.(field) = [A.(field);tmp.(field)];
    end
end

% Generates names from January to December
monthNames = month(datetime(1, 1:12, 1), 'name');

% convert from UTC
% A.Time = A.Time - hours(7);
% A.MM_mo = month(A.Time);
% A.hh_hr = hour(A.Time);
% A.mm_mn = minute(A.Time);

%% July Climate ( make a loop later )
for m = 1:12
    isJuly = (A.MM_mo == m);
    indJuly = find(isJuly);
    
    for f = 1:length(fields)
        field = fields{f};
        if length(A.(field)) ~= length(A.MM_mo)
            A.(field)(length(A.(field)):length(A.MM_mo)) = NaN;
        end
        C.(field) = A.(field)(indJuly);
    end
    
    
    %% hour of day
    str = [string(C.hh_hr) + ':' + string(C.mm_mn)];
    C.HOD = datetime(str, 'InputFormat', 'HH:mm');

    figure
    ax1 = subplot(2, 1, 1);
    plot(C.HOD, C.North, '.')
    ylabel({'North', '[m/s]'})
    ylim([-10 10])
    axis square
    hold on

    ax2 = subplot(2, 1, 2);
    plot(C.HOD, C.East, '.')
    ylabel({'East', '[m/s]'})
    ylim([-10 10])
    axis square
    hold on

    sgtitle(sprintf('%s', monthNames{m}))


    startINDs = [1;find(diff(datenum(C.HOD)) < -0.5)];
    endINDs = [startINDs(2:end); length(C.HOD)];
    
    u = unique(C.HOD);
    nt = length(u);
    for f = 1:length(fields)
        field = fields{f};
        if isa(C.(field), 'double')
            for i = 1:length(startINDs)-2
                inds = indJuly(startINDs(i))+[0:nt-1]';
                S(m).(field)(:, i) = A.(field)(inds, :);
            end
        end
    end

    N = S(m).North;
    E = S(m).East;

    Nbar = mean(N, 2, 'omitmissing');
    plot(ax1, u, Nbar, 'r', 'LineWidth', 3)
    yline(ax1, 0, 'k--', 'LineWidth', 1)

    Ebar = mean(E, 2, 'omitmissing');
    plot(ax2, u, Ebar, 'r', 'LineWidth', 3)
    yline(ax2, 0, 'k--', 'LineWidth', 1)

    TOD = u;

    savestr = sprintf('../../../../EGN499/data/%s_wind_climate.mat', monthNames{m});
    save(savestr, "N", "E", "Nbar", "Ebar", "TOD", '-mat')

end
