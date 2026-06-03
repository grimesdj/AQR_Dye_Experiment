%% load and plot RBR

clear all
close all

%% find files

fname = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', 'FullExperiment', 'raw', 'RBR', '*.rsk');
files1 = dir(fname);
fname = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', 'FullExperiment', 'raw', 'RBR_7_26_24', '*.rsk');
files2 = dir(fname);
keyboard

files = cat(1, files1, files2);
Rtimes = string.empty;




% for each mooring...
for i = 1:4

    Moor = load("../../../../Kelp_data/data/2024_PROCESSED_DATA/M" + i + "/L1/mooring_M" + i + ".mat");
    count = 1;

    Names.(Moor.Name) = string.empty;

    figure
    hold on
    grid on
    %ax1 = subplot(2, 1, 1);
    %ax2 = subplot(2, 1, 2);
    % ...find SNs...
    SNs = Moor.Instrument_SN;
    for j = 1:length(SNs)
        % ... grab a single SN...
        SN = num2str(SNs(j));
        for k = 1:length(files)
            % ...that matches rsk SN...
            if contains(files(k).name, SN)
                % ...and load it...
                fprintf('loading %s...\n', files(k).name)
                fin = fullfile(files(k).folder, files(k).name);
                [RSK_load, flg] = quick_load_rbr(fin);
                if ~flg

                    RSK(count) = RSK_load;

                    % ...and plot it...
                    plot(datetime(RSK(count).Time, 'convertfrom', 'datenum'), RSK(count).Temp, '-', 'LineWidth', 1, 'DisplayName', files(k).name)
                    drawnow
                    % count and save names
                    Names.(Moor.Name)(end+1) = string(files(k).name);
                    count = count + 1;
                else
                    fprintf('%s is not a thermister\n', files(k).name(1:end-4))
                end
            end
        end
    end    
    dum = length(Moor.Temperature_mab);
    fprintf('found %d of %d ins\n', count-1, dum)
    hold on
    xline(datetime(Moor.Time_deploy, 'InputFormat', 'yyyy/MM/dd HH:mm:ss'), 'k--', 'linewidth', 2, 'label', sprintf('%s release time', Moor.Name), 'FontSize', 18, 'labelorientation', 'horizontal')
    title(sprintf('%s Temperature', Moor.Name), 'FontSize', 20)
    legend

    % collect recorded release times
    Rtimes = [Rtimes; Moor.Time_deploy];
end


%% Check instruments that aren't in the list

fields = fieldnames(Names);
allSN = string.empty;
% Get all SNs that were listed in moorings...
for i = 1:length(fields)
    field = fields{i};
    allSN = [allSN Names.(field)];
end

allSN = allSN';
count = 1;
figure
hold on

% ... for each file...
for i = 1:length(files)
    % ...if it wasn't in the list...
    if ~any(contains(allSN, files(i).name))
        % ... load it...
        fin = fullfile(files(i).folder, files(i).name);
        [unlisted_load, flg] = quick_load_rbr(fin);
        % ... check if it was temp ....
        if ~flg
            unlisted(count) = unlisted_load;
            % ... plot.
            plot(datetime(unlisted(count).Time, 'convertfrom', 'datenum'), unlisted(count).Temp, '-', 'LineWidth', 1, 'DisplayName', files(i).name)
            drawnow
           
            count = count + 1;
        else
            fprintf('%s is not a thermister\n', files(i).name(1:end-4))
        end
    end
end

% plot recorded release times on unlisted figure
for i = 1:length(Rtimes)
    xline(datetime(Rtimes(i), 'InputFormat', 'yyyy/MM/dd HH:mm:ss'), 'k--', 'Label', sprintf('M%d', i), 'LineWidth', 2)
end