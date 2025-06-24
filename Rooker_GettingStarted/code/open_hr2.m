
releaseNum = 1


inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat"];

Data = load(inputFiles(releaseNum));

TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");
dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.time <= datenum(TRange.EndTime_UTC_(releaseNum)));

fid = fopen("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR.hr2");

% Format: 6 integers/floats for date & time, 3 floats, 9 integers
formatSpec = '%d %d %d %d %d %f %f %f %f %d %d %d %d %d %d %d %d %d';

data = textscan(fid, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fid);

dt = datetime(data{3}, data{1}, data{2}, data{4}, data{5}, data{6});

sensorValues = [data{7}, data{8}, data{9}];  % The 3 float values
a1 = data{10};
a2 = data{11};
a3 = data{12};

c1 = data{13};
c2 = data{14};
c3 = data{15};
            
C = [c1, c2, c3];
A = [a1, a2, a3];
% Plot each sensor value as a separate time series
figure;

% Sensor labels (optional)
sensorLabels = {'Beam 1', 'Beam 2', 'Beam 3'};

for i = 1:3
    figure;
    subplot(3, 1, 1);  % Create a 3-row subplot
    plot(dt(dye), sensorValues(dye, i), 'b.');
    ylabel('Magnitude');
    title(sensorLabels{i});
    grid on;
    subplot(3, 1, 2);
    plot(dt(dye), A(dye, i), 'r.');
    subplot(3, 1, 3);
    plot(dt(dye), C(dye, i), 'g.');
    ylim([0, 100])
    figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\timeseries\hr2_ts_r" + releaseNum + "_beam" + i + ".pdf"];
    %exportgraphics(gcf, figname);
    %close(gcf)
end

xlabel('Time');


sensorLabels = {'Beam 1', 'Beam 2', 'Beam 3'};

for i = 1:3
    figure;
    histogram(sensorValues(dye, i), 'BinWidth', 0.01);
    title(['Histogram of ', sensorLabels{i}]);
    xlabel('Magnitude');
    ylabel('Count');
    grid on;
    xlim([-0.8 0.8])
    figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\histogram\hr2_hist_r" + releaseNum + "_beam" + i + ".pdf"];
    %exportgraphics(gcf, figname);
    %close(gcf)
end

addpath 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_repo\AQR_Dye_Experiment\Rooker_GettingStarted\scratch'
ellipse 