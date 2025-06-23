function multiplot(releaseNum, qc)

version = 'raw';

if nargin < 2
    qc = 0;
elseif strcmp(qc,'on')
    qc = 1;
else
    qc = 0;
end

% fetch Data
inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat"];
          

Data = load(inputFiles(releaseNum));


TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.time <= datenum(TRange.EndTime_UTC_(releaseNum)));

% Get L0
    Time = Data.time;

if qc
    version = 'L0';
   minAmp = min(Data.a1, min(Data.a2, Data.a3));
    minCor = min(Data.c1, min(Data.c2, Data.c3));
    flagA = find(minAmp <= 30);
    flagC = find(minCor <= 30);
    flagind = unique([flagA' flagC']);
    flagind; %iterate through bins smh
    
    Data.a1(flagind) = NaN;
    Data.a2(flagind) = NaN;
    Data.a3(flagind) = NaN;
    Data.b1(flagind) = NaN;
    Data.b2(flagind) = NaN;
    Data.b3(flagind) = NaN;
    Data.c1(flagind) = NaN;
    Data.c2(flagind) = NaN;
    Data.c3(flagind) = NaN;
    Data.east(flagind) = NaN;
    Data.north(flagind) = NaN;
    Data.up(flagind) = NaN;
end

L0.East = Data.east(:,35);
L0.North = Data.north(:,35);
fraction_of_goood_AQD_data = sum(~isnan(L0.East(dye)))/numel(L0.East(dye))

% Get M1

M1 = fetch_M1(releaseNum);

M1.Einterp = interp1(M1.time', M1.east(:, 6), Time);
M1.Ninterp = interp1(M1.time', M1.north(:, 6), Time);

% Get Vector
Vector = load(inputFiles(releaseNum+3));
= rotation
Vector.East = interp1(Vector.time, Vector.east, Time);
Vector.North = interp1(Vector.time, Vector.north, Time);




% Plot 

tsf = figure;
ax1 = subplot(2,1,1);
plot(Time(dye), L0.East(dye, :), 'b')
ylim(ax1,[-1.5 1.5]);
ylabel(ax1,'East');
yline(-0.72/2)
yline(0.72/2)

hold on
Mplot = plot(Time(dye), M1.Einterp(dye), 'r');
plot(Time(dye), Vector.East(dye), 'g')
uistack(Mplot, 'bottom')
hold off

ax2 = subplot(2,1,2);
plot(Time(dye), L0.North(dye), 'b')

hold on
Mplot = plot(Time(dye), M1.Ninterp(dye), 'r');
plot(Time(dye), Vector.North(dye), 'g')
ylim(ax2,[-1.5 1.5]);
ylabel(ax2,'North');
yline(-0.72/2)
yline(0.72/2)
uistack(Mplot, 'bottom')

hold off
linkaxes([ax1,ax2],'x')



% Export Plot
figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\ENU\" + version + "_ENU_r" + releaseNum + ".pdf"];
exportgraphics(gcf, figname);


end

