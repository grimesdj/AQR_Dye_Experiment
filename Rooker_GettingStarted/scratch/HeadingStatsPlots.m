releaseNum = 1
inputFiles = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\raw\KELP2_AquadoppHR_raw.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat";
              "C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Vector_L0.mat"];
          

Data = load(inputFiles(releaseNum));


TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.time <= datenum(TRange.EndTime_UTC_(releaseNum+1)));

% Get L0
L0.Mag = [(Data.east(:,34).^2 + Data.north(:,34).^2).^0.5];
L0.theta = atan2d(Data.east(:,34), Data.north(:,34));
% Get M1
M1 = fetch_M1(releaseNum);
M1.Mag = [(M1.east(:, 11).^2 + M1.north(:, 11).^2).^0.5];
M1.theta = atan2d(M1.east(:,11), M1.north(:,11));
% Get Vector
releaseNum = releaseNum+3;
    AQD = Data;
    Data = load(inputFiles(releaseNum));    
    rotation
    releaseNum = releaseNum-3;

Vector = Data;
Vector.Mag = [(Data.east.^2 + Data.north.^2).^0.5];
Vector.theta = atan2d(Data.east, Data.north);

u1 = L0.theta;
DU = gradient(u1);
L0Flag = find(abs(DU)>180);
u1 = M1.theta;
DU = gradient(u1);
M1Flag = find(abs(DU)>180);
u1 = Vector.theta;
DU = gradient(u1);
VectorFlag = find(abs(DU)>180);

% Unwrap (I hope this works)

            
        L0.theta(i+1) = 2*L0.theta(i) + L0.theta(i+1);
    
        M1.theta(i+1) = 2*M1.theta(i)+ M1.theta(i+1);
    
        Vector.theta(i+1) = 2*Vector.theta(i)+ Vector.theta(i+1);
    


% Plot 

tsf = figure;
ax1 = subplot(2,1,1);
plot(L0.theta(dye, :), 'b')
% ylim(ax1,[-360 360]);
ylabel(ax1,'direction');
hold on
plot(M1.theta(dye, :), 'r')
plot(Vector.theta(dye, :), 'g')
hold off

ax2 = subplot(2,1,2);
plot(L0.Mag(dye, :), 'b')
hold on
plot(M1.Mag(dye, :), 'r')
plot(Vector.Mag(dye, :), 'g')
ylabel(ax2,'magnitude');
hold off






