% Unwrapping AQD Release 1 Data?

Data = load('/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release1/raw/KELP1_AquadoppHR_raw.mat');

TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");

dye = find(Data.time >= datenum(TRange.StartTime_UTC_(1)) & Data.time <= datenum(TRange.EndTime_UTC_(1)));

% use acceleration and jolt to filter bad data
% u1   = Data.b1;
% d1   = gradient(u1)/Data.Config.dt;
% dd1  = gradient(d1)/Data.Config.dt;
% u2   = Data.b2;
% d2   = gradient(u2)/Data.Config.dt;
% dd2  = gradient(d2)/Data.Config.dt;
% u3   = Data.b3;
% d3   = gradient(u3)/Data.Config.dt;
% dd3  = gradient(d3)/Data.Config.dt;
%
% r01  =  nanstd(u1(:));
% r02  =  nanstd(u2(:));
% r03  =  nanstd(u3(:));
% R0   = (u1./r01).^2 + (u2./r02).^2 + (u3./r03).^2;
% %
% r11  = nanstd(d1(:));
% r12  = nanstd(d2(:));
% r13  = nanstd(d3(:));
% R1   = (d1./r11).^2 + (d2./r12).^2 + (d3./r13).^2;
% %
% r21  = nanstd(dd1(:));
% r22  = nanstd(dd2(:));
% r23  = nanstd(dd3(:));
% R2   = (dd1./r21).^2 + (dd2./r22).^2 + (dd3./r23).^2;
%
% valid = (R0<0.5) & (R1<5) & (R2<10);
% Data.qcFlag = valid;


Data.qcFlag = 0*Data.b1;

 Velocities = {Data.b1(dye,:), Data.b2(dye,:), Data.b3(dye,:)};
% MaxV = max(max(Velocities{1}(:)), max(max(Velocities{2}(:)), max(Velocities{3}(:))));
% MinV = abs(min(min(Velocities{1}(:)), min(min(Velocities{2}(:)), min(Velocities{3}(:)))));
% Vrange = (MaxV+MinV)/2;


Vrange = (0.72*sind(25)+0.31*cos(25))/2;

for i = 1:3
    Velocity = Velocities{i};
    u   = Velocity;
    d   = gradient(u,1);
    dd  = gradient(d,1);
    Data.qcFlag = d >= 0.0001;% | dd >= Vrange/10;
    for j = dye

        for k = 1:75

            if Data.qcFlag(j, k)
                if Velocity(j, k) > 0
                    Velocity(j, k) = Velocity(j, k) - Vrange;
                else Velocity(j, k) < 0;
                    Velocity(j, k) = Velocity(j, k) + Vrange;
                end
            end
        end
    end
    Velocities{i} = Velocity;
end


figure, histogram(Velocities{i}(dye, 1), 'BinWidth', 0.01)









