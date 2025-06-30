% Unwrapping AQD Release 1 Data?

Data = load('/Users/jkr6136/OneDrive - UNC-Wilmington/Kelp_data/Summer2025/Rooker/Release1/L0/KELP1_AquadoppHR_L0.mat');


% use acceleration and jolt to filter bad data
% u1   = Data.Velocity_Beam1;
% d1   = gradient(u1)/Data.Config.dt;
% dd1  = gradient(d1)/Data.Config.dt;
% u2   = Data.Velocity_Beam2;
% d2   = gradient(u2)/Data.Config.dt;
% dd2  = gradient(d2)/Data.Config.dt;
% u3   = Data.Velocity_Beam3;
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


Data.qcFlag = 0*Data.Velocity_Beam1;

Velocities = [Data.Velocity_Beam1 Data.Velocity_Beam2 Data.Velocity_Beam3];
MaxV = max(Velocities(1), max(Velocities(2), Velocities(3)));
MinV = min(Velocities(1), min(Velocities(2), Velocities(3)));
Vrange = (MaxV+MinV)/2;

for i = 1:3
    Velocity = Velocities(i);
    u   = Velocity;
    d   = gradient(u)/Data.Config.dt;
    dd  = gradient(d)/Data.Config.dt;
    Data.qcFlag = d > Vrange;
    for j = 1:length(Data.Time)

        for k = 1:Data.Config.Nbins

            if ~Data.qcFlag(k, j)
                if Velocity(k, j) > 0
                    Velocity(k, j) = Velocity(k, j) - 2*Vrange;
                else Velocity(k,j) < 0;
                    Velocity(k, j) = Velocity(k, j) + 2*Vrange;
                end
            end
        end
    end
    Velocities(i) = Velocity;
end














