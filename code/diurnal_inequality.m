clear all
close all

fin = '~/OneDriveUNCW/KELP-vingilote/data/2024_PROCESSED_DATA/M1/L1/ADCP/ADCP_M1_L1_currents.mat';
load(fin)


S2 = 12;
M2 = 12.42;
K1 = 23.93;
O1 = 25.82;

t = currents.Time;
t = (t-t(1))'*24;

p = currents.Pressure';

ig = ~isnan(p);

t = t(ig);
p = p(ig);
p = detrend(p);


figure, plot(datetime(t/24+currents.Time(1),'convertFrom','datenum'),p,'-k','linewidth',2)

XX = [1+0*t, sin(2*pi/S2*t), cos(2*pi/S2*t), sin(2*pi/M2*t), cos(2*pi/M2*t), sin(2*pi/K1*t), cos(2*pi/K1*t),sin(2*pi/O1*t), cos(2*pi/O1*t)];

cff = XX\p;


A_S2 = sqrt( cff(2).^2 + cff(3).^2 )
A_M2 = sqrt( cff(4).^2 + cff(5).^2 )
A_K1 = sqrt( cff(6).^2 + cff(7).^2 )
A_O1 = sqrt( cff(8).^2 + cff(9).^2 )


hold on,
plot(datetime(t/24+currents.Time(1),'convertFrom','datenum'),XX*cff,'--r','linewidth',2)

% $$$ plot(datetime(t/24+currents.Time(1),'convertFrom','datenum'),XX(:,1:5)*cff(1:5),'-b','linewidth',2)
% $$$ 
% $$$ plot(datetime(t/24+currents.Time(1),'convertFrom','datenum'),XX(:,6:9)*cff(6:9),'-m','linewidth',2)
