% Pasting scripts that are going to be used for L1 proccesing

% use acceleration and jolt to filter bad data
u1   = A.b1;
d1   = gradient(u1)/A.config.dt;
dd1  = gradient(d1)/A.config.dt;
u2   = A.b2;
d2   = gradient(u2)/A.config.dt;
dd2  = gradient(d2)/A.config.dt;
u3   = A.b3;
d3   = gradient(u3)/A.config.dt;
dd3  = gradient(d3)/A.config.dt;
%
r01  =  nanstd(u1(:));
r02  =  nanstd(u2(:));
r03  =  nanstd(u3(:));
R0   = (u1./r01).^2 + (u2./r02).^2 + (u3./r03).^2;
%
r11  = nanstd(d1(:));
r12  = nanstd(d2(:));
r13  = nanstd(d3(:));
R1   = (d1./r11).^2 + (d2./r12).^2 + (d3./r13).^2;
%
r21  = nanstd(dd1(:));
r22  = nanstd(dd2(:));
r23  = nanstd(dd3(:));
R2   = (dd1./r21).^2 + (dd2./r22).^2 + (dd3./r23).^2;
%
valid = (R0<0.5) & (R1<5) & (R2<10);
A.qcFlag = A.qcFlag & valid;

% quick convolution running mean filter
np1 = round(0.3/A.config.binSize);% 10 cm vertical 
np2 = 31;% 5min for 1Hz data
f1 = hamming(np1);f1 = f1./sum(f1);
f2 = hamming(np2);f2 = f2./sum(f2);
% do a nan-mean filter, keep track of normalization
on = conv2(f1,f2,A.qcFlag','same');
%
A.A1 = conv2(f1,f2,(A.a1.*A.qcFlag)','same')./on;
A.A2 = conv2(f1,f2,(A.a2.*A.qcFlag)','same')./on;
A.A3 = conv2(f1,f2,(A.a3.*A.qcFlag)','same')./on;
%
%
A.C1 = conv2(f1,f2,(A.c1.*A.qcFlag)','same')./on;
A.C2 = conv2(f1,f2,(A.c2.*A.qcFlag)','same')./on;
A.C3 = conv2(f1,f2,(A.c3.*A.qcFlag)','same')./on;

V1 = conv2(f1,f2,(A.v1.*A.qcFlag)','same')./on;
V2 = conv2(f1,f2,(A.v2.*A.qcFlag)','same')./on;
V3 = conv2(f1,f2,(A.v3.*A.qcFlag)','same')./on;

% Deleted ENU convs but its the same form
A.Velocity_Beam1 = (conv2(f1,f2,(A.b1.*A.qcFlag)','same')./on)';
A.Velocity_Beam2 = (conv2(f1,f2,(A.b2.*A.qcFlag)','same')./on)';
A.Velocity_Beam3 = (conv2(f1,f2,(A.b3.*A.qcFlag)','same')./on)';


% get the time-averaged current.
% Note, this is not in depth normalized (sigma) coordinates.
% first, nan any values that don't pass QC.
V1(~A.qcFlag')=nan;
V2(~A.qcFlag')=nan;
V3(~A.qcFlag')=nan;
flag = sum(A.qcFlag,1) > 0.50*nsamples;
% Uz = time averegd; U = depth & time averaged
Uz = nanmean(V1,2); U = nanmean(Uz); Uz(~flag)=nan;
Vz = nanmean(V2,2); V = nanmean(Vz); Vz(~flag)=nan;
Wz = nanmean(V3,2); W = nanmean(Wz); Wz(~flag)=nan;



A.VelXAvg = Uz;
A.VelYAvg = Vz;
A.VelZAvg = Wz;


