% Testing the Shcherbina et al 2018 unwrapper
%
%
%

% clear
% 
% load('../Data/TiogaDockTest_23Oct2020/TiogaTest_mst.297.00000.ad2cp.00000.mat');
% 
% v_wrapped = Data.BurstHR_VelBeam2';

function [v_unwrap] = Shcherbina_Unwrap(v_wrapped, Vr)

[nbins, nt] = size(v_wrapped);
v_unwrap = v_wrapped;

filt_diff = (v_wrapped - medfilt2(v_wrapped,[5,3]))./mean(Vr);
%filt_diff1 = filt_diff;
suspect_pts = abs(filt_diff)>1.1;

for ncol = 1:nt
%ncol = 30;

si = find(suspect_pts(:,ncol));
%create difference operator
E = eye(nbins);
D = diff(E);
E = E(:,si);
v_prime = D*v_wrapped(:,ncol);

%solve least squares problem, and correct profile
r = round( (2*Vr(ncol)*D*E)\v_prime );
v_unwrap(si, ncol) = v_wrapped(si, ncol) - r*2*Vr(ncol);

end

% - plots;
% figure(1),clf
% subplot(3,1,1)
% pcolor(v_wrapped(:,1:ncol)./mean(Vr));
% shading flat
% title('Original')
% colorbar
% caxis([-1,1]*5)
% 
% subplot(3,1,2)
% pcolor(v_wrapped(:,1:ncol)./mean(Vr) - filt_diff1(:,1:ncol));
% shading flat
% title('Filtered')
% colorbar
% caxis([-1,1]*5)
% 
% subplot(3,1,3)
% pcolor(v_unwrap(:,1:ncol)./mean(Vr));
% shading flat
% title('Unwrap')
% colorbar
% caxis([-1,1]*5)
