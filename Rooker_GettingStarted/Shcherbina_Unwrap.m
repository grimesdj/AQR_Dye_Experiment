


% Testing the Shcherbina et al 2018 unwrapper
%
clear
%

 Data = load('../../../Kelp_data/Summer2025/Rooker/Release1/raw/KELP1_AquadoppHR_raw.mat');
% 
releaseNum = 1
TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");
dye = find(Data.time >= datenum(TRange.StartTime_UTC_(releaseNum)) & Data.time <= datenum(TRange.EndTime_UTC_(releaseNum)));



% 
% 

% for beam = 1:1
for bin = 1:1
 v_wrapped = Data.b1(dye, bin);
 
 Vr = ((Data.sspeed(dye, :).^2)/(4*1000^2*2.5));
 
 v_unwrap = Shcherbina_Unwrapper(v_wrapped, Vr);

%  Data.(sprintf('b%d',beam))(dye,:) = v_unwrap;
%  figure, plot(Data.(sprintf('b%d',beam))(dye,:), '.')
% end

unwraps(:,bin)= v_unwrap;
end

%save('../../../Kelp_data/data/Release1/raw/KELP1_AquadoppHR_unwrap.mat', '-struct', 'Data')
 
function [v_unwrap] = Shcherbina_Unwrapper(v_wrapped, Vr)

[nbins, nt] = size(v_wrapped);
v_unwrap = v_wrapped;

filt_diff = (v_wrapped - medfilt2(v_wrapped,[5,3]))./Vr;
filt_diff1 = filt_diff;
suspect_pts = abs(filt_diff)>0.5;

figure, plot(v_wrapped, '.')
hold on, plot(v_wrapped(suspect_pts), 'r.')
hold on, plot(v_wrapped(~suspect_pts), 'g.')

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

% plots;
% figure(1),clf
% subplot(3,1,1)
% pcolor((v_wrapped(:,1:ncol)./mean(Vr))');
% shading flat
% title('Original')
% colorbar
% caxis([-1,1])
% 
% subplot(3,1,2)
% pcolor((v_wrapped(:,1:ncol)./mean(Vr) - filt_diff1(:,1:ncol))');
% shading flat
% title('Filtered')
% colorbar
% caxis([-1,1])
% 
% subplot(3,1,3)
% pcolor((v_unwrap(:,1:ncol)./mean(Vr))');
% shading flat
% title('Unwrap')
% colorbar
% caxis([-1,1])
end

