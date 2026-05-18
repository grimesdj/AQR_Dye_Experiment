
function L0_AQD(outputFile, L0Dir, L0Name)
% 
%   USAGE: L0_AQD(outputFile, L0Dir, L0Name)
%       outputFile = folder path and filename (no extension) to raw data
%       L0Dir = directory to save finished L0 data
%       L0Name = L0 filename (without extension) to be saved
% 
%       takes raw AQD data and performs L0 QA/QC
% 
A = load([outputFile, '.mat']);
fprintf('\n============================\nDo you want to unwrap beam Velocities?')
unwrap = input('(1 = yes; 0 = no)');
if unwrap == 1
    for beam = 1:3 
        [A.(sprintf('Velocity_Beam%d',beam)), A.(sprintf('Suspect_Beam%d', beam))] = unwrap_AQD(A.(sprintf('Velocity_Beam%d',beam)), A.VRange);
    end
    
    A.Correlation_Beam1(find(A.Suspect_Beam1)) = 999;
    A.Correlation_Beam2(find(A.Suspect_Beam2)) = 999;
    A.Correlation_Beam3(find(A.Suspect_Beam3)) = 999;

% rotate unwrapped data!
    % XYZ
    disp('rotating unwrapped data to XYZ')
    shape = size(A.Velocity_Beam1);
    XYZ= A.Config.transform_matrix*[A.Velocity_Beam1(:)'; A.Velocity_Beam2(:)'; A.Velocity_Beam3(:)'];
    A.Velocity_X = reshape(XYZ(1,:)',shape);
    A.Velocity_Y = reshape(XYZ(2,:)',shape);
    A.Velocity_Z = reshape(XYZ(3,:)',shape);

    % ENU
    disp('rotating unwrapped data to ENU')
    hh = reshape(pi*(A.Heading-90)/180,1,1,A.Config.Nsamples);
    pp = reshape(pi*A.Pitch/180,1,1,A.Config.Nsamples);
    rr = reshape(pi*A.Roll/180,1,1,A.Config.Nsamples);
    H = [ cos(hh), sin(hh), 0*hh;...
         -sin(hh), cos(hh), 0*hh;...
          0*hh,      0*hh,  0*hh+1];
    P = [cos(pp), 0*pp, -sin(pp);...
          0*pp,  0*pp+1, 0*pp   ;...
         sin(pp),  0*pp,  cos(pp)];
     
    O = [1+0*rr 0*rr 0*rr;...
        0*rr cos(rr) -sin(rr);...
        0*rr sin(rr) cos(rr)];
     
     
    
    for j = 1:A.Config.Nsamples
        R   = H(:,:,j)*P(:,:,j)*O(:,:,j);
        ENU = R*[A.Velocity_X(j,:);A.Velocity_Y(j,:);A.Velocity_Z(j,:)];
        A.Velocity_East  (j,:) = ENU(1,:);
        A.Velocity_North (j,:) = ENU(2,:);
        A.Velocity_Up    (j,:) = ENU(3,:);
    end

end

%

A.dbins = (A.Config.blank + A.Config.binSize*A.Config.Nbins);
A.maxRange = (A.Pressure-A.pressureOffset).*cosd(20)-1*A.Config.binSize;
A.ylims      = [0 min(max(A.maxRange),max(A.dbins))];
dum1       = A.maxRange.*ones(1,A.Config.Nbins);
dum2       = ones(A.Config.Nsamples,1)*A.dbins;
qcFlag0    =  (dum2<=dum1);
A.Amplitude_Minimum   = min(A.Amplitude_Beam1, min(A.Amplitude_Beam2, A.Amplitude_Beam3));
A.Correlation_Minimum = min(A.Correlation_Beam1, min(A.Correlation_Beam2, A.Correlation_Beam3, 'omitnan'), 'omitnan');   
A.qcFlag              =  double( qcFlag0 & A.Amplitude_Minimum > 20 & A.Correlation_Minimum > 40 );

% QCFlag
disp('Applying qcFlag to NaN data A < 20 and C < 40')

%

A.Velocity_East = A.Velocity_East.*A.qcFlag;
A.Velocity_North = A.Velocity_North.*A.qcFlag;
A.Velocity_Up = A.Velocity_Up.*A.qcFlag;

A.Velocity_East(~A.qcFlag)=nan;
A.Velocity_North(~A.qcFlag)=nan;
A.Velocity_Up(~A.qcFlag)=nan;

A.Velocity_Beam1 = A.Velocity_Beam1.*A.qcFlag;
A.Velocity_Beam2 = A.Velocity_Beam2.*A.qcFlag;
A.Velocity_Beam3 = A.Velocity_Beam3.*A.qcFlag;

A.Velocity_Beam1(~A.qcFlag)=nan;
A.Velocity_Beam2(~A.qcFlag)=nan;
A.Velocity_Beam3(~A.qcFlag)=nan;

A.Velocity_X = A.Velocity_X.*A.qcFlag;
A.Velocity_Y = A.Velocity_Y.*A.qcFlag;
A.Velocity_Z = A.Velocity_Z.*A.qcFlag;

A.Velocity_X(~A.qcFlag)=nan;
A.Velocity_Y(~A.qcFlag)=nan;
A.Velocity_Z(~A.qcFlag)=nan;

%

disp('skipping nc file for now')

%

% Summary figure
figure
ax1 = subplot(3, 1, 1);
plot(ax1, A.Time, A.Velocity_Beam1, '.')
ylabel({'Beam 1', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor
ylim([-1 1])

ax2 = subplot(3, 1, 2);
plot(ax2, A.Time, A.Velocity_Beam2, '.')
ylabel({'Beam 2', 'Velocity, [m/s]'})
set(gca, "Xtick", [])
set(gca, 'fontsize', 18)
grid minor
ylim([-1 1])

ax3 = subplot(3, 1, 3);
plot(ax3, A.Time, A.Velocity_Beam3, '.')
ylabel({'Beam 3', 'Velocity, [m/s]'})
datetick(gca,'x','mmm-dd HH:MM','keeplimits')
set(gca, 'fontsize', 18)
linkaxes([ax1 ax2 ax3], 'x')
ylim([-1 1])

sgtitle('AQD L0 Beam Velocities', 'Fontsize', 25)

disp('Saving L0 data')
save([L0Dir,'/',L0Name,'.mat'],'-struct','A')

end