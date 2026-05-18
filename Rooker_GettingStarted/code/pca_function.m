function pca_function(Dir, FileName)
% pca_function performs principal component analysis on 2D velocity data
%   USAGE: pca_function(Dir, FileName)
%
%   Dir      = Directory of data
%   FileName = Base file name (no extension)

A = load(fullfile(Dir, [FileName '.mat']));

VX = A.Velocity_X;
VY = A.Velocity_Y;

% Mean velocities (omit NaNs)
u = mean(VX, 2, 'omitnan');
v = mean(VY, 2, 'omitnan');
uv = [u(:) v(:)];
uv = uv(~any(isnan(uv),2), :); 

% PCA using covariance
C = cov(uv);
[V, D] = eig(C);
[eigVals, idx] = sort(diag(D), 'descend');
V = V(:, idx);

% Collect Major axis
majorAxis = V(:,1);

% Find rotation angle
theta = atan2(majorAxis(2), majorAxis(1));
theta = rad2deg(theta);     
fprintf('Rotation angle = %.2f degrees\n', theta);


% Ensure right-handed coordinate system
if det(V) < 0
    V(:,2) = -V(:,2);
end

% Rotate data
uv_rot = uv * V;

% Plot input vs rotated data
figure
plot(uv(:,1), uv(:,2), '.r')
hold on
scale = sqrt(eigVals);
plot([0 scale(1)*V(1,1)], [0 scale(1)*V(2,1)], '--b', ...
     [0 scale(2)*V(1,2)], [0 scale(2)*V(2,2)], '--b');
axis equal, grid on
xlabel('u'), ylabel('v')
title('Input Data', 'Eigenvectors')

figure
plot(uv_rot(:,1), uv_rot(:,2), '.g')
axis equal, grid on
xlabel('Major Axis'), ylabel('Minor Axis')
title('PCA Rotated Data')
legend('Rotated points')

% Apply rotation to full velocity fields
shape = size(VX);
Vrot = [VX(:) VY(:)] * V;

A.PCA_X = reshape(Vrot(:,1), shape);
A.PCA_Y = reshape(Vrot(:,2), shape);

disp('Saving PCA data...')
save(fullfile(Dir, [FileName '_PCA.mat']), '-struct', 'A');
end
