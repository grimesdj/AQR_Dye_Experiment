function pca_function(L0Dir, L0Name)
% first create random variable then rotate 30-degrees counter clockwise
A = load([L0Dir, '/', L0Name, '.mat']);

%VX = fillmissing(A.Velocity_X, 'nearest');
%VY = fillmissing(A.Velocity_Y, 'nearest');
VX = A.Velocity_X;
VY = A.Velocity_Y;

u = mean(VX, 2, 'omitnan');
v = mean(VY, 2, 'omitnan');

u(isnan(u)) = [];
v(isnan(v)) = [];
N  = length(u);


%
% rotate 30-deg
% uv = [cosd(30) -sind(30);...
%       sind(30)  cosd(30)]*[u0';v0'];
% u  = uv(1,:)';
% v  = uv(2,:)';
%
%figure,
p1 = plot(u,v,'.r');
%
%
% analytic PCA
u2avg = (u'*u)/(N-1);
v2avg = (v'*v)/(N-1);
uvavg = (u'*v)/(N-1);
sig2pos = 0.5*( (u2avg + v2avg) + sqrt( (u2avg-v2avg)^2 + 4*uvavg^2 ));
sig2neg = 0.5*( (u2avg + v2avg) - sqrt( (u2avg-v2avg)^2 + 4*uvavg^2 ));
%
% eigen vectors
%     [ v_positive  , v_negative    ];
V   = [ uvavg          ,  sig2neg-v2avg;...
        (sig2pos-u2avg),  uvavg];
%
% normalized eigen vectors
norm = diag(V'*V);
V    = V*diag(1./sqrt(norm));
%
% sort eigen values/vectors
eig  = [sig2pos,sig2neg];
[eigsort,srt] = sort(eig,'descend');
V    = V(:,srt);
%
% make sure it's right-handed using cross-product (determinate>0).
isRHcoord = det(V);
if (isRHcoord<0)
    if V(1)<0
        V(:,1)=-V(:,1);
    else
        V(:,2)=-V(:,2);
    end
end
%
% unrotate the data using principal vectors
uv1 = [u, v]*V;
u1  = uv1(:,1);
v1  = uv1(:,2);
% plot vectors (scaled by eigen-value) on figure
hold on,
plot( u1, v1, '.c')
p2 = plot( [0;sqrt(eig(1))*V(1,1)], [0;sqrt(eig(1))*V(2,1)],'--b', [0;sqrt(eig(2))*V(1,2)], [0;sqrt(eig(2))*V(2,2)],'--b');
%
axis equal
grid on
legend([p1(1) p2(1)],'input data','eigen vectors')
%
xlabel(' rotated-major-axis ')
ylabel(' rotated-minor-axis ')

% apply new rotation matrix individually 

% 
shape = size(VX);

VX1 = reshape(VX, [], 1);
VY1 = reshape(VY, [], 1);

Vrot = [VX1 VY1]*V;

A.PCA_X = reshape(Vrot(:,1), shape);
A.PCA_Y = reshape(Vrot(:,2), shape);
%
figure, plot( A.PCA_X, A.PCA_Y , '.c')
%p3 = plot( [0;sqrt(eig(1))*V(1,1)], [0;sqrt(eig(1))*V(2,1)],'--b', [0;sqrt(eig(2))*V(1,2)], [0;sqrt(eig(2))*V(2,2)],'--b');
%
axis equal
grid on
%
xlabel(' major-axis ')
ylabel(' minor-axis ')

disp('Saving PCA data')
save([L0Dir,'/',L0Name,'.mat'],'-struct','A')


end