for i = 1:3
    u = sensorValues(dye, i);
% Calculate derivative surrogates
    for j = 1:length(u)
        k = i+1
        if i>1
        l = i-1
        else l = i
        end
        du(i) = (u(k)-u(l))/2;
        
    end
    for j = 1:length(u)
        k = i+1
        if i>1
        l = i-1
        else l = i
        end
        ddu(i) = (du(k)-du(l))/2;
    end
% Calculate std and absolute max

usig = std(u);
dusig = std(du);
ddusig = std(ddu);

n = length(data);

umax = sqrt(2*ln(n))*usig;
dumax = sqrt(2*ln(n))*dusig;
ddumax = sqrt(2*ln(n))*ddusig;

umin = -umax;
dumin = -dumax;
ddumin = -ddumax;

% Calculate roatation angle


% Calculate ellipse



theta = atan(sum(u .* ddu) / sum(u.^2));  

[u, v] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 50));
x = dumax * cos(u) .* sin(v);  % du
y = umax  * sin(u) .* sin(v);  % u
z = ddumax * cos(v);           % ddu

% Flatten points into 3xN
pts = [x(:)'; y(:)'; z(:)'];

% Rotation matrix around x-axis (du) by angle θ
Rx = [1      0           0;
      0  cos(theta) -sin(theta);
      0  sin(theta)  cos(theta)];

rotated_pts = Rx * pts;

% Reshape
x_rot = reshape(rotated_pts(1, :), size(x));
y_rot = reshape(rotated_pts(2, :), size(y));
z_rot = reshape(rotated_pts(3, :), size(z));

surf(x_rot, y_rot, z_rot, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
xlabel('du'), ylabel('u'), zlabel('ddu')
title('Tilted Ellipsoid in (du, u, ddu) Space')
axis equal
grid on
light, lighting gouraud, camlight headlight

end