% converting Vector XYZ to ENU matching the orientation of Aquadopp

% get aquadopp orietation
%headAQD = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat", 'head');

% get Vector data
%Vector = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat");
function rotation(theta, phi, omega)
st = sind(theta);
sp = sind(phi);
so = sind(omega);
ct = cosd(theta);
cp = cosd(phi);
co = cosd(omega);
Head = [st ct 0; -(ct) st 0; 0 0 1]
P = [ 0 cp sp; -1 0 0; 0 -(sp) cp]
R = [0 1 0; -(co) 0 so; so 0 co]
H*P*R
end