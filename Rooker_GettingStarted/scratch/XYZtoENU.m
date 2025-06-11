% converting Vector XYZ to ENU matching the orientation of Aquadopp

% get aquadopp orietation
headAQD = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat", 'head');

% get Vector data
Vector = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_Vector_L0.mat");