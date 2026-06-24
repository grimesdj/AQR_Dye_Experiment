function [L, EOFs, EC, Error, Skill,lam] = EOF( Y , num ); 
% EOF computes EOFs via SVD. 
%
% Usage:  [L, EOFs, EC, Error] = EOF( Y, num ); 
%
% Y is the data matrix to be evaluated, and whose columns contain
% the different stations and rows contain timeseries from each
% station. num is the number of singular values to return, and used
% to compute ECs and Errors, etc.
%
% The mean of Y is removed prior to SVD. Only compute first EOFs
% corresponding to the number of stations (columns) in Y.
% 
% L: Eigenvalues 
%
% EOFs: The associated eigenvectors
%
% EC: The coefficients to reconstruct Y
%
% Errors: sum of squared residuals from EC*EOFs
%
% Skill: ratio of prediction variance to observed for each station
% (column)

[N, M] = size(Y);

% remove mean
Y = Y - repmat(mean(Y,1),N,1);

% detrend
Y = detrend(Y);

% remove barotropic mode
disp('rm barotropic...')
o = ones(1, size(Y, 2));
Ybar = mean(Y, 2);
Barotropic = Ybar * o;
Y = Y - Barotropic;

if nargin == 1
    [C, lam, EOFs] = svd(Y,0);% first M EOFs 
else
    [C, lam, EOFs] = svds(Y,num);% first num EOFs
end
 
EC = C*lam;% Expansion Coeffs.
L = diag(lam).^2./(N-1);% Eigenvalues

% COV = EC'*EC./(N-1); %Covariance matrix

diffs = (Y-EC*EOFs');
Error = sqrt(sum(diffs(:)'*diffs(:)));
Skill = var(EC*EOFs')./var(Y);