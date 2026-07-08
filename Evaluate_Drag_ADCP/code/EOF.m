function [L, EOFs, EC, Error, Skill,lam, Barotropic] = EOF( Y , num, BT ); 
% EOF computes EOFs via SVD. 
%
% Usage:  [L, EOFs, EC, Error, Skill, lam, Barotropic] = EOF( Y, num ); 
%
% Y is the data matrix to be evaluated, and whose columns contain
% the different stations and rows contain timeseries from each
% station. num is the number of singular values to return, and used
% to compute ECs and Errors, etc.
%
% The mean of Y is removed prior to SVD. Only compute first EOFs
% corresponding to the number of stations (columns) in Y.
%
% BT is a binary flag that determines whether Barotropic Velocity gets removed prior to EOF analysis (defaults to off)
%
%
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

if nargin < 3 || isempty(BT)
    BT = 0;
end

[N, M] = size(Y);

% remove time mean
Y = Y - mean(Y, 1);

if BT
    % remove barotropic mode
    disp('rm barotropic...')
    o = ones(1, size(Y, 2));
    Ybar = mean(Y, 2);
    Barotropic = Ybar * o;
    Y = Y - Barotropic;
else
    Barotropic = zeros(size(Y));
end


% detrend
Y = detrend(Y);


if nargin <2 || isempty(num)
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