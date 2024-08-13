function counts = despike_fluoro(raw_counts,dt);
%
% USAGE: counts = despike_fluoro(raw_counts,dt);
%
% raw_counts is a vector of instrument counts/raw ppb; dt is the sample interval in seconds. 
% function looks for isolated spikes of .02*dt*(1+dye/4) up or down, and
% any incidence of slopes >.05*dt*(1+dye/12) and nans them (threshold of
% approx. 1 counts and 3 counts per second, respectively at dye=0;
% and approx. 4 & 7 counts per second, respectively at dye=10;
% and approx. 30 & 40 counts per second, respectively at dye = 100). 

tmp = raw_counts;
[m,n] = size(tmp);

% $$$ % get order of dye
% $$$ ord = [.01 .1 1 10 25 50 75 100];
% $$$ ortmp2 = repmat(ord,m,1);
% $$$ [r,c] = size(ord);
% $$$ 
% $$$ tmp2 = repmat(tmp,1,c);
% $$$ [~,col_tmpord] = min(abs((tmp2./ortmp2-1)),[],2);
% $$$ scale = ord(col_tmpord(1:end-1))';

% remove isolated up spikes
delta0 = 0.2;
i1 = find(diff(tmp(2:end))./dt./(1+tmp(2:end-1))>delta0);% was /10>.02
i2 = find(diff(tmp(2:end))./dt./(1+tmp(1:end-2))<-delta0);
i3 = intersect(i1+1,i2)+1;

% remove isolated down spikes
i4 = find(diff(tmp(1:end-1))./dt./(1+tmp(1:end-2))<delta0);
i5 = find(diff(tmp(1:end-1))./dt./(1+tmp(2:end-1))>delta0);
i6 = intersect(i4+1,i5);

i = union(i3,i6);
tmp(i+1) = nan;

% remove residuals
delta1 = 0.5;
i1 = find(diff(tmp)./dt./(1+tmp(1:end-1))>delta1);% was /10>.05
i2 = find(diff(tmp)./dt./(1+tmp(1:end-1))<-delta1);
i = union(i1,i2);
tmp(i+1) = nan;

counts = tmp;