%
% 1) load M1-M4 temperature versus ... (mab?)
% 2) estimate dT/dz, then low-pass filter (10-hrs) to (T_lp, dTdz_lp)
% 3) estimate z' = (T-T_lp)/dTdz_lp;
% 4) plot timeseries of z' at depth with largest mean(dTdz_lp); or mid-water column
rootDir = '/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/2024_PROCESSED_DATA/';
Salinity = 33.6;
g = 9.81;
for jj = 1:4
    % 1) load...
    fin = sprintf([rootDir,'M%d/L1/mooring_M%d.mat'],jj,jj);
    tmp = load(fin);
    % use pressure to map from mab to depth relative to mean tide elevation
    if jj<3
        fin = sprintf([rootDir,'M%d/L1/ADCP/ADCP_M%d_L1.mat'],jj,jj);
        load(fin,'currents');
        p = mean(currents.Pressure);
        p = (p-tmp.Temperature_mab); 
    else
        p = (tmp.Pressure_mab+mean(tmp.Pressure,'omitnan')-tmp.Temperature_mab);
    end
    % 2) estimate density and dT/dz
    tmp.Density = gsw_rho(Salinity+0*tmp.Temperature,tmp.Temperature,p);%    eval(sprintf('M%d.Density = gsw_rho(Salinity,M%d.Temperature,p);',jj));
    tmp.dTdz    = diff(tmp.Temperature,1,1)./diff(tmp.Temperature_mab);
    tmp.N2      = -g/mean(tmp.Density,'all','omitnan').*( diff(tmp.Density,1,1)./diff(tmp.Temperature_mab)  );
    %
    %
    dt = (tmp.Time(2)-tmp.Time(1))*86400;
    Nt = 10*3600/dt;
    flt= hamming(Nt); flt = flt./sum(flt);
    T_lp = 
    %
    %
    eval(sprintf('M%d = tmp;',jj));
end
figure, hold on, tmp = M2; for jj=1:size(tmp.Temperature,1), plot(tmp.Time,tmp.Temperature(1:jj,:)','k',tmp.Time,tmp.Temperature(jj,:)','r'), set(gca,'xlim',tmp.Time(1)+[0 1]), pause(1),end



