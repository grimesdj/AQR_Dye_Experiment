function coord = GGA2LL(gps_table)
% Converts to Lat/Lon from GGA formatted coordinates
% coord = GGA_2LL(gps_table)
tmp_lat = gps_table(:,3);
tmp_lon = gps_table(:,5);
if istable(tmp_lat)
    tmp_lat = table2array(tmp_lat);
    tmp_lon = table2array(tmp_lon);    
end

lat_degree = floor(tmp_lat/100);
lat_decimal = round((tmp_lat-lat_degree*100)/60*10^8)/10^8; 
coord.lat = lat_degree+lat_decimal;

long_degree= floor(tmp_lon/100); 
long_decimal = round((tmp_lon-long_degree*100)/60*10^8)/10^8; 
coord.lon = long_degree+long_decimal;
coord.lon = -coord.lon;
end