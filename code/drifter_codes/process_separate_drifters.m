[hour, minute, second, velocity, numsats, height, lat, long, quality, geoidsep, pdop, hdop, vdop] = readGPS('GPSUSER_133200089_20240417_184949.TXT');
t = datenum(2024*ones(size(hour)),4*ones(size(hour)),10*ones(size(hour)),hour,minute,second);
t = t(1:9:end);
lat = lat(1:9:end);
long = -long(1:9:end);

tstart = [datenum(2024,4,10,19,17,00) datenum(2024,4,10,19,28,00) datenum(2024,4,10,19,51,00)];
tstop = [datenum(2024,4,10,19,18,00) datenum(2024,4,10,19,30,00) datenum(2024,4,10,19,52,30)];
for i = 1:3
    
    ind = find(t>tstart(i) & t<tstop(i));
    figure
    geoscatter(lat(ind),long(ind),20)
    
    drifter2(i).t = t(ind);
    drifter2(i).lat = lat(ind);
    drifter2(i).lon = long(ind);
end


figure
for i = 1:3
    geoscatter(drifter1(i).lat,drifter1(i).lon,20); hold on
    geoscatter(drifter2(i).lat,drifter2(i).lon,20); 
end