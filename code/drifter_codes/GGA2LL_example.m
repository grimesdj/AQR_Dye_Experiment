tmptable = readtable('GPSUSER_133201468_20220414_121839.TXT');

gt31 = GGA_2LL(tmptable);

figure;
geobasemap('satellite');
hold on;
geoscatter(gt31.lat,gt31.lon);