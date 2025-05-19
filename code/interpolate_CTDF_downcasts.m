function [time_grid,latitude_grid,longitude_grid,pres_grid,dye_grid,temp_grid,cond_grid] = interpolate_CTDF_downcasts(time,pres,temp,cond,dye_raw,gps);
%
% USAGE: [time_grid,latitude_grid,longitude_grid,pres_grid,dye_grid,temp_grid,cond_grid] = interpolate_CTDF_downcasts(time,pres,temp,cond,dye_raw,gps);
%
%

% break record into individual casts
% 1) first smooth pressure record to find zero-up crossings
tau = 5;% smooth timescale in seconds
dt  = 1/8;% sample rate Hz
N   = tau/dt;
f   = hanning(N); f = f./sum(f);
pres_avg = conv(pres,f,'same');
dpres     = gradient(pres);
d2pres    = gradient(dpres);
dpres_avg = gradient(pres_avg);
%
%
% locate the 'casts'
cast          = pres_avg>=0.5;
downcast      = find(cast & dpres_avg>0);
upcast        = find(cast & dpres_avg<0);
cast          = find(cast);
%
[startINDs,endINDs] = Segment(cast,(1/dt));
cast_start     = cast(startINDs);
cast_end       = cast(endINDs);
valid          = (cast_end-cast_start)>10;
cast_start(~valid)=[];
cast_end(~valid)  =[];
Ncasts         = length(cast_start);
%
%
% $$$ figure, plot(datetime(time,'convertfrom','datenum'),pres_avg)
% $$$ hold on,xline(datetime(time(cast_start),'convertfrom','datenum'),'-b')
% $$$ hold on,xline(datetime(time(cast_end),'convertfrom','datenum'),'-r')
%
% now create a uniform pressure grid with 10cm resolution to the maximum depth (15 m)
pres_grid = [0:0.1:15]';
Npres     = length(pres_grid);
time_grid = nan(1,Ncasts);
latitude_grid  = nan(1,Ncasts);
longitude_grid = nan(1,Ncasts);
temp_grid = nan(Npres,Ncasts);
dye_grid  = nan(Npres,Ncasts);
cond_grid = nan(Npres,Ncasts);
for jj = 1:Ncasts
    this_cast = find(downcast>=cast_start(jj) & downcast<=cast_end(jj));
    [~,inds]  = unique(pres(downcast(this_cast)));
    time_grid(1,jj) = time(downcast(this_cast(inds(1))));
% $$$     latlon          = gps(1).latitude + sqrt(-1)*gps(1).longitude;
% $$$     latlon_grid     = interp1(gps(1).time,latlon,time_grid(1,jj),'linear','extrap');
    latitude_grid(1,jj)  = interp1(gps(1).time,gps(1).latitude,time_grid(1,jj));% real(latlon_grid);
    longitude_grid(1,jj) = interp1(gps(1).time,gps(1).longitude,time_grid(1,jj));% imag(latlon_grid);
    temp_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),temp(downcast(this_cast(inds))),pres_grid);
    try cond_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),cond(downcast(this_cast(inds))),pres_grid);
    catch
        cond_grid(:,jj) = nan;
    end
    try dye_grid(:,jj)  = interp1(pres(downcast(this_cast(inds))),dye_raw(downcast(this_cast(inds))),pres_grid);
    catch
        dye_grid(:,jj)  = nan;
    end
end
%
dye_grid( dye_grid<0 | dye_grid>max(dye_raw) ) = nan;
%
mask = cond_grid<30 | temp_grid>25 | temp_grid<5;
dye_grid(mask) = nan;
temp_grid(mask) = nan;
cond_grid(mask)= nan;