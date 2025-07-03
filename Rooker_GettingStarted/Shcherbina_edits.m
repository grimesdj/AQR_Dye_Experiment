 v_wrapped = Data.b1(dye, 1)';
 
 Vr = ((Data.sspeed(dye, :)'.^2)/(4*1000^2*2.5));

 
 
[nbins, nt] = size(v_wrapped);
v_vals = v_wrapped(:);  % Flatten all values (bin x time) into a column

% Create time index for each point (repeat 1:nt for each bin)
time_idx = repmat(1:nt, nbins, 1);
time_vals = time_idx(:);  % Flatten into column

v_unwrap = v_wrapped;

filt_diff = (v_wrapped - medfilt2(v_wrapped,[5,3]))./Vr;
filt_diff1 = filt_diff;
suspect_pts = abs(filt_diff)>0.5;


% Flatten suspect point mask
suspect_flags = suspect_pts(:);  % true for suspect, false for good

% Plot
figure
hold on
plot(time_vals(~suspect_flags), v_vals(~suspect_flags), 'g.')
plot(time_vals(suspect_flags), v_vals(suspect_flags), 'r.')
xlabel('Time step')
ylabel('Velocity')
title('All velocities (stacked bins), colored by suspicion')
