% seperate CTD casts into dye releases
clear all
close all


load("../../../../Kelp_data/data/2024_PROCESSED_DATA/VesselCTDFData.mat")

figure, imagesc(datetime(CTD1.time_grid, 'ConvertFrom', 'datenum'), CTD1.pres_grid, log2(CTD1.dye_grid))

jumps = find(abs(diff(CTD1.time_grid))>1e-3);

tdx = [];
sig = [];
mu  = [];

disp('finding transect indicies')
for i = 1:length(jumps)
    tdx = [tdx jumps(i) jumps(i)+1];
    
end
xline(datetime(CTD1.time_grid(tdx), 'ConvertFrom', 'datenum'))

disp('seperating casts into transects')

tdiff = diff(tdx);
count = 0;

for i = 1:length(tdiff)
    if tdiff(i)>1
    count = count+1;
    dyename = sprintf('dye_trans%d', count);
    tempname = sprintf('temp_trans%d', count);

    CTD1.(dyename)  = CTD1.dye_grid(:, tdx(i):tdx(i+1));
    CTD1.(tempname) = CTD1.temp_grid(:,tdx(i):tdx(i+1));

    % Clean up data
    valid = ~isnan(CTD1.(dyename)) & ~isnan(CTD1.(tempname));
    T = CTD1.(tempname)(valid);
    D = CTD1.(dyename)(valid);
    % Normalize dye values to make them sum to 1 (i.e., turn into weights)
    w = D / nansum(D);
    % Compute KDE of T, weighted by dye
    [f, x] = ksdensity(T, 'Weights', w);  % <-- this is key 
    % Plot
    
    figure, imagesc(log2(CTD1.(dyename)))
    title(sprintf('Transect #%d', count))
    
    figure, plot(x, f, 'LineWidth', 2)
    xlabel('Temperature (°C)')
    ylabel('Weighted PDF of Dye')
    title(sprintf('Transect #%d', count))


    fnorm = f / trapz(x,f);
    
    avg = trapz(x, x.*fnorm);
    mu = [mu avg];
    
    var = trapz(x, (x-avg).^2 .* fnorm);
    
    sig = [sig sqrt(var)];

    xline(mu(count), 'g')
    xline(mu(count)+sig(count), 'r')
    xline(mu(count)-sig(count), 'r')
    
    end  
end


figure
    ax1 = subplot(2, 1, 1);
    ax2 = subplot(2, 1, 2);
    plot(ax1, mu)
    title(ax1, 'AVG')
    plot(ax2, sig)
    title(ax2, 'STD')
    xlabel('transect #')
    linkaxes([ax1 ax2], 'x')

    K = mean(diff(sig))



%CTD1.dye_trans1 = CTD1.dye_grid(:,find(CTD1.time_grid == CTD1.time_grid(transidx(1))):find(CTD1.time_grid==CTD1.time_grid(transidx(4))))