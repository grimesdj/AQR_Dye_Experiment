% seperate CTD casts into dye releases
clear all
close all


load("../../../../Kelp_data/data/2024_PROCESSED_DATA/VesselCTDFData.mat")
load("../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat")


% NaN all dye readings that are more than 2 std from Release Temp
%%
[temp_pdf, temp_x] = ksdensity(R1.Temperature(:,3));
% Compute mean and standard deviation from the PDF
dx = temp_x(2) - temp_x(1);
Tbar = trapz(temp_x, temp_x .* temp_pdf);               % Mean
mu2 = trapz(temp_x, (temp_x.^2) .* temp_pdf);         % E[x^2]
sigma = sqrt(mu2 - Tbar^2);                             % Standard deviation
% Plot PDF
figure;
plot(temp_x, temp_pdf, 'b', 'LineWidth', 2)
xlabel('Temperature (°C)')
ylabel('Probability Density')
title('PDF of Release 1 Temperature')
grid on
% Add mean and ±1 std lines
hold on
xline(Tbar, 'g-', 'LineWidth', 1.5, 'Label', 'Mean')
xline(Tbar + sigma, 'r--', 'LineWidth', 1, 'Label', '+1σ')
xline(Tbar - sigma, 'r--', 'LineWidth', 1, 'Label', '-1σ')
hold off
% Optional: show text summary
fprintf('Mean Release temperature: %.2f °C\n', Tbar)
fprintf('Release Temp Standard deviation: %.2f °C\n', sigma)
%%
Tplus = Tbar+(3*sigma);
Tminus = Tbar-(3*sigma);

CTD1.dye_grid(find(CTD1.temp_grid>=Tplus | CTD1.temp_grid<=Tminus)) = NaN;

figure, imagesc(datenum(CTD1.time_grid), CTD1.pres_grid, log2(CTD1.dye_grid))
hold on, contour(datenum(CTD1.time_grid), CTD1.pres_grid, CTD1.temp_grid, [Tminus Tplus], 'LineWidth', 1.5, 'LineColor', 'k')
hold on, contour(datenum(CTD1.time_grid), CTD1.pres_grid, CTD1.temp_grid, [Tbar Tbar], 'LineWidth', 1.5, 'LineColor', 'r')
datetick('x', 'keeplimits')

jumps = find(abs(diff(CTD1.time_grid))>1e-3);

tdx = [];
sig = [];
mu  = [];
timestamp = [];

disp('finding transect indicies')
for i = 1:length(jumps)
    tdx = [tdx jumps(i) jumps(i)+1];
    
end
xline(datenum(CTD1.time_grid(tdx)))

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
    timestamp = [timestamp CTD1.time_grid(round(mean(tdx(i):tdx(i+1))))];

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
    %close(gcf)

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
    %close(gcf)
    end  
end


figure
    ax1 = subplot(2, 1, 1);
    ax2 = subplot(2, 1, 2);
    plot(ax1, mu)
    title(ax1, 'PDF AVG')
    plot(ax2, sig)
    title(ax2, 'PDF STD')
    xlabel('transect #')
    linkaxes([ax1 ax2], 'x')


sig_interp = interp1(timestamp, sig, CTD1.time_grid, 'linear', 'extrap');

p = polyfit(datenum(CTD1.time_grid), sig_interp, 1);
K = p(1)              
intercept = p(2);      
sig_fit = polyval([K intercept], datenum(CTD1.time_grid));
figure, plot(CTD1.time_grid, sig_interp, 'w')
hold on, plot(CTD1.time_grid, sig_fit, 'r--')
legend('Interpolated \sigma', 'Quadratic Fit', 'Location', 'best')

%CTD1.dye_trans1 = CTD1.dye_grid(:,find(CTD1.time_grid == CTD1.time_grid(transidx(1))):find(CTD1.time_grid==CTD1.time_grid(transidx(4))))